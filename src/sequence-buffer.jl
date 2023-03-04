const DEFAULT_BUF_SIZE = UInt64(1024 * 1024 * 30)

@inline megabytes(n::T) where {T<:Unsigned} = (T(1024) * T(1024)) * n
@inline megabytes(n::Signed) = megabytes(unsigned(n))

@inline gigabytes(n::T) where {T<:Unsigned} = (T(1024) * T(1024) * T(1024)) * n
@inline gigabytes(n::Signed) = gigabytes(unsigned(n))

mutable struct DatastoreBuffer{T<:ReadDatastore}
    data_store::T
    bufferdata::Vector{UInt8}
    buffer_position::Int
end

function DatastoreBuffer(ds::T, buffer_size = DEFAULT_BUF_SIZE) where {T<:ReadDatastore}
    return DatastoreBuffer{T}(ds, Vector{UInt8}(undef, buffer_size), typemax(Int))
end

function DatastoreBuffer(ds::ShortReads, buffer_size = DEFAULT_BUF_SIZE)
    if _bytes_per_read(ds) > buffer_size
        error("Desired buffer size is too small")
    end
    return DatastoreBuffer{typeof(ds)}(ds, Vector{UInt8}(undef, buffer_size), typemax(Int))
end

@inline buffer_len(sb::DatastoreBuffer) = length(buffer_array(sb))
@inline buffer_position(sb::DatastoreBuffer) = sb.buffer_position
@inline bufferpos!(sb::DatastoreBuffer, pos) = sb.buffer_position = pos
@inline datastore(sb::DatastoreBuffer) = sb.data_store
@inline stream(sb::DatastoreBuffer) = stream(datastore(sb))
@inline buffer_array(sb::DatastoreBuffer) = sb.bufferdata

@inline Base.eltype(sb::DatastoreBuffer) = Base.eltype(datastore(sb))
@inline Base.firstindex(sb::DatastoreBuffer) = firstindex(datastore(sb))
@inline Base.lastindex(sb::DatastoreBuffer) = lastindex(datastore(sb))
@inline Base.length(sb::DatastoreBuffer) = length(datastore(sb))
@inline Base.eachindex(sb::DatastoreBuffer) = eachindex(datastore(sb))
@inline Base.IteratorSize(sb::DatastoreBuffer) = Base.IteratorSize(datastore(sb))
@inline Base.IteratorEltype(sb::DatastoreBuffer) = Base.IteratorEltype(datastore(sb))
@inline Base.checkbounds(sb::DatastoreBuffer, i) = Base.checkbounds(datastore(sb), i)

@inline function _load_sequence_data!(seq::LongSequence{A}, sb::DatastoreBuffer, offset::Integer) where {A<:DNAAlphabet}
    bufdata = buffer_array(sb)
    seqdata = seq.data
    GC.@preserve bufdata begin
        for i in eachindex(seqdata)
            seqdata[i] = unsafe_load(convert(Ptr{UInt64}, pointer(bufdata, offset + 1)))
            offset = offset + sizeof(UInt64)
        end
    end
    return seq    
end

# Short Reads specific buffering.

@inline function _check_for_buffer_refresh!(sb::DatastoreBuffer{<:ShortReads{<:DNAAlphabet}}, file_offset::Integer)
    # IF the desired data is contained in the buffer, there's no need to refresh.
    # Otherwise, there is!
    buf_start = buffer_position(sb)
    buf_size = buffer_len(sb)
    file_end = file_offset + _bytes_per_read(datastore(sb))
    if file_offset < buf_start || file_end > buf_start + buf_size
        bufferpos!(sb, file_offset)
        seek(stream(sb), file_offset)
        readbytes!(stream(sb), buffer_array(sb), buf_size)
    end
end

@inline function Base.getindex(sb::DatastoreBuffer{<:ShortReads{<:DNAAlphabet}}, idx::Integer)
    @boundscheck checkbounds(sb, idx)
    file_offset = _offset_of_sequence(datastore(sb), idx)
    _check_for_buffer_refresh!(sb, file_offset)
    
    buffer_offset = file_offset - buffer_position(sb)
    sequence_length = unsafe_load(convert(Ptr{UInt64}, pointer(buffer_array(sb), buffer_offset + 1)))
    buffer_offset = buffer_offset + sizeof(UInt64)
    seq = eltype(sb)(undef, sequence_length)
    
    return _load_sequence_data!(seq, sb, buffer_offset)
end

@inline function inbounds_load_sequence!(sb::DatastoreBuffer{<:ShortReads{A}}, i::Integer, seq::LongSequence{A}) where {A<:DNAAlphabet}
    file_offset = _offset_of_sequence(datastore(sb), i)
    _check_for_buffer_refresh!(sb, file_offset)
    
    buffer_offset = file_offset - buffer_position(sb)
    sequence_length = unsafe_load(convert(Ptr{UInt64}, pointer(buffer_array(sb), buffer_offset + 1)))
    buffer_offset = buffer_offset + sizeof(UInt64)
    resize!(seq, sequence_length)
    
    return _load_sequence_data!(seq, sb, buffer_offset)
end

# Long reads specific buffering

@inline function _check_for_buffer_refresh!(sb::DatastoreBuffer{LongReads{A}}, filepos::ReadPosSize) where {A<:DNAAlphabet}
    # If the buffer is too small to even fit the sequence. It will need to be
    # resized to be made bigger.
    n_bytes_req = cld(filepos.sequence_size, div(8, BioSequences.bits_per_symbol(A())))
    buf_start = buffer_position(sb)
    buf_size = buffer_len(sb)
    if buf_size < n_bytes_req
        @info "Resizing buffer!"
        resize!(sb.bufferdata, n_bytes_req)
        buf_size = n_bytes_req
    end
    file_start = filepos.offset
    file_end = file_start + n_bytes_req
    if file_start < buf_start || file_end > buf_start + buf_size
        bufferpos!(sb, file_start)
        seek(stream(sb), file_start)
        readbytes!(stream(sb), buffer_array(sb), buf_size)
    end
end

@inline function inbounds_load_sequence!(sb::DatastoreBuffer{LongReads{A}}, i::Integer, seq::LongSequence{A}) where {A<:DNAAlphabet}
    file_index = _inbounds_index_of_sequence(datastore(sb), i)
    _check_for_buffer_refresh!(sb, file_index)
    resize!(seq, file_index.sequence_size)
    buffer_offset = file_index.offset - buffer_position(sb)
    return _load_sequence_data!(seq, sb, buffer_offset)
end

@inline function Base.getindex(sb::DatastoreBuffer{LongReads{A}}, idx::Integer) where {A<:DNAAlphabet}
    @boundscheck checkbounds(sb, idx)
    file_index = _inbounds_index_of_sequence(datastore(sb), idx)
    _check_for_buffer_refresh!(sb, file_index)
    seq = eltype(sb)(undef, file_index.sequence_size)
    buffer_offset = file_index.offset - buffer_position(sb)
    return _load_sequence_data!(seq, sb, buffer_offset)
end

@inline function load_sequence!(sb::DatastoreBuffer{DS}, i::Integer, seq::T) where {T<:LongSequence,DS<:ReadDatastore{T}}
    checkbounds(sb, i)
    return inbounds_load_sequence!(sb, i, seq)
end

@inline function Base.iterate(sb::DatastoreBuffer, state = 1)
    @inbounds if firstindex(sb) ≤ state ≤ lastindex(sb)
        return sb[state], state + 1
    else
        return nothing
    end
end

@inline Base.summary(io::IO, sb::DatastoreBuffer) = print(io, "Buffered ", summary(datastore(sb)))

@inline Base.show(io::IO, sb::DatastoreBuffer) = summary(io, sb)