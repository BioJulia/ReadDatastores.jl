const DEFAULT_BUF_SIZE = UInt64(1024 * 1024 * 30)
const DEFAULT_CHUNK_SIZE = UInt64(1024 * 1024 * 4)

mutable struct SequenceBuffer{T<:ReadDatastore}
    data_store::T
    bufferdata::Vector{UInt8}
    chunk_size::Int
    buffer_position::Int
end

function SequenceBuffer(ds::T) where {T<:ReadDatastore}
    return SequenceBuffer{T}(ds, Vector{UInt8}(undef, DEFAULT_BUF_SIZE), DEFAULT_CHUNK_SIZE, typemax(Int))
end

@inline chunksize(sb::SequenceBuffer) = sb.chunk_size
@inline buffersize(sb::SequenceBuffer) = length(bufferdata(sb))
@inline bufferpos(sb::SequenceBuffer) = sb.buffer_position
@inline bufferpos!(sb::SequenceBuffer, pos) = sb.buffer_position = pos
@inline datastore(sb::SequenceBuffer) = sb.data_store
@inline stream(sb::SequenceBuffer) = stream(datastore(sb))
@inline bufferdata(sb::SequenceBuffer) = sb.bufferdata

#=
function _load_sequence(data::Vector{UInt8}, offset::Integer)
    sequence_length = unsafe_load(convert(Ptr{UInt64}, pointer(data, offset + 1)))
    offset = offset + sizeof(UInt64)
    seq = LongSequence{DNAAlphabet{4}}(sequence_length)
    return _load_sequence_data!(seq, data, offset)
end
=#

Base.eltype(sb::SequenceBuffer) = Base.eltype(datastore(sb))
Base.firstindex(sb::SequenceBuffer) = firstindex(datastore(sb))
Base.lastindex(sb::SequenceBuffer) = lastindex(datastore(sb))
Base.length(sb::SequenceBuffer) = length(datastore(sb))
Base.eachindex(sb::SequenceBuffer) = eachindex(datastore(sb))
Base.IteratorSize(sb::SequenceBuffer) = Base.IteratorSize(datastore(sb))
Base.IteratorEltype(sb::SequenceBuffer) = Base.IteratorEltype(datastore(sb))

Base.checkbounds(sb::SequenceBuffer, i::Integer) = Base.checkbounds(datastore(sb), i)

function _check_for_buffer_refresh(sb::SequenceBuffer, filepos::Integer)
    if filepos < bufferpos(sb) || (filepos + chunksize(sb)) > bufferpos(sb) + buffersize(sb)
        bufferpos!(sb, filepos)
        seek(stream(sb), filepos)
        readbytes!(stream(sb), bufferdata(sb), buffersize(sb))
    end
end

function _check_for_buffer_refresh(sb::SequenceBuffer{<:LongReads}, filepos::ReadPosSize)
    # The default chunk size is 4 megs, which should be more than enough for long reads,
    # but we place this guard here just in case.
    if chunksize(sb) < filepos.sequence_size
        message = string("Reading a sequence from buffered datastore ",
                         name(datastore(sb)),
                         " failed!\nThe size of the buffer chunk is smaller than a sequence.\nIncrease the buffer chunk size to fix this.")
        throw(ErrorException(message))
    end
    _check_for_buffer_refresh(sb, filepos.offset) 
end

function _load_sequence_data!(seq::LongSequence, sb::SequenceBuffer, offset::Integer)
    bufdata = bufferdata(sb)
    seqdata = BioSequences.encoded_data(seq)
    GC.@preserve bufdata begin
        for i in eachindex(seqdata)
            seqdata[i] = unsafe_load(convert(Ptr{UInt64}, pointer(bufdata, offset + 1)))
            offset = offset + sizeof(UInt64)
        end
    end
    return seq    
end

function _load_sequence_from_pos(sb::SequenceBuffer, pos::Integer)
    offset = pos - bufferpos(sb)
    sequence_length = unsafe_load(convert(Ptr{UInt64}, pointer(bufferdata(sb), offset + 1)))
    offset = offset + sizeof(UInt64)
    seq = LongSequence{DNAAlphabet{4}}(sequence_length)
    return _load_sequence_data!(seq, sb, offset)
end

function _load_sequence_from_pos!(sb::SequenceBuffer, pos::Integer, seq::LongSequence)
    offset = pos - bufferpos(sb)
    sequence_length = unsafe_load(convert(Ptr{UInt64}, pointer(bufferdata(sb), offset + 1)))
    offset = offset + sizeof(UInt64)
    resize!(seq, sequence_length)
    return _load_sequence_data!(seq, sb, offset)
end

function _load_sequence_from_pos(sb::SequenceBuffer, pos::ReadPosSize)
    offset = pos.offset - bufferpos(sb)
    seq = LongSequence{DNAAlphabet{4}}(pos.sequence_size)
    return _load_sequence_data!(seq, sb, offset)
end

function _load_sequence_from_pos!(sb::SequenceBuffer, pos::ReadPosSize, seq::LongSequence)
    offset = pos.offset - bufferpos(sb)
    resize!(seq, pos.sequence_size)
    return _load_sequence_data!(seq, sb, offset)
end


@inline function Base.getindex(sb::SequenceBuffer, idx::Integer)
    @boundscheck checkbounds(sb, idx)
    # Get the position in the file that the read data begins.
    position_in_file = _inbounds_index_of_sequence(datastore(sb), idx)
    # Knowing the position in the file, check if the buffer will need
    # refreshing. To avoid refresh, the buffer needs to fully contain the read
    # desired.
    _check_for_buffer_refresh(sb, position_in_file)
    # Now we know the buffer was refreshed if needed, get the sequence from the
    # buffer. This function will translate `position_in_file` to an offset in
    # the buffer, and load the sequence.
    return _load_sequence_from_pos(sb, position_in_file)
end

@inline function inbounds_load_sequence!(sb::SequenceBuffer{DS}, i::Integer, seq::T) where {T<:LongSequence,DS<:ReadDatastore{T}}
    # Get the position in the file that the read data begins.
    position_in_file = _inbounds_index_of_sequence(datastore(sb), i)
    # Knowing the position in the file, check if the buffer will need
    # refreshing. To avoid refresh, the buffer needs to fully contain the read
    # desired.
    _check_for_buffer_refresh(sb, position_in_file)
    # Now we know the buffer was refreshed if needed, get the sequence from the
    # buffer. This function will translate `position_in_file` to an offset in
    # the buffer, and load the sequence.
    return _load_sequence_from_pos!(sb, position_in_file, seq)
end

@inline function load_sequence!(sb::SequenceBuffer{DS}, i::Integer, seq::T) where {T<:LongSequence,DS<:ReadDatastore{T}}
    checkbounds(sb, i)
    return inbounds_load_sequence!(sb, i, seq)
end

@inline function Base.iterate(sb::SequenceBuffer, state = 1)
    @inbounds if firstindex(sb) ≤ state ≤ lastindex(sb)
        return sb[state], state + 1
    else
        return nothing
    end
end