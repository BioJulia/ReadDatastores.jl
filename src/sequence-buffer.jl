const DataStore = Union{PairedReadDatastore, LongReadDatastore}
const DEFAULT_BUF_SIZE = UInt64(1024 * 1024 * 30)
const DEFAULT_CHUNK_SIZE = UInt64(1024 * 1024 * 4)

mutable struct SequenceBuffer{T<:DataStore}
    data_store::T
    buffer::Vector{UInt8}
    chunk_size::Int
    buffer_position::Int
end

function SequenceBuffer(ds::T) where {T<:DataStore}
    return SequenceBuffer{T}(ds, Vector{UInt8}(undef, DEFAULT_BUF_SIZE), DEFAULT_CHUNK_SIZE, typemax(Int))
end

@inline chunksize(sb::SequenceBuffer) = sb.chunk_size
@inline buffersize(sb::SequenceBuffer) = length(buffer(sb))
@inline bufferpos(sb::SequenceBuffer) = sb.buffer_position
@inline bufferpos!(sb::SequenceBuffer, pos) = sb.buffer_position = pos
@inline datastore(sb::SequenceBuffer) = sb.data_store
@inline stream(sb::SequenceBuffer) = stream(datastore(sb))
@inline buffer(sb::SequenceBuffer) = sb.buffer

"""
    _load_sequence_data!(seq::LongSequence, buffer::Vector{UInt8}, offset::Integer)

This internal method reads the packed data of a sequence from a buffer of bytes.
"""
function _load_sequence_data!(seq::LongSequence, buffer::Vector{UInt8}, offset::Integer)
    seqdata = BioSequences.encoded_data(seq)
    for i in eachindex(seqdata)
        seqdata[i] = unsafe_load(convert(Ptr{UInt64}, pointer(buffer, offset + 1)))
        offset = offset + sizeof(UInt64)
    end
    return seq    
end

function _load_sequence(data::Vector{UInt8}, offset::Integer)
    sequence_length = unsafe_load(convert(Ptr{UInt64}, pointer(data, offset + 1)))
    offset = offset + sizeof(UInt64)
    seq = LongSequence{DNAAlphabet{4}}(sequence_length)
    return _load_sequence_data!(seq, data, offset)
end

Base.eltype(sb::SequenceBuffer) = LongSequence{DNAAlphabet{4}}
Base.firstindex(sb::SequenceBuffer) = 1
Base.lastindex(sb::SequenceBuffer) = lastindex(datastore(sb))
Base.length(sb::SequenceBuffer) = length(datastore(sb))
Base.eachindex(sb::SequenceBuffer) = Base.OneTo(lastindex(sb))
Base.IteratorSize(sb::SequenceBuffer) = Base.HasLength()
Base.IteratorEltype(sb::SequenceBuffer) = Base.HasEltype()

@inline function Base.getindex(sb::SequenceBuffer, idx::Integer)
    @boundscheck checkbounds(datastore(sb), idx)
    
    # Get the position in the io stream that the read data begins.
    roif = unsafe_read_offset_in_file(datastore(sb), idx)
    
    # Knowing the position in the iostream, check if the buffer will need
    # refreshing.
    if roif < bufferpos(sb) || (roif + chunksize(sb)) > bufferpos(sb) + buffersize(sb)
        bufferpos!(sb, roif)
        seek(stream(sb), roif)
        readbytes!(stream(sb), buffer(sb), buffersize(sb))
    end
    
    return _load_sequence(buffer(sb), roif - bufferpos(sb))
end

@inline function Base.iterate(sb::SequenceBuffer, state = 1)
    @inbounds if firstindex(sb) ≤ state ≤ lastindex(sb)
        return sb[state], state + 1
    else
        return nothing
    end
end