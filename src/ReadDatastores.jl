module ReadDatastores

export
    ReadDatastore,
    PairedReads,
    PairedReadOrientation,
    FwRv,
    RvFw,
    LongReads,
    LinkedReads,
    UCDavis10x,
    SequenceBuffer,
    
    name,
    maxseqlen,
    orientation,
    load_sequence!,
    buffer,
    read_tag
    
    
using BioSequences, FASTX

const ReadDatastoreMAGIC = 0x05D5

@enum Filetype::UInt16 begin
    PairedDS = 1
    LongDS = 2
    LinkedDS = 3
end

###
### Writing helpers
###

@inline function writestring(fd::IO, s::AbstractString)
    write(fd, s) + write(fd, '\0')
end

@inline function write_flat_vector(fd::IO, v::Vector{T}) where {T}
    n = UInt64(length(v))
    return write(fd, n) + write(fd, v)
end

@inline function read_flat_vector!(fd::IO, v::Vector{T}) where {T}
    n = read(fd, UInt64)
    v = resize!(n)
    unsafe_read(fd, pointer(v), n * sizeof(T))
    return v
end

@inline function read_flat_vector(fd::IO, ::Type{T}) where {T}
    n = read(fd, UInt64)
    v = Vector{T}(undef, n)
    unsafe_read(fd, pointer(v), n * sizeof(T))
    return v
end
    
###
### Abstract ReadDatastore type
###
abstract type ReadDatastore{S<:BioSequence} end

###
### Concrete ReadDatastore type definitions
###
include("paired-reads.jl")
include("long-reads.jl")
include("linked-reads.jl")
include("sequence-buffer.jl")

###
### ReadDatastore generics
###
@inline buffer(ds::ReadDatastore) = SequenceBuffer(ds)

"Get the length of the longest sequence in the datastore"
@inline maxseqlen(ds::ReadDatastore) = ds.max_read_len

"Get the name of the datastore"
@inline name(ds::ReadDatastore) = ds.name

"Get a reference to the underlying filestream the datastore is using"
@inline stream(ds::ReadDatastore) = ds.stream

Base.firstindex(ds::ReadDatastore) = 1
Base.lastindex(ds::ReadDatastore) = length(ds)
Base.eachindex(ds::ReadDatastore) = Base.OneTo(lastindex(ds))
Base.IteratorSize(ds::ReadDatastore) = Base.HasLength()
Base.IteratorEltype(ds::ReadDatastore) = Base.HasEltype()
Base.eltype(ds::ReadDatastore{T}) where {T} = T
Base.close(ds::ReadDatastore) = close(stream(ds))

@inline function Base.checkbounds(ds::ReadDatastore, i::Integer)
    if firstindex(ds) ≤ i ≤ lastindex(ds)
        return true
    end
    throw(BoundsError(ds, i))
end

@inline function Base.iterate(ds::ReadDatastore, state = 1)
    @inbounds if firstindex(ds) ≤ state ≤ lastindex(ds)
        return ds[state], state + 1
    else
        return nothing
    end
end

function Base.open(f::Function, ::Type{T}, filepath::AbstractString) where {T<:ReadDatastore}
    ds = open(T, filepath)
    try
        f(ds)
    finally
        close(ds)
    end
end

## Loading a sequence from a datastore - generic internals

function _inbounds_index_of_sequence end

"""
    _load_sequence_data!(prds::PairedReads, seq::LongSequence{DNAAlphabet{4}})

Reads the sequence data from the current position of a datastore's stream, into
a destination sequence `seq`.

This function knows how much data to read from the stream, based on the sequence's
data blob (i.e. `sizeof(BioSequences.encoded_data(seq))`).

!!! warning
    Uses `unsafe_read`, makes the following assumptions - does not check any of them.
    
    - The stream is at the correct position, such that an `unsafe_read` will
      read the entirety of the sequence before hitting the end of file.
      Hitting EOF should result in an error being thrown anyway so...
    - The destination sequence is currently the correct size and has been resized
      appropriately before this function has been called.
    - The data in the stream of the datastore is encoded appropriately for the
      sequence type.
"""
function _load_sequence_data!(ds::ReadDatastore{T}, seq::T) where {T<:LongSequence}
    seqdata = BioSequences.encoded_data(seq)
    GC.@preserve seqdata unsafe_read(stream(ds), pointer(seqdata), sizeof(seqdata))
    return seq
end

function _load_sequence_from_file_pos!(ds::ReadDatastore{T}, pos::Integer, seq::T) where {T<:LongSequence}
    seek(stream(ds), pos)
    seqlen = read(stream(ds), UInt64)
    resize!(seq, seqlen)
    return _load_sequence_data!(ds, seq)
end

function _load_sequence_from_file_pos!(ds::ReadDatastore{T}, pos::ReadPosSize, seq::T) where {T<:LongSequence}
    seek(stream(ds), pos.offset)
    resize!(seq, pos.sequence_size)
    return _load_sequence_data!(ds, seq)
end
    

@inline function inbounds_load_sequence!(ds::ReadDatastore{T}, i::Integer, seq::T) where {T<:LongSequence}
    pos = _inbounds_index_of_sequence(ds, i)
    return _load_sequence_from_file_pos!(ds, pos, seq)
end

@inline function load_sequence!(ds::ReadDatastore{T}, i::Integer, seq::T) where {T<:LongSequence}
    checkbounds(ds, i)
    return inbounds_load_sequence!(ds, i, seq)
end

end # module
