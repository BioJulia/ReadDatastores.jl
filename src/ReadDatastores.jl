module ReadDatastores

export 
    PairedReadDatastore,
    PairedReadOrientation,
    FwRv,
    RvFw,
    LongReadDatastore,
    SequenceBuffer
    
using BioSequences, FASTX

const SeqDataStoreMAGIC = 0x05D5

@enum Filetype::UInt16 begin
    PairedDS = 1
    LongDS = 2
end

###
### Writing helpers
###

@inline function writestring(fd::IO, s::AbstractString)
    write(fd, s) + write(fd, '\0')
end

@inline function write_flat_vector(fd::IO, v::Vector{T}) where {T}
    n = UInt64(length(v))
    write(fd, n)
    write(fd, v)
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
    

include("paired-reads.jl")
include("long-reads.jl")
include("sequence-buffer.jl")

end # module
