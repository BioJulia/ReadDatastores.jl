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
    DatastoreBuffer,
    name,
    max_read_length,
    orientation,
    load_sequence!,
    buffer,
    read_tag,
    stream,
    @openreads,
    @reads_str

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
abstract type ReadDatastore{S<:LongSequence} end

###
### ReadDatastore custom exceptions
###
struct MissingMagicError <: Exception
    filename::String
end

function Base.showerror(io::IO, e::MissingMagicError)
    print(io, "MissingMagicError: the file ")
    print(io, e.filename)
    print(io, " does not appear to be a valid read datastore file")
    print(io, ", it does not begin with the expected magic bytes.") 
end

struct DatastoreTypeError{D<:ReadDatastore} <: Exception
    filename::String
    dstype::Filetype
end

function Base.showerror(io::IO, e::DatastoreTypeError{D}) where {D}
    print(io, "DatastoreTypeError: ")
    print(io, e.filename)
    print(io, " contains a ")
    if e.dstype === PairedDS
        print(io, "paired read datastore")
    elseif e.dstype === LongDS
        print(io, "long read datastore")
    elseif e.dstype === LinkedDS
        print(io, "linked read datastore")
    else
        print(io, "unknown datastore type")
    end
    print(io, " and cannot be opened as a ")
    print(io, D)
end

struct DatastoreVersionError{D<:ReadDatastore} <: Exception
    version::Int
end

struct DatastoreEncodingError{D<:ReadDatastore} <: Exception
    filename::String
    bpe::Int
end

function Base.showerror(io::IO, e::DatastoreEncodingError{D}) where {D}
    print(io, "DatastoreEncodingError: ")
    print(io, e.filename)
    print(io, " encodes reads using ")
    print(io, e.bpe)
    print(io, " bits per element and cannot be opened as a ")
    print(io, D)
end

function __validate_datastore_header(io::IOStream, dstype::Type{D}) where {D<:ReadDatastore}
    magic = read(io, UInt16)
    dstype = reinterpret(Filetype, read(io, UInt16))
    version = read(io, UInt16)
    
    if magic !== ReadDatastoreMAGIC
        throw(MissingMagicError(io.name))
    end

    if dstype !== __ds_type_code(D)
        throw(DatastoreTypeError{D}(io.name, dstype))
    end

    if version !== __ds_version_code(D)
        throw(DatastoreVersionError{D}(version))
    end

    bps = read(io, UInt64)

    if bps != BioSequences.bits_per_symbol(Alphabet(eltype(D)))
        throw(DatastoreEncodingError{D}(io.name, bps))
    end
end

###
### Concrete ReadDatastore type definitions
###
include("short-reads.jl")
include("paired-reads.jl")
include("long-reads.jl")
include("linked-reads.jl")
include("sequence-buffer.jl")

function Base.showerror(io::IO, e::DatastoreVersionError{<:PairedReads})
    print(io, "DatastoreVersionError: ")
    print(io, "file format version of paired read datastore file (v")
    print(io, Int(e.version))
    print(io, ") is deprecated: this version of ReadDatastores.jl supports v")
    print(io, Int(PairedDS_Version))
end

function Base.showerror(io::IO, e::DatastoreVersionError{<:LongReads})
    print(io, "DatastoreVersionError: ")
    print(io, "file format version of long read datastore file (v")
    print(io, Int(e.version))
    print(io, ") is deprecated: this version of ReadDatastores.jl supports v")
    print(io, Int(LongDS_Version))
end

function Base.showerror(io::IO, e::DatastoreVersionError{<:LinkedReads})
    print(io, "DatastoreVersionError: ")
    print(io, "file format version of linked read datastore file (v")
    print(io, Int(e.version))
    print(io, ") is deprecated: this version of ReadDatastores.jl supports v")
    print(io, Int(LinkedDS_Version))
end

###
### ReadDatastore generics
###
@inline buffer(ds::ReadDatastore, buffer_size = DEFAULT_BUF_SIZE) = DatastoreBuffer(ds, buffer_size)

# TODO: Phase out
#"Get the length of the longest sequence in the datastore"
#@inline maxseqlen(ds::ReadDatastore) = ds.max_read_len

"Get the name of the datastore"
@inline name(ds::ReadDatastore) = ds.name

"Get a reference to the underlying filestream the datastore is using"
@inline stream(ds::ReadDatastore) = ds.stream

@inline Base.firstindex(ds::ReadDatastore) = 1
@inline Base.lastindex(ds::ReadDatastore) = length(ds)
@inline Base.eachindex(ds::ReadDatastore) = Base.OneTo(lastindex(ds))
@inline Base.IteratorSize(ds::ReadDatastore) = Base.HasLength()
@inline Base.IteratorEltype(ds::ReadDatastore) = Base.HasEltype()
@inline Base.eltype(::Type{<:ReadDatastore{T}}) where {T} = T
@inline Base.close(ds::ReadDatastore) = close(stream(ds))

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

function deduce_datastore_type(filename::String)::DataType
    open(filename, "r") do io
        seekstart(io)
        mn = read(io, UInt16)
        @assert mn === ReadDatastoreMAGIC
        tp = reinterpret(Filetype, read(io, UInt16))
        vn = read(io, UInt16)
        bpn = Int64(read(io, UInt64))
        if tp === PairedDS
            #@assert vn === PairedDS_Version
            if vn !== PairedDS_Version
                throw(DatastoreVersionError{PairedReads{DNAAlphabet{bpn}}}(vn))
            end
            out = PairedReads{DNAAlphabet{bpn}}
        elseif tp === LongDS
            @assert vn === LongDS_Version
            out = LongReads{DNAAlphabet{bpn}}
        elseif tp === LinkedDS
            @assert vn === LinkedDS_Version
            out = LinkedReads{DNAAlphabet{bpn}}
        else
            error("Unrecognized datastore type")
        end
        return out
    end
end

"""
    @openreads filename::String

A convenience macro for opening read datastores.

The open method for ReadDatastores requires the specific type of the ReadDatastore
as the first argument.

This can be cumbersome as the specific type of the datastore may be forgotten or
not known in the first place if the user recieved the datastore from a colleague.

This macro allows you to open a datastore using only the filename.

The macro opens the file, peeks at the header and deduces the type of datastore
contained in the file. It then returns a complete and correctly formed `open`
command for the datastore. This allows the user to forget about the specific
datastore type whilst still maintaining type certainty.
"""
macro openreads(filename::String)
    dstype = deduce_datastore_type(filename)
    return :(open($dstype, $filename))
end

"""
@openreads filename::String name

A convenience macro for opening read datastores and assigning it a name other
than the datastore's default name.

The open method for ReadDatastores requires the specific type of the ReadDatastore
as the first argument.

This can be cumbersome as the specific type of the datastore may be forgotten or
not known in the first place if the user recieved the datastore from a colleague.

This macro allows you to open a datastore using only the filename.

The macro opens the file, peeks at the header and deduces the type of datastore
contained in the file. It then returns a complete and correctly formed `open`
command for the datastore. This allows the user to forget about the specific
datastore type whilst still maintaining type certainty.
"""
macro openreads(filename::String, name::Symbol)
    dstype = deduce_datastore_type(filename)
    return :(open($dstype, $filename, $(QuoteNode(name))))
end

macro openreads(filename::String, binding::Symbol, expression::Expr)
    dstype = deduce_datastore_type(filename)
    quote
        open($dstype, $filename) do $(esc(binding))
            $(esc(expression))
        end
    end
end

macro openreads(filename::String, name::Symbol, binding::Symbol, expression::Expr)
    dstype = deduce_datastore_type(filename)
    quote
        open($dstype, $filename, $(QuoteNode(name))) do $(esc(binding))
            $(esc(expression))
        end
    end
end

macro openreads(filename::String, name::Symbol, fn::Symbol)
    dstype = deduce_datastore_type(filename)
    return :(open($fn, $dstype, $filename, $(QuoteNode(name))))
end

macro reads_str(filename::String)
    dstype = deduce_datastore_type(filename)
    return open(dstype, filename)
end

macro reads_str(filename::String, flag)
    dstype = deduce_datastore_type(filename)
    if flag === "s"
        return open(dstype, filename)
    elseif flag === "d"
        return :(open($dstype, $filename))
    else
        error("Invalid flag option")
    end
end

end # module
