struct ReadPosSize
    offset::UInt64
    sequence_size::UInt64
end

Base.:(==)(x::ReadPosSize, y::ReadPosSize) = x.offset == y.offset && x.sequence_size == y.sequence_size

struct LongReads{A<:DNAAlphabet} <: ReadDatastore{LongSequence{A}}
    filename::String
    name::Symbol
    default_name::Symbol
    read_to_file_positions::Vector{ReadPosSize}
    stream::IO
end

index(lrds::LongReads) = lrds.read_to_file_positions

"Get the length of the longest sequence in the datastore"
function max_read_length(lrds::LongReads)
    max = zero(UInt64)
    @inbounds for i in index(lrds)
        ss = i.sequence_size
        if i.sequence_size > max
            max = i.sequence_size
        end
    end
    return max
end

const LongDS_Version = 0x0003

__ds_type_code(::Type{<:LongReads})   = LongDS
__ds_version_code(::Type{<:LongReads}) = LongDS_Version

###
### LongReads Header
###

# | Field                         | Value  | Type        |
# |:-----------------------------:|:------:|:-----------:|
# | Magic number                  | 0x05D5 | UInt16      | 2
# | Datastore type                | 0x0002 | UInt16      | 2
# | Version number                | 0x0001 | UInt16      | 2
# | Number of bits used per nuc   | 2 or 4 | UInt64      | 8
# | Index position in file        | N/A    | UInt64      | 8
# | Default name of the datastore | N/A    | String      | N


"""
    LongReads{A}(rdr::FASTQ.Reader, outfile::String, name::String, min_size::Integer) where {A<:DNAAlphabet}

Construct a Long Read Datastore from a FASTQ file reader.

# Arguments
- `rdr::FASTQ.Reader`: The reader of the fastq formatted file.
- `outfile::String`: A prefix for the datastore's filename, the full filename
  will include a ".loseq" extension, which will be added automatically.
- `name::String`: A string denoting a default name for your datastore.
  Naming datastores is useful for downstream applications.
- `min_size::Integer`: A minimum read length (in base pairs). When building
  the datastore, if any read sequence is shorter than this cutoff, then the read
  will be discarded.

# Examples
```jldoctest
julia> using FASTX, ReadDatastores

julia> longrdr = open(FASTQ.Reader, "test/human_nanopore_tester_2D.fastq")
FASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)

julia> ds = LongReads{DNAAlphabet{2}}(longrdr, "human-nanopore-tester", "nanopore-test", 0)
[ Info: Building long read datastore from FASTQ file
[ Info: Writing long reads to datastore
[ Info: Done writing paired read sequences to datastore
[ Info: 0 reads were discarded due to a too short sequence
[ Info: Writing index to datastore
[ Info: Built long read datastore with 10 reads
Long Read Datastore 'nanopore-test': 10 reads

```
"""
LongReads{A}(rdr::FASTQ.Reader, outfile::String, name::Union{String,Symbol}, min_size::Integer) where {A<:DNAAlphabet} = LongReads{A}(rdr, outfile, name, convert(UInt64, min_size))

function LongReads{A}(rdr::FASTQ.Reader, outfile::String, name::Union{String,Symbol}, min_size::UInt64) where {A<:DNAAlphabet}
    discarded = 0
    
    read_to_file_position = Vector{ReadPosSize}()
    ofs = open(outfile * ".loseq", "w")
    
    bps = UInt64(BioSequences.bits_per_symbol(A()))
    
    write(ofs, ReadDatastoreMAGIC, LongDS, LongDS_Version, bps, zero(UInt64))
    
    writestring(ofs, String(name))
    
    record = FASTQ.Record()
    seq = LongSequence{A}(min_size)
    
    @info "Building long read datastore from FASTQ file"
    
    @info "Writing long reads to datastore"
    
    while !eof(rdr)
        read!(rdr, record)
        seq_len = FASTQ.seqlen(record)
        if seq_len < min_size
            discarded = discarded + 1
            continue
        end
        offset = position(ofs)
        resize!(seq, seq_len)
        copyto!(seq, record)
        push!(read_to_file_position, ReadPosSize(offset, seq_len))
        write(ofs, seq.data)
    end
    
    @info "Done writing long read sequences to datastore"
    @info string(discarded, " reads were discarded due to a too short sequence")
    
    fpos = UInt64(position(ofs))
    
    @info "Writing index to datastore"
    
    write_flat_vector(ofs, read_to_file_position)
    
    # Go to the top and dump the number of reads and the position of the index.
    seek(ofs, sizeof(ReadDatastoreMAGIC) + sizeof(Filetype) + sizeof(LongDS_Version) + sizeof(bps))
    write(ofs, fpos)
    close(ofs)
    
    @info string("Built long read datastore with ", length(read_to_file_position), " reads") 
    
    stream = open(outfile * ".loseq", "r+")
    return LongReads{A}(outfile * ".loseq", Symbol(name), Symbol(name), read_to_file_position, stream)
end

function Base.open(::Type{LongReads{A}}, filename::String, name::Union{String,Symbol,Nothing} = nothing) where {A<:DNAAlphabet}
    fd = open(filename, "r")
    
    __validate_datastore_header(fd, LongReads{A})
    
    fpos = read(fd, UInt64)
    default_name = Symbol(readuntil(fd, '\0'))
    seek(fd, fpos)
    read_to_file_position = read_flat_vector(fd, ReadPosSize)
    return LongReads{A}(filename, isnothing(name) ? default_name : Symbol(name), default_name, read_to_file_position, fd)
end

###
### Getting a sequence
###

Base.length(lrds::LongReads) = length(lrds.read_to_file_positions)

Base.summary(io::IO, lrds::LongReads) = print(io, "Long Read Datastore '", lrds.name, "': ", length(lrds), " reads")

function Base.show(io::IO, lrds::LongReads)
    summary(io, lrds)
end

@inline _inbounds_index_of_sequence(lrds::LongReads, idx::Integer) = @inbounds lrds.read_to_file_positions[idx]

@inline function inbounds_load_sequence!(lrds::LongReads{A}, pos::ReadPosSize, seq::LongSequence{A}) where {A<:DNAAlphabet}
    seek(stream(lrds), pos.offset)
    resize!(seq, pos.sequence_size)
    return _load_sequence_data!(lrds, seq)
end

@inline function inbounds_load_sequence!(lrds::LongReads{A}, idx::Integer, seq::LongSequence{A}) where {A<:DNAAlphabet}
    pos_size = _inbounds_index_of_sequence(lrds, idx)
    return inbounds_load_sequence!(lrds, pos_size, seq)
end

@inline function load_sequence!(lrds::LongReads{A}, idx::Integer, seq::LongSequence{A}) where {A<:DNAAlphabet}
    checkbounds(lrds, idx)
    return inbounds_load_sequence!(lrds, idx, seq)
end

@inline function Base.getindex(lrds::LongReads, idx::Integer)
    @boundscheck checkbounds(lrds, idx)
    pos_size = _inbounds_index_of_sequence(lrds, idx)
    seq = eltype(lrds)(pos_size.sequence_size)
    return inbounds_load_sequence!(lrds, pos_size, seq)
end
