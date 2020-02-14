@enum PairedReadOrientation::UInt64 FwRv=1 RvFw=2

struct PairedReads{A<:DNAAlphabet} <: ReadDatastore{LongSequence{A}}
    filename::String         # Filename datastore was opened from.
    name::Symbol             # Name of the datastore. Useful for other applications.
    defaultname::Symbol      # Default name, useful for other applications.
    max_read_len::UInt64     # Maximum size of any read in this datastore.
    chunksize::UInt64        # Number of chunks (UInt64) of sequence data per read.
    fragsize::UInt64         # Fragment size of library.
    readpos_offset::UInt64
    size::UInt64
    orientation::PairedReadOrientation
    stream::IOStream
end

"Get the orientation of the read pairs"
@inline orientation(prds::PairedReads) = prds.orientation
@inline BioSequences.bits_per_symbol(prds::PairedReads{A}) where {A<:DNAAlphabet} = BioSequences.bits_per_symbol(A())

###
### PairedReads Header format on disk.
###

# | Field             | Value  | Type        |
# |:-----------------:|:------:|:-----------:|
# | Magic number      | 0x05D5 | UInt16      | 2
# | Datastore type    | 0x0001 | UInt16      | 2
# | Version number    | 0x0001 | UInt16      | 2
# | Default name      | N/A    | String      | 7
# | Maximum read size | N/A    | UInt64      | 8
# | Chunk size        | N/A    | UInt64      | 8
# | Fragment size     | N/A    | UInt64      | 8
# | Orientation       | N/A    | Orientation | 8
# | Bits Per Nuc      | N/A    | UInt64      | 8
# | Number of Reads   | N/A    | UInt64      | 8

const PairedDS_Version = 0x0002

"""
    PairedReads(rdrx::FASTQ.Reader, rdry::FASTQ.Reader, outfile::String, name::String, minsize::Integer, maxsize::Integer, fragsize::Integer, orientation::PairedReadOrientation)

Construct a Paired Read Datastore from a pair of FASTQ file readers.

Paired-end sequencing reads typically come in the form of two FASTQ files, often
named according to a convention `*_R1.fastq` and `*_R2.fastq`.
One file contains all the "left" sequences of each pair, and the other contains
all the "right" sequences of each pair. The first read pair is made of the first
record in each file.

# Arguments
- `rdrx::FASTQ.Reader`: The reader of the `*_R1.fastq` file.
- `rdxy::FASTQ.Reader`: The reader of the `*_R2.fastq` file.
- `outfile::String`: A prefix for the datastore's filename, the full filename
  will include a ".prseq" extension, which will be added automatically.
- `name::String`: A string denoting a default name for your datastore.
  Naming datastores is useful for downstream applications.
- `minsize::Integer`: A minimum read length (in base pairs). When building
  the datastore, if any pair of reads has one or both reads shorter than this
  cutoff, then the pair will be discarded.
- `maxsize::Integer`: A maximum read length (in base pairs). When building
  the datastore, if any read has a greater length, it will be resized to this
  maximum length and added to the datastore.
- `fragsize::Integer`: The average fragment length of the paired end library
  that was sequenced. This value is entirely optional, but may be important for
  downstream applications.
- `orientation::PairedReadOrientation`: The orientation of the reads. Set it to
  `FwRv` for building a datastore from a sequenced paired end library, and set
  it to `RvFw` if you are building the datastore from reads sequenced from a
  long mate pair library.

# Examples
```jldoctest
julia> using FASTX, ReadDatastores

julia> fwq = open(FASTQ.Reader, "test/ecoli_tester_R1.fastq")
FASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)

julia> rvq = open(FASTQ.Reader, "test/ecoli_tester_R2.fastq")
FASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)

julia> ds = PairedReads(fwq, rvq, "ecoli-test-paired", "my-ecoli-test", 250, 300, 0, FwRv)
[ Info: Building paired read datastore from FASTQ files
[ Info: Writing paired reads to datastore
[ Info: Done writing paired read sequences to datastore
[ Info: 0 read pairs were discarded due to a too short sequence
[ Info: 14 reads were truncated to 300 base pairs
[ Info: Created paired sequence datastore with 10 sequence pairs
Paired Read Datastore 'my-ecoli-test': 20 reads (10 pairs)

```
"""
function PairedReads{A}(rdrx::FASTQ.Reader, rdry::FASTQ.Reader, outfile::String, name::Union{String,Symbol},
                     minsize::Integer, maxsize::Integer, fragsize::Integer, orientation::PairedReadOrientation) where {A<:DNAAlphabet}
    return PairedReads{A}(rdrx, rdry, outfile, name, convert(UInt64, minsize),
                       convert(UInt64, maxsize), convert(UInt64, fragsize), orientation)
end

function PairedReads{A}(rdrx::FASTQ.Reader, rdry::FASTQ.Reader,
                     outfile::String, name::Union{Symbol,String},
                     minsize::UInt64, maxsize::UInt64,
                     fragsize::UInt64, orientation::PairedReadOrientation) where {A<:DNAAlphabet}
    
    # Create and allocate the sequence and record objects.
    lread = FASTQ.Record()
    rread = FASTQ.Record()
    lseq = LongSequence{A}(maxsize)
    rseq = LongSequence{A}(maxsize)
    
    #chunksize::UInt64 = BioSequences.seq_data_len(DNAAlphabet{4}, maxsize)
    chunksize::UInt64 = length(BioSequences.encoded_data(lseq))
    bps = UInt64(BioSequences.bits_per_symbol(A()))
    
    fd = open(outfile * ".prseq", "w")
    
    # Write magic no, datastore type, version number.
    sizepos = write(fd, ReadDatastoreMAGIC, PairedDS, PairedDS_Version) +
    # Write the default name of the datastore.
    writestring(fd, String(name)) +
    # Write the read size, and chunk size.
    write(fd, maxsize, chunksize, fragsize, orientation, bps)
    # Write space for size variable (or number of read pairs).
    readpos = write(fd, UInt64(0)) + sizepos
    
    pairs = discarded = truncated = 0
    
    @info "Building paired read datastore from FASTQ files"
    @info "Writing paired reads to datastore"
    
    while !eof(rdrx) && !eof(rdry)
        # Read in the two records.
        read!(rdrx, lread)
        read!(rdry, rread)
        
        llen = UInt64(FASTQ.seqlen(lread))
        rlen = UInt64(FASTQ.seqlen(rread))
        # If either read is too short, discard them both.
        if llen < minsize || rlen < minsize
            discarded += 1
            continue
        end
        
        if llen > maxsize
            truncated += 1
            ln = maxsize
        else
            ln = llen
        end
        
        if rlen > maxsize
            truncated += 1
            rn = maxsize
        else
            rn = rlen
        end
        
        pairs += 1
        
        # Copy FASTQ records to sequence variables, thus encoding
        # them in 2 bit format.
        copyto!(lseq, 1, lread, 1, ln)
        copyto!(rseq, 1, rread, 1, rn)
        
        # Write the two sizes reads one after the other: 
        # For each sequence, write the size, and then the datachunk.
        
        write(fd, ln)
        write(fd, lseq.data)
        
        write(fd, rn)
        write(fd, rseq.data)
    end
    
    nreads::UInt64 = pairs * 2
    
    seek(fd, sizepos)
    write(fd, nreads)
    
    close(fd)
    
    @info "Done writing paired read sequences to datastore"
    @info string(discarded, " read pairs were discarded due to a too short sequence")
    @info string(truncated, " reads were truncated to ", maxsize, " base pairs")
    @info string("Created paired sequence datastore with ", pairs, " sequence pairs")
    
    stream = open(outfile * ".prseq", "r+")
    return PairedReads{A}(outfile * ".prseq", Symbol(name), Symbol(name), maxsize, chunksize, fragsize,
                       readpos, nreads, orientation, stream)
end

function Base.open(::Type{PairedReads{A}}, filename::String, name::Union{String,Symbol,Nothing} = nothing) where {A<:DNAAlphabet}
    fd = open(filename, "r")
    magic = read(fd, UInt16)
    dstype = reinterpret(Filetype, read(fd, UInt16))
    version = read(fd, UInt16)
    
    @assert magic == ReadDatastoreMAGIC
    @assert dstype == PairedDS
    @assert version == PairedDS_Version
    
    default_name = Symbol(readuntil(fd, '\0'))
    
    max_read_len = read(fd, UInt64)
    chunksize = read(fd, UInt64)
    fragsize = read(fd, UInt64)
    orientation = reinterpret(PairedReadOrientation, read(fd, UInt64))
    bps = read(fd, UInt64)
    @assert bps == BioSequences.bits_per_symbol(A())
    nreads = read(fd, UInt64)
    readpos_offset = position(fd)
    
    return PairedReads{A}(filename, isnothing(name) ? default_name : Symbol(name), default_name,
                       max_read_len, chunksize, fragsize,
                       readpos_offset, nreads, orientation, fd)
end

@inline Base.length(prds::PairedReads) = prds.size

Base.summary(io::IO, prds::PairedReads) = print(io, "Paired Read Datastore '", prds.name, "': ", length(prds), " reads (", div(length(prds), 2), " pairs)")

function Base.show(io::IO, prds::PairedReads)
    summary(io, prds)
end

bytes_per_read(prds::PairedReads) = (prds.chunksize + 1) * sizeof(UInt64)

@inline _inbounds_index_of_sequence(prds::PairedReads, idx::Integer) = prds.readpos_offset + (bytes_per_read(prds) * (idx - 1))

@inline function Base.getindex(prds::PairedReads, idx::Integer)
    @boundscheck checkbounds(prds, idx)
    seq = eltype(prds)(prds.max_read_len)
    return inbounds_load_sequence!(prds, idx, seq)
end
