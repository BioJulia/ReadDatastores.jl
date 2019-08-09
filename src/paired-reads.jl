@enum PairedReadOrientation::UInt64 FwRv=1 RvFw=2

struct PairedReadDatastore
    filename::String         # Filename datastore was opened from.
    name::String             # Name of the datastore. Useful for other applications.
    defaultname::String      # Default name, useful for other applications.
    readsize::UInt64         # Maximum size of any read in this datastore.
    chunksize::UInt64        # Number of chunks of sequence data per read.
    fragsize::UInt64         # Fragment size of library.
    readpos_offset::UInt64
    size::UInt64
    orientation::PairedReadOrientation
    stream::IOStream
end

orientation(prds::PairedReadDatastore) = prds.orientation
maxseqlen(prds::PairedReadDatastore) = prds.readsize
name(prds::PairedReadDatastore) = prds.name

const PRDS = PairedReadDatastore

###
### PairedReadDatastore Header
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
# | Number of Reads   | N/A    | UInt64      | 8

const PairedDS_Version = 0x0001

function PairedReadDatastore(rdrx::FASTQ.Reader, rdry::FASTQ.Reader,
                            outfile::String, name::String,
                            minsize::UInt64, maxsize::UInt64,
                            fragsize::UInt64, orientation::PairedReadOrientation)
    
    # Create and allocate the sequence and record objects.
    lread = FASTQ.Record()
    rread = FASTQ.Record()
    lseq = LongSequence{DNAAlphabet{4}}(maxsize)
    rseq = LongSequence{DNAAlphabet{4}}(maxsize)
    
    #chunksize::UInt64 = BioSequences.seq_data_len(DNAAlphabet{4}, maxsize)
    chunksize::UInt64 = length(BioSequences.encoded_data(lseq))
    
    fd = open(outfile, "w")
    
    # Write magic no, datastore type, version number.
    sizepos = write(fd, SeqDataStoreMAGIC, PairedDS, PairedDS_Version) +
    # Write the default name of the datastore.
    writestring(fd, name) +
    # Write the read size, and chunk size.
    write(fd, maxsize, chunksize, fragsize, orientation)
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
    
    stream = open(outfile, "r+")
    return PairedReadDatastore(outfile, name, name,
                              maxsize, chunksize, fragsize,
                              readpos, nreads, orientation, stream)
end

function Base.open(::Type{PairedReadDatastore}, filename::String)
    fd = open(filename, "r")
    magic = read(fd, UInt16)
    dstype = reinterpret(Filetype, read(fd, UInt16))
    version = read(fd, UInt16)
    
    @assert magic == SeqDataStoreMAGIC
    @assert dstype == PairedDS
    @assert version == PairedDS_Version
    
    default_name = readuntil(fd, '\0')
    
    readsize = read(fd, UInt64)
    chunksize = read(fd, UInt64)
    fragsize = read(fd, UInt64)
    orientation = reinterpret(PairedReadOrientation, read(fd, UInt64))
    nreads = read(fd, UInt64)
    readpos_offset = position(fd)
    
    @info magic
    @info dstype
    @info version
    @info readsize
    @info chunksize
    @info fragsize
    @info orientation
    @info nreads
    @info readpos_offset
    
    return PairedReadDatastore(filename, default_name, default_name,
                              readsize, chunksize, fragsize,
                              readpos_offset, nreads, orientation, fd)
end

@inline Base.length(prds::PairedReadDatastore) = prds.size

Base.summary(io::IO, prds::PairedReadDatastore) = print(io, "Paired Sequence Datastore '", prds.name, "': ", length(prds), " reads")

function Base.show(io::IO, prds::PairedReadDatastore)
    summary(io, prds)
end

bytes_per_read(prds::PRDS) = (prds.chunksize + 1) * sizeof(UInt64)
@inline unsafe_read_offset_in_file(prds::PRDS, idx::Integer) = prds.readpos_offset + (bytes_per_read(prds) * (idx - 1))

@inline function unsafe_load_read!(prds::PRDS, idx::Integer, seq::LongSequence{DNAAlphabet{4}})
    seek(prds.stream, unsafe_read_offset_in_file(prds, idx))
    seqlen = read(prds.stream, UInt64)
    resize!(seq, seqlen)
    unsafe_read(prds.stream, pointer(seq.data), length(seq.data) * sizeof(UInt64))
    return seq
end

@inline function load_read!(prds::PRDS, idx::Integer, seq::LongSequence{DNAAlphabet{4}})
    checkbounds(prds, idx)
    return unsafe_load_read!(prds, idx, seq)
end

@inline function Base.checkbounds(prds::PRDS, i::Integer)
    if 1 ≤ i ≤ lastindex(prds)
        return true
    end
    throw(BoundsError(prds, i))
end

firstindex(prds::PRDS) = 1
lastindex(prds::PRDS) = length(prds)
eachindex(prds::PRDS) = Base.OneTo(lastindex(prds))

@inline function Base.getindex(prds::PRDS, idx::Integer)
    @boundscheck checkbounds(prds, idx)
    seq = LongDNASeq(prds.readsize)
    return unsafe_load_read!(prds, idx, seq)
end

Base.IteratorSize(prds::PRDS) = Base.HasLength()
Base.IteratorEltype(prds::PRDS) = Base.HasEltype()
Base.eltype(prds::PRDS) = LongSequence{DNAAlphabet{4}}

@inline function Base.iterate(prds::PRDS, state = 1)
    @inbounds if firstindex(prds) ≤ state ≤ lastindex(prds)
        return prds[state], state + 1
    else
        return nothing
    end
end