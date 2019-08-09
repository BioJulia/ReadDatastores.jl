struct ReadPosSize
    offset::UInt64
    sequence_size::UInt64
end

Base.:(==)(x::ReadPosSize, y::ReadPosSize) = x.offset == y.offset && x.sequence_size == y.sequence_size

struct LongReadDatastore
    filename::String
    name::String
    default_name::String
    read_to_file_positions::Vector{ReadPosSize}
    stream::IO
end

const LRDS = LongReadDatastore

index(lrds::LRDS) = lrds.read_to_file_positions

const LongDS_Version = 0x0001

###
### LongReadDatastore Header
###

# | Field                         | Value  | Type        |
# |:-----------------------------:|:------:|:-----------:|
# | Magic number                  | 0x05D5 | UInt16      | 2
# | Datastore type                | 0x0002 | UInt16      | 2
# | Version number                | 0x0001 | UInt16      | 2
# | Number of reads               | N/A    | UInt64      | 8
# | Index position in file        | N/A    | UInt64      | 8
# | Default name of the datastore | N/A    | String      | N

function LongReadDatastore(rdr::FASTQ.Reader, outfile::String, name::String, min_size::UInt64)
    n_reads = one(UInt64)
    discarded = 0
    
    read_to_file_position = ReadPosSize[ReadPosSize(0, 0)]
    ofs = open(outfile, "w")
    
    write(ofs, SeqDataStoreMAGIC, LongDS, LongDS_Version, n_reads, zero(UInt64))
    
    writestring(ofs, name)
    
    record = FASTQ.Record()
    seq = LongSequence{DNAAlphabet{4}}(min_size)
    
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
        n_reads = n_reads + 1
    end
    
    @info "Done writing paired read sequences to datastore"
    @info string(discarded, " reads were discarded due to a too short sequence")
    
    fpos = UInt64(position(ofs))
    
    @info "Writing index to datastore"
    
    write_flat_vector(ofs, read_to_file_position)
    
    # Go to the top and dump the number of reads and the position of the index.
    seek(ofs, sizeof(SeqDataStoreMAGIC) + sizeof(Filetype) + sizeof(LongDS_Version))
    write(ofs, n_reads)
    write(ofs, fpos)
    close(ofs)
    
    @info string("Built long read datastore with ", n_reads - 1, " reads") 
    
    stream = open(outfile, "r+")
    return LongReadDatastore(outfile, name, name, read_to_file_position, stream)
end

function Base.open(::Type{LongReadDatastore}, filename::String)
    fd = open(filename, "r")
    magic = read(fd, UInt16)
    dstype = reinterpret(Filetype, read(fd, UInt16))
    version = read(fd, UInt16)
    
    @assert magic == SeqDataStoreMAGIC
    @assert dstype == LongDS
    @assert version == LongDS_Version
    
    n_reads = read(fd, UInt64)
    
    fpos = read(fd, UInt64)
    
    default_name = readuntil(fd, '\0')
    
    seek(fd, fpos)
    
    read_to_file_position = read_flat_vector(fd, ReadPosSize)
    
    return LongReadDatastore(filename, default_name, default_name, read_to_file_position, fd)
end
    