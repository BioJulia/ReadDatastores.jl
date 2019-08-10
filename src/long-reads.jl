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
# | Index position in file        | N/A    | UInt64      | 8
# | Default name of the datastore | N/A    | String      | N

function LongReadDatastore(rdr::FASTQ.Reader, outfile::String, name::String, min_size::UInt64)
    discarded = 0
    
    read_to_file_position = Vector{ReadPosSize}()
    ofs = open(outfile, "w")
    
    write(ofs, SeqDataStoreMAGIC, LongDS, LongDS_Version, zero(UInt64))
    
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
    end
    
    @info "Done writing paired read sequences to datastore"
    @info string(discarded, " reads were discarded due to a too short sequence")
    
    fpos = UInt64(position(ofs))
    
    @info "Writing index to datastore"
    
    write_flat_vector(ofs, read_to_file_position)
    
    # Go to the top and dump the number of reads and the position of the index.
    seek(ofs, sizeof(SeqDataStoreMAGIC) + sizeof(Filetype) + sizeof(LongDS_Version))
    write(ofs, fpos)
    close(ofs)
    
    @info string("Built long read datastore with ", length(read_to_file_position), " reads") 
    
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
    
    fpos = read(fd, UInt64)
    
    default_name = readuntil(fd, '\0')
    
    seek(fd, fpos)
    
    read_to_file_position = read_flat_vector(fd, ReadPosSize)
    
    return LongReadDatastore(filename, default_name, default_name, read_to_file_position, fd)
end

###
### Getting a sequence
###

Base.length(lrds::LRDS) = length(lrds.read_to_file_positions)

firstindex(lrds::LRDS) = 1
lastindex(lrds::LRDS) = length(lrds)
eachindex(lrds::LRDS) = Base.OneTo(lastindex(lrds))

@inline function Base.checkbounds(lrds::LRDS, i::Integer)
    if firstindex(lrds) ≤ i ≤ lastindex(lrds)
        return true
    end
    throw(BoundsError(lrds, i))
end

@inbounds inbounds_position_and_size(lrds::LRDS, idx::Integer) = @inbounds lrds.read_to_file_positions[idx]

@inbounds function position_and_size(lrds::LRDS, idx::Integer)
    checkbounds(lrds, idx)
    return inbounds_position_and_size(lrds, idx)
end

@inline function unsafe_load_read!(lrds::LRDS, pos_size::ReadPosSize, seq::LongSequence{DNAAlphabet{4}})
    seek(lrds.stream, pos_size.offset)
    resize!(seq, pos_size.sequence_size)
    unsafe_read(lrds.stream, pointer(seq.data), length(seq.data) * sizeof(UInt64))
    return seq
end

@inline function inbounds_load_read!(lrds::LRDS, idx::Integer, seq::LongSequence{DNAAlphabet{4}})
    pos_size = inbounds_position_and_size(lrds, idx)
    return unsafe_load_read!(lrds, pos_size, seq)
end

@inline function load_read!(lrds::LRDS, idx::Integer, seq::LongSequence{DNAAlphabet{4}})
    checkbounds(lrds, idx)
    return inbounds_load_read!(lrds, idx, seq)
end

@inline function Base.getindex(lrds::LRDS, idx::Integer)
    @boundscheck checkbounds(lrds, idx)
    pos_size = inbounds_position_and_size(lrds, idx)
    seq = LongDNASeq(pos_size.sequence_size)
    return unsafe_load_read!(lrds, pos_size, seq)
end

Base.IteratorSize(lrds::LRDS) = Base.HasLength()
Base.IteratorEltype(lrds::LRDS) = Base.HasEltype()
Base.eltype(lrds::LRDS) = LongSequence{DNAAlphabet{4}}

@inline function Base.iterate(lrds::LRDS, state = 1)
    @inbounds if firstindex(lrds) ≤ state ≤ lastindex(lrds)
        return lrds[state], state + 1
    else
        return nothing
    end
end