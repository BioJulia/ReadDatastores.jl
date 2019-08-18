abstract type LinkedReadsFormat end

struct UCDavisTenX <: LinkedReadsFormat end

const UCDavis10x = UCDavisTenX()

const LinkedTag = UInt32
mutable struct LinkedReadData
    seq1::LongSequence{DNAAlphabet{4}}
    seq2::LongSequence{DNAAlphabet{4}}
    seqlen1::UInt64
    seqlen2::UInt64
    tag::LinkedTag
end

Base.isless(a::LinkedReadData, b::LinkedReadData) = a.tag < b.tag
LinkedReadData(len) = LinkedReadData(LongDNASeq(len), LongDNASeq(len), zero(UInt64), zero(UInt64), zero(LinkedTag))

const LinkedDS_Version = 0x0001

function _extract_tag_and_sequences!(current_data::LinkedReadData, fwrec::FASTQ.Record, rvrec::FASTQ.Record, readsize::UInt64, ::UCDavisTenX)
    fwid = FASTQ.identifier(fwrec)
    newtag = zero(UInt32)
    @inbounds for i in 1:16
        newtag <<= 2
        letter = fwid[i]
        if letter == 'C'
            newtag = newtag + 1
        elseif letter == 'G'
            newtag = newtag + 2
        elseif letter == 'T'
            newtag = newtag + 3
        elseif letter != 'A'
            newtag = zero(UInt32)
            break
        end
    end
    current_data.tag = newtag
    current_data.seqlen1 = UInt64(min(readsize, FASTQ.seqlen(fwrec)))
    current_data.seqlen2 = UInt64(min(readsize, FASTQ.seqlen(rvrec)))
    copyto!(current_data.seq1, 1, fwrec, 1, current_data.seqlen1)
    copyto!(current_data.seq2, 1, rvrec, 1, current_data.seqlen2)
end

struct LinkedReads <: ReadDatastore{LongSequence{DNAAlphabet{4}}}
    filename::String          # Filename datastore was opened from.
    name::String              # Name of the datastore. Useful for other applications.
    defaultname::String       # Default name, useful for other applications.
    readsize::UInt64          # Maximum size of any read in this datastore.
    chunksize::UInt64
    readpos_offset::UInt64    #  
    read_tags::Vector{UInt32} #
    stream::IO
end

@inline stream(lrds::LinkedReads) = lrds.stream

function LinkedReads(fwq::FASTQ.Reader, rvq::FASTQ.Reader, outfile::String, name::String, format::LinkedReadsFormat, readsize::UInt64, chunksize::Int = 1000000)
    n_pairs = 0
    chunk_files = String[]
    
    @info string("Building tag sorted chunks of ", chunksize, " pairs")
    
    fwrec = FASTQ.Record()
    rvrec = FASTQ.Record()
    chunk_data = [LinkedReadData(readsize) for _ in 1:chunksize]
    datachunksize = length(BioSequences.encoded_data(first(chunk_data).seq1))
    
    while !eof(fwq) && !eof(rvq)
        # Read in `chunksize` read pairs.
        chunkfill = 0
        while !eof(fwq) && !eof(rvq) && chunkfill < chunksize
            read!(fwq, fwrec)
            read!(rvq, rvrec)
            cd_i = chunk_data[chunkfill + 1]
            _extract_tag_and_sequences!(cd_i, fwrec, rvrec, readsize, format)
            if cd_i.tag != zero(UInt32)
                chunkfill = chunkfill + 1
            end
        end
        # Sort the linked reads by tag
        sort!(view(chunk_data, 1:chunkfill))
        # Dump the data to a chunk dump
        chunk_fd = open(string("sorted_chunk_", length(chunk_files), ".data"), "w")
        @info string("Dumping ", chunkfill, " tag-sorted read pairs to chunk ", length(chunk_files))
        for j in 1:chunkfill
            cd_j = chunk_data[j]
            write(chunk_fd, cd_j.tag)
            write(chunk_fd, cd_j.seqlen1)
            write(chunk_fd, BioSequences.encoded_data(cd_j.seq1))
            write(chunk_fd, cd_j.seqlen2)
            write(chunk_fd, BioSequences.encoded_data(cd_j.seq2))
        end
        close(chunk_fd)
        push!(chunk_files, string("sorted_chunk_", length(chunk_files), ".data"))
        @info "Dumped"
        n_pairs = n_pairs + chunkfill
        @info string("Processed ", n_pairs, " read pairs so far") 
    end
    close(fwq)
    close(rvq)
    @info "Finished building tag sorted chunks"
    
    @info "Performing merge from disk"
    read_tag = zeros(UInt32, n_pairs)
    @info string("Leaving space for ", n_pairs, " read_tag entries")
    
    output = open(outfile, "w")
    # Write magic no, datastore type, version number.
    write(output, ReadDatastoreMAGIC, LinkedDS, LinkedDS_Version)
    # Write the default name of the datastore.
    writestring(output, name)
    # Write the read size, and chunk size.
    write(output, readsize, datachunksize)
    
    read_tag_offset = position(output)
    write_flat_vector(output, read_tag)
    
    # Multiway merge of chunks...
    chunk_fds = map(x -> open(x, "r"), chunk_files)
    openfiles = length(chunk_fds)
    next_tags = Vector{UInt32}(undef, openfiles)
    for i in eachindex(chunk_fds)
        next_tags[i] = read(chunk_fds[i], LinkedTag)
    end
    
    #filebuffer = Vector{UInt8}(undef, sizeof(LinkedTag) + 2(sizeof(UInt64) + datachunksize))
    #filebuffer = Vector{UInt8}(undef, 2(sizeof(UInt64) + datachunksize))
    filebuffer = Vector{UInt8}(undef, 16(datachunksize + 1))
    empty!(read_tag)
    while openfiles > 0
        mintag = minimum(next_tags)
        # Copy from each file until then next tag is > mintag
        for i in eachindex(chunk_fds)
            fd = chunk_fds[i]
            while !eof(fd) && next_tags[i] == mintag
                # Read from chunk file, and write to final file.
                readbytes!(fd, filebuffer)
                write(output, filebuffer)
                push!(read_tag, mintag)
                if eof(fd)
                    openfiles = openfiles - 1
                    next_tags[i] = typemax(LinkedTag)
                    @info string("Chunk ", i, " is finished")
                else
                    next_tags[i] = read(fd, LinkedTag)
                end
            end
        end
    end
    
    @info "Finished merge from disk"
    seek(output, read_tag_offset)
    @info string("Writing down ", n_pairs, " read tag entries")
    write_flat_vector(output, read_tag)
    readspos = position(output)
    close(output)
    for fd in chunk_fds
        close(fd)
    end
    for fname in chunk_files
        rm(fname)
    end
    @info string("Created linked sequence datastore with ", n_pairs, " sequence pairs")
    return LinkedReads(outfile, name, name, readsize, datachunksize, readspos, read_tag, open(outfile, "r"))
end

function Base.open(::Type{LinkedReads}, filename::String)
    fd = open(filename, "r")
    magic = read(fd, UInt16)
    dstype = reinterpret(Filetype, read(fd, UInt16))
    version = read(fd, UInt16)
    
    @assert magic == ReadDatastoreMAGIC
    @assert dstype == LinkedDS
    @assert version == LinkedDS_Version
    
    default_name = readuntil(fd, '\0')
    readsize = read(fd, UInt64)
    chunksize = read(fd, UInt64)
    
    read_tags = read_flat_vector(fd, UInt32)
    
    return LinkedReads(filename, default_name, default_name, readsize, chunksize, position(fd), read_tags, fd)
end

bytes_per_read(lrds::LinkedReads) = (lrds.chunksize + 1) * sizeof(UInt64)
@inline _inbounds_index_of_sequence(lrds::LinkedReads, idx::Integer) = lrds.readpos_offset + (bytes_per_read(lrds) * (idx - 1))

@inline Base.length(lrds::LinkedReads) = length(lrds.read_tags) * 2

@inline function Base.getindex(lrds::LinkedReads, idx::Integer)
    @boundscheck checkbounds(lrds, idx)
    seq = eltype(lrds)(lrds.readsize)
    return inbounds_load_sequence!(lrds, idx, seq)
end