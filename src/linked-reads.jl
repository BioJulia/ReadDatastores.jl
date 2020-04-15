abstract type LinkedReadsFormat end

struct UCDavisTenX <: LinkedReadsFormat end
struct RawTenX <: LinkedReadsFormat end

const UCDavis10x = UCDavisTenX()
const Raw10x = RawTenX()

const LinkedTag = UInt32
mutable struct LinkedReadData{A<:DNAAlphabet}
    seq1::LongSequence{A}
    seq2::LongSequence{A}
    seqlen1::UInt64
    seqlen2::UInt64
    tag::LinkedTag
end

Base.isless(a::LinkedReadData, b::LinkedReadData) = a.tag < b.tag
LinkedReadData{A}(len) where {A<:DNAAlphabet} = LinkedReadData{A}(LongSequence{A}(len), LongSequence{A}(len), zero(UInt64), zero(UInt64), zero(LinkedTag))

const LinkedDS_Version = 0x0003

function _extract_tag_and_sequences!(current_data::LinkedReadData, fwrec::FASTQ.Record, rvrec::FASTQ.Record, max_read_len::UInt64, ::UCDavisTenX)
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
    current_data.seqlen1 = UInt64(min(max_read_len, FASTQ.seqlen(fwrec)))
    current_data.seqlen2 = UInt64(min(max_read_len, FASTQ.seqlen(rvrec)))
    copyto!(current_data.seq1, 1, fwrec, 1, current_data.seqlen1)
    copyto!(current_data.seq2, 1, rvrec, 1, current_data.seqlen2)
end

struct LinkedReads{A<:DNAAlphabet} <: ShortReads{A}
    filename::String          # Filename datastore was opened from.
    name::Symbol              # Name of the datastore. Useful for other applications.
    defaultname::Symbol       # Default name, useful for other applications.
    max_read_len::UInt64      # Maximum size of any read in this datastore.
    chunksize::UInt64
    readpos_offset::UInt64    #  
    read_tags::Vector{UInt32} #
    stream::IO
end

"""
    LinkedReads{A}(fwq::FASTQ.Reader, rvq::FASTQ.Reader, outfile::String, name::String, format::LinkedReadsFormat, max_read_len::Integer, chunksize::Int = 1000000) where {A<:DNAAlphabet}

Construct a Linked Read Datastore from a pair of FASTQ file readers.

Paired sequencing reads typically come in the form of two FASTQ files, often
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
- `format::LinkedReadsFormat`: Specify the linked reads format of your fastq files.
  If you have plain 10x files, set format to `Raw10x`. If you have 10x reads
  output in the UCDavis format, set the format to `UCDavis10x`.
- `chunksize::Int = 1000000`: How many read pairs to process per disk batch
  during the tag sorting step of construction.

# Examples
```jldoctest
julia> using FASTX, ReadDatastores

julia> fqa = open(FASTQ.Reader, "test/10x_tester_R1.fastq")
FASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)

julia> fqb = open(FASTQ.Reader, "test/10x_tester_R2.fastq")
FASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)

julia> ds = LinkedReads{DNAAlphabet{2s}}(fqa, fqb, "10xtest", "ucdavis-test", UCDavis10x, 250)
[ Info: Building tag sorted chunks of 1000000 pairs
[ Info: Dumping 83 tag-sorted read pairs to chunk 0
[ Info: Dumped
[ Info: Processed 83 read pairs so far
[ Info: Finished building tag sorted chunks
[ Info: Performing merge from disk
[ Info: Leaving space for 83 read_tag entries
[ Info: Chunk 0 is finished
[ Info: Finished merge from disk
[ Info: Writing 83 read tag entries
[ Info: Created linked sequence datastore with 83 sequence pairs
Linked Read Datastore 'ucdavis-test': 166 reads (83 pairs)

```
"""
function LinkedReads{A}(fwq::FASTQ.Reader, rvq::FASTQ.Reader, outfile::String, name::Union{String,Symbol}, format::LinkedReadsFormat, max_read_len::Integer, chunksize::Int = 1000000) where {A<:DNAAlphabet}
    return LinkedReads{A}(fwq, rvq, outfile, name, format, convert(UInt64, max_read_len), chunksize)
end

function LinkedReads{A}(fwq::FASTQ.Reader, rvq::FASTQ.Reader, outfile::String, name::Union{String,Symbol}, format::LinkedReadsFormat, max_read_len::UInt64, chunksize::Int = 1000000) where {A<:DNAAlphabet}
    n_pairs = 0
    chunk_files = String[]
    
    @info string("Building tag sorted chunks of ", chunksize, " pairs")
    
    fwrec = FASTQ.Record()
    rvrec = FASTQ.Record()
    chunk_data = [LinkedReadData{A}(max_read_len) for _ in 1:chunksize]
    datachunksize = length(BioSequences.encoded_data(first(chunk_data).seq1))
    
    while !eof(fwq) && !eof(rvq)
        # Read in `chunksize` read pairs.
        chunkfill = 0
        while !eof(fwq) && !eof(rvq) && chunkfill < chunksize
            read!(fwq, fwrec)
            read!(rvq, rvrec)
            cd_i = chunk_data[chunkfill + 1]
            _extract_tag_and_sequences!(cd_i, fwrec, rvrec, max_read_len, format)
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
    
    bps = UInt64(BioSequences.bits_per_symbol(A()))
    output = open(outfile * ".lrseq", "w")
    # Write magic no, datastore type, version number, and bits per symbol.
    write(output, ReadDatastoreMAGIC, LinkedDS, LinkedDS_Version, bps)
    # Write the default name of the datastore.
    writestring(output, String(name))
    # Write the read size, chunk size.
    write(output, max_read_len, datachunksize)
    
    read_tag_offset = position(output)
    write_flat_vector(output, read_tag)
    
    # Multiway merge of chunks...
    chunk_fds = map(x -> open(x, "r"), chunk_files)
    openfiles = length(chunk_fds)
    next_tags = Vector{UInt32}(undef, openfiles)
    for i in eachindex(chunk_fds)
        next_tags[i] = read(chunk_fds[i], LinkedTag)
    end
    
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
                    @info string("Chunk ", i - 1, " is finished")
                else
                    next_tags[i] = read(fd, LinkedTag)
                end
            end
        end
    end
    
    @info "Finished merge from disk"
    seek(output, read_tag_offset)
    @info string("Writing ", n_pairs, " read tag entries")
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
    return LinkedReads{A}(outfile * ".lrseq", Symbol(name), Symbol(name), max_read_len, datachunksize, readspos, read_tag, open(outfile * ".lrseq", "r"))
end

function Base.open(::Type{LinkedReads{A}}, filename::String, name::Union{Symbol,String,Nothing} = nothing) where {A<:DNAAlphabet}
    fd = open(filename, "r")
    magic = read(fd, UInt16)
    dstype = reinterpret(Filetype, read(fd, UInt16))
    version = read(fd, UInt16)
    
    @assert magic == ReadDatastoreMAGIC
    @assert dstype == LinkedDS
    @assert version == LinkedDS_Version
    
    bps = read(fd, UInt64)
    @assert bps == UInt64(BioSequences.bits_per_symbol(A()))
    
    default_name = Symbol(readuntil(fd, '\0'))
    max_read_len = read(fd, UInt64)
    chunksize = read(fd, UInt64)
    
    read_tags = read_flat_vector(fd, UInt32)
    
    return LinkedReads{A}(filename, isnothing(name) ? default_name : Symbol(name), default_name, max_read_len, chunksize, position(fd), read_tags, fd)
end

@inline _read_data_begin(prds::LinkedReads) = prds.readpos_offset
@inline _bytes_per_read(prds::LinkedReads) = (prds.chunksize + 1) * sizeof(UInt64)
@inline max_read_length(prds::LinkedReads) = prds.max_read_len

@inline Base.length(lrds::LinkedReads) = length(lrds.read_tags) * 2

Base.summary(io::IO, lrds::LinkedReads) = print(io, "Linked Read Datastore '", lrds.name, "': ", length(lrds), " reads (", div(length(lrds), 2), " pairs)")

@inline inbounds_read_tag(lrds::LinkedReads, idx::Integer) = @inbounds lrds.read_tags[div(idx + 1, 2)]

"""
    read_tag(lr::LinkedReads, i::Integer)

Get the tag for the i'th read
"""
function read_tag(lr::LinkedReads, i::Integer)
    checkbounds(lr, i)
    inbounds_readtag(lr, i)
end

