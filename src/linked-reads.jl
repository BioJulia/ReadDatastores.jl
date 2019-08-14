
@enum LinkedReadsFormat UCDavis=1 Raw=2

const LinkedTag = UInt32

mutable struct LinkedReadData
    seq1::LongSequence{DNAAlphabet{4}}
    seq2::LongSequence{DNAAlphabet{4}}
    tag::LinkedTag
end

Base.isless(a::LinkedReadData, b::LinkedReadData) = a.tab < b.tag
LinkedReadData(len) = LinkedReadData(LongDNASeq(len), LongDNASeq(len), zero(LinkedTag))
    

function build_from_fastq(fwq::FASTQ.Reader, rvq::FASTQ.Reader, outfile::String, name::String, format::LinkedReadsFormat, readsize::Int, chunksize::Int)
    read_tag = Vector{UInt32}()
    
    @info string("Building tag sorted chunks of ", chunksize, " pairs")
    chunk_files = Vector{String}()
    
    
    fwrec = FASTQ.Record()
    rvrec = FASTQ.Record()
    
    chunk_data = [LinkedReadData(readsize) for _ in 1:chunksize]
    #chunk_files = Vector{IO}()
    
    while !eof(fwq) && !eof(rvq)
        # Read in `chunksize` read pairs.
        n_pairs = 0
        while !eof(fwq) && !eof(rvq) && n_pairs < chunksize
            read!(fwq, fwrec)
            read!(rvq, rvrec)
            
        end
        
        # Sort the linked reads by tag
        sort!(chunk_data)
        # Dump the data to a chunk dump
        chunk_fd = open(string("sorted_chunk_", length(chunk_files), ".data"), "w")
        
        
        
    end
    
end