@testset "Linked read datastores" begin
    
    function ucdavis_tag(rec::FASTQ.Record)
        id = FASTQ.identifier(rec)
        tag = zero(UInt32)
        @inbounds for i in 1:16
            tag <<= 2
            letter = id[i]
            if letter == 'C'
                tag = tag + 1
            elseif letter == 'G'
                tag = tag + 2
            elseif letter == 'T'
                tag = tag + 3
            elseif letter != 'A'
                tag = zero(UInt32)
                break
            end
        end
        return tag
    end
    
    function get_fastq_seqs(r1, r2, maxlen)
        r1s = open(FASTQ.Reader, r1) do rdr
            collect(rdr)
        end
        r2s = open(FASTQ.Reader, r2) do rdr
            collect(rdr)
        end
        tags = map(r1s) do rec
            ucdavis_tag(rec)
        end
        keep = tags .!= zero(UInt32)
        r1s = FASTQ.sequence.(r1s[keep])
        r2s = FASTQ.sequence.(r2s[keep])
        tags = tags[keep]
        for seq in r1s
            if length(seq) > maxlen
                resize!(seq, maxlen)
            end
        end
        for seq in r1s
            if length(seq) > maxlen
                resize!(seq, maxlen)
            end
        end
        order = sortperm(tags)
        return r1s[order], r2s[order]
    end
    
    
    
    function check_round_trip(R1, R2, maxlen, chunksize = 1000000)
        r1_seqs, r2_seqs = get_fastq_seqs(R1, R2, maxlen)
        
        fqa = open(FASTQ.Reader, R1)
        fqb = open(FASTQ.Reader, R2)
        
        ds = LinkedReads(fqa, fqb, "10xtest", Symbol("ucdavis-test"), UCDavis10x, maxlen, chunksize)
        ds2 = open(LinkedReads, "10xtest.lrseq")
        
        ds_seqs = collect(ds)
        ds2_seqs = collect(ds2)
        
        return ds_seqs[1:2:end] == ds2_seqs[1:2:end] == r1_seqs && ds_seqs[2:2:end] == ds2_seqs[2:2:end] == r2_seqs
    end
    
    function check_show(ds, msg)
        buf = IOBuffer()
        show(buf, ds)
        return String(take!(buf)) == msg
    end
    
    @test check_round_trip("10x_tester_R1.fastq", "10x_tester_R2.fastq", 250, 10)
    
    ds = open(LinkedReads, "10xtest.lrseq")
    @test ReadDatastores.name(ds) == Symbol("ucdavis-test")
    @test ReadDatastores.maxseqlen(ds) == 250
    @test check_show(ds, "Linked Read Datastore 'ucdavis-test': 166 reads (83 pairs)")
    @test firstindex(ds) == 1
    @test lastindex(ds) == 166
    @test Base.IteratorSize(ds) == Base.HasLength()
    @test Base.IteratorEltype(ds) == Base.HasEltype()
    @test Base.eltype(ds) == LongSequence{DNAAlphabet{4}}
    @test_throws BoundsError ds[200]
    @test_throws BoundsError buffer(ds)[200]
    @test_throws BoundsError load_sequence!(ds, 200, dna"")
    @test collect(ds) == collect(buffer(ds)) == open(LinkedReads, "10xtest.lrseq") do ds
        collect(ds)
    end
end