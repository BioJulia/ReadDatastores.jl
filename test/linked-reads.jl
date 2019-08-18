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
    
    
    
    function check_round_trip(R1, R2, maxlen)
        r1_seqs, r2_seqs = get_fastq_seqs(R1, R2, maxlen)
        
        fqa = open(FASTQ.Reader, R1)
        fqb = open(FASTQ.Reader, R2)
        
        ds = LinkedReads(fqa, fqb, "10xtest.lrds", "ucavis-test", UCDavis10x, UInt64(maxlen))
        ds2 = open(LinkedReads, "10xtest.lrds")
        
        ds_seqs = collect(ds)
        ds2_seqs = collect(ds2)
        
        return ds_seqs[1:2:end] == ds2_seqs[1:2:end] == r1_seqs && ds_seqs[2:2:end] == ds2_seqs[2:2:end] == r2_seqs
    end
    
    #function check_show(ds, msg)
    #    buf = IOBuffer()
    #    show(buf, ds)
    #    return String(take!(buf)) == msg
    #end
    
    @test check_round_trip("10x_tester_R1.fastq", "10x_tester_R2.fastq", 250)
    
    #=
    ds = open(PairedReads, "ecoli-pe.prds")
    @test ReadDatastores.name(ds) == "ecoli-pe"
    @test ReadDatastores.maxseqlen(ds) == 300
    @test ReadDatastores.orientation(ds) == FwRv
    @test check_show(ds, "Paired Read Datastore 'ecoli-pe': 20 reads")
    @test firstindex(ds) == 1
    @test lastindex(ds) == 20
    @test Base.IteratorSize(ds) == Base.HasLength()
    @test Base.IteratorEltype(ds) == Base.HasEltype()
    @test Base.eltype(ds) == LongSequence{DNAAlphabet{4}}
    @test_throws BoundsError ds[100]
    @test_throws BoundsError buffer(ds)[100]
    @test_throws BoundsError load_sequence!(ds, 100, dna"")
    @test collect(ds) == collect(buffer(ds)) == open(PairedReads, "ecoli-pe.prds") do ds
        collect(ds)
    end
    =#
    
end