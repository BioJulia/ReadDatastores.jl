@testset "Paired read datastores" begin
    function get_fastq_seqs(file)
        seqs = map(open(FASTQ.Reader, file) do rdr
            collect(rdr)
        end) do rec
            FASTQ.sequence(LongDNASeq, rec)
        end
        for seq in seqs
            if length(seq) > 300
                resize!(seq, 300)
            end
        end
        return seqs
    end
    
    
    
    function check_round_trip(R1, R2)
        r1_seqs = get_fastq_seqs(R1)
        r2_seqs = get_fastq_seqs(R2)
        
        fqa = open(FASTQ.Reader, R1)
        fqb = open(FASTQ.Reader, R2)
        
        ds = PairedReads(fqa, fqb,
                         "ecoli-pe.prds", "ecoli-pe",
                         UInt64(250), UInt64(300), UInt64(0), FwRv)
        
        ds_seqs = collect(ds)
        
        return ds_seqs[1:2:end] == r1_seqs && ds_seqs[2:2:end] == r2_seqs
    end
    
    @test check_round_trip("ecoli_tester_R1.fastq", "ecoli_tester_R2.fastq")

    ds = open(PairedReads, "ecoli-pe.prds")
    @test_throws BoundsError ds[100]
    @test_throws BoundsError buffer(ds)[100]
    @test_throws BoundsError load_sequence!(ds, 100, dna"")
    @test collect(ds) == collect(buffer(ds))
    
    
end