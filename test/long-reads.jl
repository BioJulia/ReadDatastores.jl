@testset "Long read datastores" begin
    function get_fastq_seqs(file)
        seqs = map(open(FASTQ.Reader, file) do rdr
            collect(rdr)
        end) do rec
            FASTQ.sequence(LongDNASeq, rec)
        end
        return seqs
    end
    
    function check_round_trip(FQ)
        seqs = get_fastq_seqs(FQ)
        fq = open(FASTQ.Reader, FQ)
        ds = LongReadDatastore(fq, "human-nanopore.prds", "human-nanopore", UInt64(0))
        ds_seqs = collect(ds)
        return ds_seqs == seqs
    end
    
    @test check_round_trip("human_nanopore_tester_2D.fastq")
end