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
        ds = LongReadDatastore(fq, "human-nanopore.lrds", "human-nanopore", UInt64(0))
        ds2 = open(LongReadDatastore, "human-nanopore.lrds")
        ds_seqs = collect(ds)
        ds2_seqs = collect(ds2)
        return ds_seqs == seqs == ds2_seqs
    end
    
    @testset "Human oxford nanopore 2D consensus reads tester" begin
        @test check_round_trip("human_nanopore_tester_2D.fastq")
    
        ds = open(LongReadDatastore, "human-nanopore.lrds")
        ds2 = open(LongReadDatastore, "human-nanopore.lrds")
    
        @test ReadDatastores.index(ds) == ReadDatastores.index(ds2)
        @test ReadDatastores.index(ds)[1] == ReadDatastores.index(ds2)[1]
        @test firstindex(ds) == firstindex(ds2) == 1
        @test lastindex(ds) == lastindex(ds2) == 10
        
        @test Base.IteratorSize(ds) == Base.IteratorSize(ds2) == Base.HasLength()
        @test Base.IteratorEltype(ds) == Base.IteratorEltype(ds2) == Base.HasEltype()
        @test Base.eltype(ds) == Base.eltype(ds2) == LongSequence{DNAAlphabet{4}}

        
    end
end