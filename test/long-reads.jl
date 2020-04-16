@testset "Long read datastores" begin
    function get_fastq_seqs(file)
        seqs = map(open(FASTQ.Reader, file) do rdr
            collect(rdr)
        end) do rec
            FASTQ.sequence(LongDNASeq, rec)
        end
        return seqs
    end
    
    function check_round_trip(::Type{A}, FQ) where {A<:DNAAlphabet}
        seqs = get_fastq_seqs(FQ)
        fq = open(FASTQ.Reader, FQ)
        ds = LongReads{A}(fq, "human-nanopore", "human-nanopore", 0)
        ds2 = open(LongReads{A}, "human-nanopore.loseq")
        ds_seqs = collect(ds)
        ds2_seqs = collect(ds2)
        return ds_seqs == seqs == ds2_seqs
    end
    
    function check_show(ds, msg)
        buf = IOBuffer()
        show(buf, ds)
        return String(take!(buf)) == msg
    end
    
    function run_tests(::Type{A}) where {A<:DNAAlphabet}
        @test check_round_trip(A, "human_nanopore_tester_2D.fastq")
        ds = open(LongReads{A}, "human-nanopore.loseq")
        ds2 = open(LongReads{A}, "human-nanopore.loseq")
        @test ReadDatastores.deduce_datastore_type("human-nanopore.loseq") == LongReads{A}
        @test ReadDatastores.index(ds) == ReadDatastores.index(ds2)
        @test ReadDatastores.index(ds)[1] == ReadDatastores.index(ds2)[1]
        @test firstindex(ds) == firstindex(ds2) == 1
        @test lastindex(ds) == lastindex(ds2) == 10
        
        @test check_show(ds, "Long Read Datastore 'human-nanopore': 10 reads")
        
        @test Base.IteratorSize(ds) == Base.IteratorSize(ds2) == Base.HasLength()
        @test Base.IteratorEltype(ds) == Base.IteratorEltype(ds2) == Base.HasEltype()
        @test Base.eltype(ds) == Base.eltype(ds2) == LongSequence{A}
        
        @test collect(ds) == collect(buffer(ds)) == collect(buffer(ds, 1))
        @test ds[5] == buffer(ds)[5] == load_sequence!(ds, 5, LongSequence{A}()) == load_sequence!(buffer(ds), 5, LongSequence{A}())
        
        bds = buffer(ds)
        @test eltype(bds) == eltype(ds)
        @test firstindex(bds) == firstindex(ds)
        @test eachindex(bds) == eachindex(ds)
        @test Base.IteratorSize(bds) == Base.IteratorSize(ds)
        @test Base.IteratorEltype(bds) == Base.IteratorEltype(ds)
    end
    
    @testset "Human oxford nanopore 2D consensus reads tester" begin
        run_tests(DNAAlphabet{4})
        run_tests(DNAAlphabet{2})
    end
end