@testset "Paired read datastores" begin
    function get_fastq_seqs(::Type{A}, file) where {A<:DNAAlphabet}
        seqs = map(open(FASTQ.Reader, file) do rdr
            collect(rdr)
        end) do rec
            FASTQ.sequence(LongSequence{A}, rec)
        end
        for seq in seqs
            if length(seq) > 300
                resize!(seq, 300)
            end
        end
        return seqs
    end
    
    function check_round_trip(::Type{A}, R1, R2) where {A<:DNAAlphabet}
        r1_seqs = get_fastq_seqs(A, R1)
        r2_seqs = get_fastq_seqs(A, R2)
        
        fqa = open(FASTQ.Reader, R1)
        fqb = open(FASTQ.Reader, R2)
        
        ds = PairedReads{A}(fqa, fqb, "ecoli-pe", "ecoli-pe", 250, 300, 0, FwRv)
        
        ds_seqs = collect(ds)
        
        return ds_seqs[1:2:end] == r1_seqs && ds_seqs[2:2:end] == r2_seqs
    end
    
    function check_show(ds, msg)
        buf = IOBuffer()
        show(buf, ds)
        return String(take!(buf)) == msg
    end
    
    function run_checks(::Type{A}) where {A<:DNAAlphabet}
        @test check_round_trip(A, "ecoli_tester_R1.fastq", "ecoli_tester_R2.fastq")
        ds = open(PairedReads{A}, "ecoli-pe.prseq")
        @test ReadDatastores.deduce_datastore_type("ecoli-pe.prseq") == PairedReads{A}
        @test ReadDatastores.name(ds) == Symbol("ecoli-pe")
        @test ReadDatastores.max_read_length(ds) == 300
        @test ReadDatastores.orientation(ds) == FwRv
        @test check_show(ds, "Paired Read Datastore 'ecoli-pe': 20 reads (10 pairs)")
        @test firstindex(ds) == 1
        @test lastindex(ds) == 20
        @test Base.IteratorSize(ds) == Base.HasLength()
        @test Base.IteratorEltype(ds) == Base.HasEltype()
        @test Base.eltype(ds) == LongSequence{A}
        @test_throws BoundsError ds[100]
        @test_throws BoundsError buffer(ds)[100]
        @test_throws ErrorException buffer(ds, 1)
        @test_throws ErrorException buffer(ds, ReadDatastores._bytes_per_read(ds) - 1)
        @test_throws BoundsError load_sequence!(ds, 100, LongSequence{A}())
        
        
        
        @test collect(ds) == collect(buffer(ds)) == open(PairedReads{A}, "ecoli-pe.prseq") do ds
            collect(ds)
        end
        @test ds[5] == buffer(ds)[5] == load_sequence!(ds, 5, LongSequence{A}()) == load_sequence!(buffer(ds), 5, LongSequence{A}())
    end
    
    run_checks(DNAAlphabet{4})
    run_checks(DNAAlphabet{2})
    
    @test_throws ReadDatastores.DatastoreVersionError{PairedReads{DNAAlphabet{2}}} open(PairedReads{DNAAlphabet{2}}, "ecoli-paired-old.prseq")
    @test_throws ReadDatastores.DatastoreEncodingError{PairedReads{DNAAlphabet{2}}} open(PairedReads{DNAAlphabet{4}}, "ecoli-paired-old.prseq")
end