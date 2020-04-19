module TestReadDatastores

using Test, FASTX, ReadDatastores, BioSequences

include("paired-reads.jl")
include("long-reads.jl")
include("linked-reads.jl")

@testset "Error Messages" begin
    buf = IOBuffer()
    
    Base.showerror(buf, ReadDatastores.MissingMagicError("myreads.prseq"))
    @test String(take!(buf)) == "MissingMagicError: the file myreads.prseq does not appear to be a valid read datastore file, it does not begin with the expected magic bytes."
    
    Base.showerror(buf, ReadDatastores.DatastoreTypeError{PairedReads{DNAAlphabet{2}}}("myreads.prseq", ReadDatastores.LongDS))
    @test String(take!(buf)) == "DatastoreTypeError: myreads.prseq contains a long read datastore and cannot be opened as a PairedReads{DNAAlphabet{2}}"
end

end # module
