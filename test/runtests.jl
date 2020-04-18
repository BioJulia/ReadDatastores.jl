module TestReadDatastores

using Test, FASTX, ReadDatastores, BioSequences

include("paired-reads.jl")
include("long-reads.jl")
include("linked-reads.jl")

@testset "Error Messages" begin
    buf = IOBuffer()
    
    Base.showerror(buf, ReadDatastores.MissingMagicError("myreads.prseq"))
    @test String(take!(buf)) == "MissingMagicError: the file myreads.prseq does not appear to be a valid read datastore file, it does not begin with the expected magic bytes."
    
    
end

end # module
