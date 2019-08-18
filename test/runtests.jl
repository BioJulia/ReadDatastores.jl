module TestReadDatastores

using Test, FASTX, ReadDatastores, BioSequences

include("paired-reads.jl")
include("long-reads.jl")
include("linked-reads.jl")

end # module
