module ReadDatastores

export 
    PairedReadDatastore,
    PairedReadOrientation,
    FwRv,
    RvFw
    
using BioSequences, FASTX

const SeqDataStoreMAGIC = 0x05D5

@enum Filetype::UInt16 begin
    PairedDS = 1
end

###
### Writing helpers
###

@inline function writestring(fd::IO, s::AbstractString)
    write(fd, s) + write(fd, '\0')
end

include("paired-reads.jl")

end # module
