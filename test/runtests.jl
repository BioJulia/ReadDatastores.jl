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
    @test String(take!(buf)) == "DatastoreTypeError: myreads.prseq contains a long read datastore and cannot be opened as a ReadDatastores.PairedReads{BioSequences.DNAAlphabet{2}}"
    
    Base.showerror(buf, ReadDatastores.DatastoreVersionError{PairedReads{DNAAlphabet{2}}}(UInt16(2)))
    @test String(take!(buf)) == "DatastoreVersionError: file format version of paired read datastore file (v2) is deprecated: this version of ReadDatastores.jl supports v$(Int(ReadDatastores.PairedDS_Version))"
    
    Base.showerror(buf, ReadDatastores.DatastoreVersionError{LongReads{DNAAlphabet{2}}}(UInt16(2)))
    @test String(take!(buf)) == "DatastoreVersionError: file format version of long read datastore file (v2) is deprecated: this version of ReadDatastores.jl supports v$(Int(ReadDatastores.LongDS_Version))"
    
    Base.showerror(buf, ReadDatastores.DatastoreVersionError{LinkedReads{DNAAlphabet{2}}}(UInt16(2)))
    @test String(take!(buf)) == "DatastoreVersionError: file format version of linked read datastore file (v2) is deprecated: this version of ReadDatastores.jl supports v$(Int(ReadDatastores.LinkedDS_Version))"
    
    Base.showerror(buf, ReadDatastores.DatastoreEncodingError{PairedReads{DNAAlphabet{2}}}("myreads.prseq", 4))
    @test String(take!(buf)) == "DatastoreEncodingError: myreads.prseq encodes reads using 4 bits per element and cannot be opened as a ReadDatastores.PairedReads{BioSequences.DNAAlphabet{2}}"
end

end # module
