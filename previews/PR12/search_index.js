var documenterSearchIndex = {"docs":
[{"location":"additional-methods/#Additional-methods-1","page":"Additional Methods","title":"Additional methods","text":"","category":"section"},{"location":"additional-methods/#","page":"Additional Methods","title":"Additional Methods","text":"All ReadDatastores permit construction from fastq files, iteration, and indexing.","category":"page"},{"location":"additional-methods/#","page":"Additional Methods","title":"Additional Methods","text":"But some specific ReadDatastore types also have additional specific methods. They are listed here.","category":"page"},{"location":"additional-methods/#","page":"Additional Methods","title":"Additional Methods","text":"read_tag\nname\nmax_read_length\norientation\nstream","category":"page"},{"location":"additional-methods/#ReadDatastores.read_tag","page":"Additional Methods","title":"ReadDatastores.read_tag","text":"read_tag(lr::LinkedReads, i::Integer)\n\nGet the tag for the i'th read\n\n\n\n\n\n","category":"function"},{"location":"additional-methods/#ReadDatastores.name","page":"Additional Methods","title":"ReadDatastores.name","text":"Get the name of the datastore\n\n\n\n\n\n","category":"function"},{"location":"additional-methods/#ReadDatastores.orientation","page":"Additional Methods","title":"ReadDatastores.orientation","text":"Get the orientation of the read pairs\n\n\n\n\n\n","category":"function"},{"location":"build-datastores/#","page":"Building & loading datastores","title":"Building & loading datastores","text":"CurrentModule = ReadDatastores\nDocTestSetup = quote\n    using BioSequences\n    using ReadDatastores\n    using FASTX\n    fwq = open(FASTQ.Reader, \"test/ecoli_tester_R1.fastq\")\n    rvq = open(FASTQ.Reader, \"test/ecoli_tester_R2.fastq\")\n    PairedReads{DNAAlphabet{2}}(fwq, rvq, \"test/ecoli-test-paired\", \"my-ecoli-test\", 250, 300, 0, FwRv)\nend","category":"page"},{"location":"build-datastores/#Building-Read-Datastores-1","page":"Building & loading datastores","title":"Building Read Datastores","text":"","category":"section"},{"location":"build-datastores/#","page":"Building & loading datastores","title":"Building & loading datastores","text":"To build a read datastore you first need to decide what sort of read data you have.","category":"page"},{"location":"build-datastores/#","page":"Building & loading datastores","title":"Building & loading datastores","text":"The process of building these datastores is consistent, but for each datastore there are datastore specific options the constructor will accept:","category":"page"},{"location":"build-datastores/#","page":"Building & loading datastores","title":"Building & loading datastores","text":"PairedReads\nLongReads\nLinkedReads","category":"page"},{"location":"build-datastores/#ReadDatastores.PairedReads","page":"Building & loading datastores","title":"ReadDatastores.PairedReads","text":"PairedReads{A}(rdrx::FASTQ.Reader, rdry::FASTQ.Reader, outfile::String, name::Union{String,Symbol}, minsize::Integer, maxsize::Integer, fragsize::Integer, orientation::PairedReadOrientation) where {A<:DNAAlphabet}\n\nConstruct a Paired Read Datastore from a pair of FASTQ file readers.\n\nPaired-end sequencing reads typically come in the form of two FASTQ files, often named according to a convention *_R1.fastq and *_R2.fastq. One file contains all the \"left\" sequences of each pair, and the other contains all the \"right\" sequences of each pair. The first read pair is made of the first record in each file.\n\nArguments\n\nrdrx::FASTQ.Reader: The reader of the *_R1.fastq file.\nrdxy::FASTQ.Reader: The reader of the *_R2.fastq file.\noutfile::String: A prefix for the datastore's filename, the full filename will include a \".prseq\" extension, which will be added automatically.\nname::String: A string denoting a default name for your datastore. Naming datastores is useful for downstream applications.\nminsize::Integer: A minimum read length (in base pairs). When building the datastore, if any pair of reads has one or both reads shorter than this cutoff, then the pair will be discarded.\nmaxsize::Integer: A maximum read length (in base pairs). When building the datastore, if any read has a greater length, it will be resized to this maximum length and added to the datastore.\nfragsize::Integer: The average fragment length of the paired end library that was sequenced. This value is entirely optional, but may be important for downstream applications.\norientation::PairedReadOrientation: The orientation of the reads. Set it to FwRv for building a datastore from a sequenced paired end library, and set it to RvFw if you are building the datastore from reads sequenced from a long mate pair library.\n\nExamples\n\njulia> using FASTX, ReadDatastores\n\njulia> fwq = open(FASTQ.Reader, \"test/ecoli_tester_R1.fastq\")\nFASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)\n\njulia> rvq = open(FASTQ.Reader, \"test/ecoli_tester_R2.fastq\")\nFASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)\n\njulia> ds = PairedReads{DNAAlphabet{2}}(fwq, rvq, \"ecoli-test-paired\", \"my-ecoli-test\", 250, 300, 0, FwRv)\n[ Info: Building paired read datastore from FASTQ files\n[ Info: Writing paired reads to datastore\n[ Info: Done writing paired read sequences to datastore\n[ Info: 0 read pairs were discarded due to a too short sequence\n[ Info: 14 reads were truncated to 300 base pairs\n[ Info: Created paired sequence datastore with 10 sequence pairs\nPaired Read Datastore 'my-ecoli-test': 20 reads (10 pairs)\n\n\n\n\n\n\n","category":"type"},{"location":"build-datastores/#ReadDatastores.LongReads","page":"Building & loading datastores","title":"ReadDatastores.LongReads","text":"LongReads{A}(rdr::FASTQ.Reader, outfile::String, name::String, min_size::Integer) where {A<:DNAAlphabet}\n\nConstruct a Long Read Datastore from a FASTQ file reader.\n\nArguments\n\nrdr::FASTQ.Reader: The reader of the fastq formatted file.\noutfile::String: A prefix for the datastore's filename, the full filename will include a \".loseq\" extension, which will be added automatically.\nname::String: A string denoting a default name for your datastore. Naming datastores is useful for downstream applications.\nmin_size::Integer: A minimum read length (in base pairs). When building the datastore, if any read sequence is shorter than this cutoff, then the read will be discarded.\n\nExamples\n\njulia> using FASTX, ReadDatastores\n\njulia> longrdr = open(FASTQ.Reader, \"test/human_nanopore_tester_2D.fastq\")\nFASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)\n\njulia> ds = LongReads{DNAAlphabet{2}}(longrdr, \"human-nanopore-tester\", \"nanopore-test\", 0)\n[ Info: Building long read datastore from FASTQ file\n[ Info: Writing long reads to datastore\n[ Info: Done writing paired read sequences to datastore\n[ Info: 0 reads were discarded due to a too short sequence\n[ Info: Writing index to datastore\n[ Info: Built long read datastore with 10 reads\nLong Read Datastore 'nanopore-test': 10 reads\n\n\n\n\n\n\n","category":"type"},{"location":"build-datastores/#ReadDatastores.LinkedReads","page":"Building & loading datastores","title":"ReadDatastores.LinkedReads","text":"LinkedReads{A}(fwq::FASTQ.Reader, rvq::FASTQ.Reader, outfile::String, name::String, format::LinkedReadsFormat, max_read_len::Integer, chunksize::Int = 1000000) where {A<:DNAAlphabet}\n\nConstruct a Linked Read Datastore from a pair of FASTQ file readers.\n\nPaired sequencing reads typically come in the form of two FASTQ files, often named according to a convention *_R1.fastq and *_R2.fastq. One file contains all the \"left\" sequences of each pair, and the other contains all the \"right\" sequences of each pair. The first read pair is made of the first record in each file.\n\nArguments\n\nrdrx::FASTQ.Reader: The reader of the *_R1.fastq file.\nrdxy::FASTQ.Reader: The reader of the *_R2.fastq file.\noutfile::String: A prefix for the datastore's filename, the full filename will include a \".prseq\" extension, which will be added automatically.\nname::String: A string denoting a default name for your datastore. Naming datastores is useful for downstream applications.\nformat::LinkedReadsFormat: Specify the linked reads format of your fastq files. If you have plain 10x files, set format to Raw10x. If you have 10x reads output in the UCDavis format, set the format to UCDavis10x.\nchunksize::Int = 1000000: How many read pairs to process per disk batch during the tag sorting step of construction.\n\nExamples\n\njulia> using FASTX, ReadDatastores\n\njulia> fqa = open(FASTQ.Reader, \"test/10x_tester_R1.fastq\")\nFASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)\n\njulia> fqb = open(FASTQ.Reader, \"test/10x_tester_R2.fastq\")\nFASTX.FASTQ.Reader{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(BioGenerics.Automa.State{TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}}(TranscodingStreams.TranscodingStream{TranscodingStreams.Noop,IOStream}(<mode=idle>), 1, 1, false), nothing)\n\njulia> ds = LinkedReads{DNAAlphabet{2s}}(fqa, fqb, \"10xtest\", \"ucdavis-test\", UCDavis10x, 250)\n[ Info: Building tag sorted chunks of 1000000 pairs\n[ Info: Dumping 83 tag-sorted read pairs to chunk 0\n[ Info: Dumped\n[ Info: Processed 83 read pairs so far\n[ Info: Finished building tag sorted chunks\n[ Info: Performing merge from disk\n[ Info: Leaving space for 83 read_tag entries\n[ Info: Chunk 0 is finished\n[ Info: Finished merge from disk\n[ Info: Writing 83 read tag entries\n[ Info: Created linked sequence datastore with 83 sequence pairs\nLinked Read Datastore 'ucdavis-test': 166 reads (83 pairs)\n\n\n\n\n\n\n","category":"type"},{"location":"build-datastores/#Loading-Read-Datastores-1","page":"Building & loading datastores","title":"Loading Read Datastores","text":"","category":"section"},{"location":"build-datastores/#","page":"Building & loading datastores","title":"Building & loading datastores","text":"Once you have built a datastore, you can open it in other projects again using a Base.open method:","category":"page"},{"location":"build-datastores/#","page":"Building & loading datastores","title":"Building & loading datastores","text":"julia> ds = open(PairedReads{DNAAlphabet{2}}, \"test/ecoli-test-paired.prseq\", \"my-ecoli-pe\")\nPaired Read Datastore 'my-ecoli-pe': 20 reads (10 pairs)","category":"page"},{"location":"build-datastores/#","page":"Building & loading datastores","title":"Building & loading datastores","text":"The open method takes a ReadDatastore type, the filename of the datastore, and, optionally a name for the datastore, which may be omitted, in which case the datastore will use the default name that was specified on construction.","category":"page"},{"location":"build-datastores/#","page":"Building & loading datastores","title":"Building & loading datastores","text":"You can also use do block notation when opening a datastore and that will ensure that the underlying stream of the reader will close.","category":"page"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"CurrentModule = ReadDatastores\nDocTestSetup = quote\n    using ReadDatastores, BioSequences\nend","category":"page"},{"location":"indexing/#Indexing-1","page":"Indexing & Iteration","title":"Indexing","text":"","category":"section"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"Base.getindex is defined for ReadDatastores:","category":"page"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"julia> ds = open(PairedReads{DNAAlphabet{2}}, \"test/ecoli-pe.prseq\", \"my-ecoli-pe\")\nPaired Read Datastore 'my-ecoli-pe': 20 reads (10 pairs)\n\njulia> ds[5]\n300nt DNA Sequence:\nNACATGCACTTCAACGGCATTACTGGTGACCTCTTCGTC…TTCTATCAACGCAAAAGGGTTACACAGATAATCGTCAGC","category":"page"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"Indexing a read datastore creates a new sequence. If you want to load a sequence from a datastore and into an existing sequence, then you can use the load_sequence! method.","category":"page"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"julia> seq = dna\"\"\n0nt DNA Sequence:\n< EMPTY SEQUENCE >\n\njulia> load_sequence!(ds, 6, seq)\n300nt DNA Sequence:\nATTACTGCGATTACTGCTGCGAATTTTTTCATGTTTATT…GTCCACTGGTTTACACAAGGTCGTAAGGGAAAAGAGGCG\n\njulia> seq\n300nt DNA Sequence:\nATTACTGCGATTACTGCTGCGAATTTTTTCATGTTTATT…GTCCACTGGTTTACACAAGGTCGTAAGGGAAAAGAGGCG\n","category":"page"},{"location":"indexing/#Iteration-1","page":"Indexing & Iteration","title":"Iteration","text":"","category":"section"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"The ReadDatastore types also support the Base.iterate interface:","category":"page"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"julia> collect(ds)\n20-element Array{LongSequence{DNAAlphabet{4}},1}:\n NGGGCTTTAAAATCCACTTTTTCCATATCGATAGTCACG…ATTTCTTCGATTCTTCTTTGTCACCGCAGCCAGCAAGAG\n GTGGGTTTTTATCGGCTGGCACATGTGTTGGGACAATTT…GGCTTTCAATACGCTGTTTTCCCTCGTTGTTTCATCTGT\n NTGAACTCCACATCCTGCGGATCGTAAACCGTCACCTCT…TCTTCCAGGCAGGCCGCCAGGGTATCACCTTCCAGACCA\n GATGAATCTGGCGGTTATTAACGGTAACAATAACCAGCA…AGACGGCAAACCGGCTGCAGGCGGTAGGTTGTTGCAGGT\n NACATGCACTTCAACGGCATTACTGGTGACCTCTTCGTC…TTCTATCAACGCAAAAGGGTTACACAGATAATCGTCAGC\n ATTACTGCGATTACTGCTGCGAATTTTTTCATGTTTATT…GTCCACTGGTTTACACAAGGTCGTAAGGGAAAAGAGGCG\n NCGGTTGAGTTCAAAGGCAAAGATTTGCTTGCGCTGTCG…TTTTCCGGCGGCGAGAAAAAGCGCAACGATTTTTTGCAA\n TTCGTCCCTGATATAGCACATGAACGTAATCAGGCTTGA…AATCTTCCGGCATCTTCAGGAGAGCGATTTTCTCTTCCA\n NACGACACATTACCGGAAATTCAGGCCGACCCGGACAGG…GTTGAACAACACGGTGGTACAATTCAGGTCGCAAGCCAG\n TCCACCACCAGAATATCGATATTATCGTGCGTCATCCTT…TCACGCCCGCGCCGCTTTCGCTGGCCGTCACGCTAATCA\n NCGTAACTTTATTCATATCTCTTCCCCCTCCCTGTACTT…TCTGTTACCGCATGGCGGCAGTGCGCTGGTCGATATGAC\n ATCGGGTAGGGGACGGAACGAATACGACAGTCAATATTC…AAGACTTTATCGTGCGGTCCGAACCGACTTTGTGGCGGC\n NGCCCTGGAGCTGGTGAAAGAAGGTCGAGCGCAAGCCTG…CAATCCTCGCGTGGCGTTGCTCAATATTGGTGAAGAAGA\n GAAAGGAACATCCTGACAACACCTTCCATCGTCTTTAAT…ATAAAGGCAAATTGCACCACCATGATGCTGTCCCAATCA\n NGTCTGGTGGTGCCTCTTTACTTAAGGAATTTCATCCTG…TAACGATGCCAGGCACCTGCGAAACTTTCCTGCACCAGC\n GACCGTTTTTCCCCAATCCGAGAACGCATAAATCCAGAC…TTTCTTCCCGGTAATGATACGTCACTATTGGAGTGGCCC\n NAGAGGCCACAGCGCGCCCATAATGGCGACTGAAAGCCA…TTTCACCGCGGTGACCGGAATCAGGGCAAATTCGACATG\n AAAAGGATCGCCGACCTTAACCATTCTGAATGTGATTGG…CTGGTGCCTGTCATATTTCGAACTCTGGGGGGACAGCAT\n NTGAGCAAATATGCCCGACCCAGCCTCATGACAGCGATA…AACCGAAAAAAAAGTAATCGTCGGCATGTCCGGCGGTGT\n AGGCTTTAAATTTGATCTCTTTGTTGCACAGAATATCCG…GCCAGGAAGAAACGGAGGAACCGACACCGCCGGCCATGC\n","category":"page"},{"location":"indexing/#Buffers-1","page":"Indexing & Iteration","title":"Buffers","text":"","category":"section"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"When iterating over a ReadDatastore either using Base.getindex or load_sequence!, you can sacrifice some memory for a buffer, to reduce the number of times the hard disk is read, and speed up sequential access and iteration. You can use the buffer method to wrap a datastore in such a buffer.","category":"page"},{"location":"indexing/#","page":"Indexing & Iteration","title":"Indexing & Iteration","text":"julia> bds = buffer(ds)\nBuffered Paired Read Datastore 'my-ecoli-pe': 20 reads (10 pairs)\n\njulia> for i in eachindex(bds)\n           load_sequence!(bds, i, seq)\n           println(length(seq))\n       end\n298\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n300\n","category":"page"},{"location":"#ReadDatastores-1","page":"Home","title":"ReadDatastores","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"(Image: Latest Release) (Image: MIT license)  (Image: Stable documentation) (Image: Latest documentation) (Image: Lifecycle) (Image: Chat)","category":"page"},{"location":"#Description-1","page":"Home","title":"Description","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Not your papa's FASTQ files.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"ReadDatastores provides a set of datastore types for storing and randomly accessing sequences from read datasets from disk. Each datastore type is optimised to the type of read data stored.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Using these data-stores grants greater performance than using text files that store reads (see FASTX.jl, XAM.jl, etc.) since the sequences are stored in BioSequences.jl succinct bit encodings already, and preset formats/layouts of the binary files means no need to constantly validate the input.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"A paired read datastore is provided for paired-end reads and long mate-pairs (Illumina MiSeq etc).\nA long read datastore is provided for long-reads (Nanopore, PacBio etc.)\nA linked read datastore is provided for shorter reads that are linked or grouped using some additional (typically proximity based) tag (10x).","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Also included is the ability to buffer these datastores, sacrificing some RAM, for faster iteration / sequential access of the reads in the datastore. ","category":"page"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"ReadDatastores is made available to install through BioJulia's package registry.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Julia by default only watches the \"General\" package registry, so before you start, you should add the BioJulia package registry.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Start a julia terminal, hit the ] key to enter pkg mode (you should see the prompt change from julia> to pkg> ), then enter the following command:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"registry add https://github.com/BioJulia/BioJuliaRegistry.git","category":"page"},{"location":"#","page":"Home","title":"Home","text":"After you've added the registry, you can install ReadDatastores from the julia REPL. Press ] to enter pkg mode again, and enter the following:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"add ReadDatastores","category":"page"},{"location":"#","page":"Home","title":"Home","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"},{"location":"#Testing-1","page":"Home","title":"Testing","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"ReadDatastores is tested against Julia 1.X on Linux, OS X, and Windows.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: ) (Image: )","category":"page"},{"location":"#Contributing-1","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Take a look at the contributing files detailed contributor and maintainer guidelines, and code of conduct.","category":"page"},{"location":"#Financial-contributions-1","page":"Home","title":"Financial contributions","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"We also welcome financial contributions in full transparency on our open collective. Anyone can file an expense. If the expense makes sense for the development of the community, it will be \"merged\" in the ledger of our open collective by the core contributors and the person who filed the expense will be reimbursed.","category":"page"},{"location":"#Backers-and-Sponsors-1","page":"Home","title":"Backers & Sponsors","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Thank you to all our backers and sponsors!","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Love our work and community? Become a backer.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: backers)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Does your company use BioJulia? Help keep BioJulia feature rich and healthy by sponsoring the project Your logo will show up here with a link to your website.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: )","category":"page"},{"location":"#Questions?-1","page":"Home","title":"Questions?","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"If you have a question about contributing or using BioJulia software, come on over and chat to us on Gitter, or you can try the Bio category of the Julia discourse site.","category":"page"},{"location":"read-datastores/#ReadDatastore-types-1","page":"Datastore types","title":"ReadDatastore types","text":"","category":"section"},{"location":"read-datastores/#","page":"Datastore types","title":"Datastore types","text":"ReadDatastores.jl implements several different types that index read sequence data on disk, permitting efficient random access and iteration, without having to keep every read in memory at once.","category":"page"},{"location":"read-datastores/#","page":"Datastore types","title":"Datastore types","text":"Currently there are three types of read datastore:","category":"page"},{"location":"read-datastores/#","page":"Datastore types","title":"Datastore types","text":"PairedReads handles short (somewhat regular in length) paired reads such as those sequenced from a paired-end library, or a long mate pair library.","category":"page"},{"location":"read-datastores/#","page":"Datastore types","title":"Datastore types","text":"LongReads handles long reads of inconsistent length, such as you would get from single-molecule real-time sequencing of a long read library.","category":"page"},{"location":"read-datastores/#","page":"Datastore types","title":"Datastore types","text":"LinkedReads handles short (somewhat regular in length) paired reads, however unlike typical paired-end reads pairs are also linked to each other, through some proximity tagging mechanism (like 10x).","category":"page"},{"location":"read-datastores/#","page":"Datastore types","title":"Datastore types","text":"All 3 types are subtypes of the ReadDatastore abstract type, and their API and behaviour is consistent, except in cases where the read data itself demands a divergence in behaviour (e.g. you can't ask for the tag of a read from a PairedReads dataset, but you can with a LinkedReads dataset).","category":"page"},{"location":"read-datastores/#","page":"Datastore types","title":"Datastore types","text":"The three types of read processed by the three types of datastore are not specific to any one technology, sequencing machine, or company. Rather the datastores  were chosen by considering the different characteristics of read datasets currently produced.","category":"page"}]
}
