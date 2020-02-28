```@meta
CurrentModule = ReadDatastores
DocTestSetup = quote
    using BioSequences
    using ReadDatastores
    using FASTX
    println(pwd())
    fwq = open(FASTQ.Reader, "../test/ecoli_tester_R1.fastq")
    rvq = open(FASTQ.Reader, "../test/ecoli_tester_R2.fastq")
    PairedReads{DNAAlphabet{2}}(fwq, rvq, "ecoli-test-paired", "my-ecoli-test", 250, 300, 0, FwRv)
end
```

# Building Read Datastores

To build a read datastore you first need to decide what sort of read data you have.

The process of building these datastores is consistent, but for each datastore
there are datastore specific options the constructor will accept:

```@docs
PairedReads
LongReads
LinkedReads
```

# Loading Read Datastores

Once you have built a datastore, you can open it in other projects again using
a `Base.open` method:

```jldoctest
julia> ds = open(PairedReads{DNAAlphabet{2}}, "ecoli-test-paired.prseq", "my-ecoli-pe")
Paired Read Datastore 'my-ecoli-pe': 20 reads (10 pairs)
```

The open method takes a `ReadDatastore` type, the filename of the datastore,
and, optionally a name for the datastore, if you omit the name, the datastore
will use the default name that was specified on construction.

You can also use `do` block notation when opening a datastore and that will
ensure that the underlying stream of the reader will close.
