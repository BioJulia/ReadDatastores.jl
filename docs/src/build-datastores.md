```@meta
CurrentModule = ReadDatastores
DocTestSetup = quote
    using ReadDatastores
end
```

# Building Read Datastores

To build a read datastore you first need to decide what sort of read data you have.

If you have short (somewhat regular in length) paired reads such as that were
sequenced from a paired-end library or a long mate pair library, you should
build a `PairedReads` datastore.

If you have long reads of inconsistent length, such as you would get from
single-molecule real-time sequencing of a long read library, you should build
a `LongReads` datastore.

If you have short (somewhat regular in length) paired reads, and thr pairs are
linked linked through some proximity tagging mechanism (like 10x), then you
should build a LinkedReads datastore.

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
julia> ds = open(PairedReads, "test/ecoli-pe.prseq", "my-ecoli-pe")
Paired Read Datastore 'my-ecoli-pe': 20 reads (10 pairs)
```

The open method takes a `ReadDatastore` type, the filename of the datastore,
and, optionally a name for the datastore, which may be omitted, in which case
the datastore will use the default name that was specified on construction. 
