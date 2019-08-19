```@meta
CurrentModule = ReadDatastores
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