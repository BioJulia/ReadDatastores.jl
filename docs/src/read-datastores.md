# `ReadDatastore` types

ReadDatastores.jl implements several different types that index read sequence
data on disk, permitting efficient random access and iteration, without having
to keep every read in memory at once.

Currently there are three types of read datastore:

`PairedReads` handles short (somewhat regular in length) paired reads such as
those sequenced from a paired-end library, or a long mate pair library.

`LongReads` handles long reads of inconsistent length, such as you would get from
single-molecule real-time sequencing of a long read library.

`LinkedReads` handles short (somewhat regular in length) paired reads,
however unlike typical paired-end reads pairs are also linked to each other,
through some proximity tagging mechanism (like 10x).

All 3 types are subtypes of the `ReadDatastore` abstract type, and their API and
behaviour is consistent, except in cases where the read data itself demands a
divergence in behaviour (e.g. you can't ask for the tag of a read from a `PairedReads`
dataset, but you can with a `LinkedReads` dataset).

The three types of read processed by the three types of datastore are not specific
to any one technology, sequencing machine, or company. Rather the datastores 
were chosen by considering the different characteristics of read datasets
currently produced.