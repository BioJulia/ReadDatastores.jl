```@meta
CurrentModule = ReadDatastores
```

# PairedReadDatastores

`PairedReadDatastores` are designed to store paired sequences, such sequences
are generated from sequencing paired-end libraries or long mate pair libraries.

Typically such read datasets come in the form of two FASTQ files, often named
according to a convention `*_R1.fastq` and `*_R2.fastq`. One file contains all
the "left" sequences of each pair, and the other contains all the "right" sequences
of each pair. So the first read pair consistes of the first record in each file. 