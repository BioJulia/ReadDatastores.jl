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

## Using the macros and string literals

If you try to open a datastore with the `open` method as above, and you provide
the wrong type you will get an error.

This will protect you from trying to open a long read datastore as a paired read
datastore and such, but it's not always convenient.

For example, if you have a paired read datastore storing reads in 2-bit format
and you tried to open it as a `PairedReads{DNAAlphabet{4}}` type you will still
get an error.

This is obviously correct behaviour, you don't want to be loading sequences
using a different encoding to the one they were saved with!

However, in a practical setting this will get annoying: maybe you want to use
some long reads you put in a datastore a while ago but don't remember if your
datastore file is a `LongReads{DNAAlphabet{4}}` or a `LongReads{DNAAlphabet{2}}`.
Or maybe you get a `somefile.prseq` file from a colleague, and from the extension,
you deduce it is paired reads but even then that's not guaranteed.

To solve this problem a few convenience macros are provided for you, so you can
load datastores without specifying the datastore type, yet still avoid type
uncertainty in your generated code.

The macro `@dsopen` accepts a filepath, and optionally a datastore name.
The macro is evaluated and expanded before you julia code is compiled.
During that time, the header of the datastore file is peeked at, and the correct
ReadDatastore subtype is deduced, and a correctly typed `open` statement is generated.

For example if the file `myreads.prseq` was a 2-bit encoded paired read datastore,
and you had the following statement in your script: `ds = @dsopen "myreads.prseq"`

The statement would be expanded to: `ds = open(PairedReads{DNAAlphabet{2}}, "myreads.prseq")`

You can also open a datastore using a string literal e.g. `reads"myreads.prseq"`.
When you do this, type type of the datastore is detected as with `@dsopen`,
however, rather than returning the expression
`ds = open(PairedReads{DNAAlphabet{2}}, "myreads.prseq")`, as `@dsopen` does,
the `open` is executed and the macro returns the value of `ds`.

!!! note
    In order for `@dsopen` to work properly in any script the datastore file
    must exist prior to running the script. This is unavoidable because
    macros are evaluated and expanded first before the resulting expanded code
    is compiled and run.
    
    So creating a datastore file and loading it again using `@dsopen` within the
    same script will not work, and `@dsopen` will try to peek at the file and
    deduce its contents before the script can generate the datastore file.
    You will get an error telling you the file can't be found / does not exist.
    
    In an interactive setting, in which statements are entered, compiled and run
    by the REPL one at a time, this should rarely be a problem.
    
**In summary: If in doubt about the datastore type of a file, simply use `@dsopen`**