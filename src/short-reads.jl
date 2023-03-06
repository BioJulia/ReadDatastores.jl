abstract type ShortReads{A<:DNAAlphabet} <: ReadDatastore{LongSequence{A}} end

@inline BioSequences.bits_per_symbol(prds::ShortReads{A}) where {A<:DNAAlphabet} = BioSequences.bits_per_symbol(A())

@inline _offset_of_sequence(sr::ShortReads, idx::Integer) = _read_data_begin(sr) + (_bytes_per_read(sr) * (idx - 1))

@inline function inbounds_load_sequence!(ds::ShortReads{A}, i::Integer, seq::LongSequence{A}) where {A<:DNAAlphabet}
    pos = _offset_of_sequence(ds, i)
    seek(stream(ds), pos)
    seqsize = read(stream(ds), UInt64)
    resize!(seq, seqsize)
    return _load_sequence_data!(ds, seq)
end

@inline function load_sequence!(sr::ShortReads{A}, idx::Integer, seq::LongSequence{A}) where {A<:DNAAlphabet}
    checkbounds(sr, idx)
    return inbounds_load_sequence!(sr, idx, seq)
end

@inline function Base.getindex(sr::ShortReads{A}, idx::Integer) where {A<:DNAAlphabet}
    @boundscheck checkbounds(sr, idx)
    seq = eltype(sr)(undef, max_read_length(sr))
    return inbounds_load_sequence!(sr, idx, seq)
end

function Base.show(io::IO, sr::ShortReads)
    summary(io, sr)
end