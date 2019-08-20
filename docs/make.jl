using Documenter, ReadDatastores

makedocs(
    format = Documenter.HTML(),
    sitename = "ReadDatastores.jl",
    pages = [
        "Home"         => "index.md",
        "Building & loading datastores" => "build-datastores.md",
        "Indexing & Iteration" => "indexing.md"
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/ReadDatastores.jl.git",
    deps = nothing,
    make = nothing
)