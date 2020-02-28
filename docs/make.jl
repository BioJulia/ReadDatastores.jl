using Documenter, ReadDatastores

makedocs(
    format = Documenter.HTML(),
    sitename = "ReadDatastores.jl",
    pages = [
        "Home"         => "index.md",
        "Datastore types" => "read-datastores.md",
        "Building & loading datastores" => "build-datastores.md",
        "Indexing & Iteration" => "indexing.md",
        "Additional Methods" => "additional-methods.md"
    ],
    authors = "Ben J. Ward, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/ReadDatastores.jl.git",
    deps = nothing,
    make = nothing,
    push_preview = true
)