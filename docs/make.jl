using WeirdMatrices
using Documenter

makedocs(;
    modules=[WeirdMatrices],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    repo="https://github.com/jagot/WeirdMatrices.jl/blob/{commit}{path}#L{line}",
    sitename="WeirdMatrices.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jagot.github.io/WeirdMatrices.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jagot/WeirdMatrices.jl",
)
