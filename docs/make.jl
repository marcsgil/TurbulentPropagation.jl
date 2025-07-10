using Documenter, TurbulentPropagation

makedocs(
    sitename="TurbulentPropagation.jl",
    pages=[
        "index.md",
        "spectra.md",
    ],
    draft=true
)

deploydocs(
    repo="https://github.com/marcsgil/TurbulentPropagation.jl.git",
)