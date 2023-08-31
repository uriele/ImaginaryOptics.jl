using ImaginaryOptics
using Documenter

DocMeta.setdocmeta!(ImaginaryOptics, :DocTestSetup, :(using ImaginaryOptics); recursive=true)

makedocs(;
    modules=[ImaginaryOptics],
    authors="Marco Menarini <menarini.marco@gmail.com>",
    repo="https://github.com/marcom/ImaginaryOptics.jl/blob/{commit}{path}#{line}",
    sitename="ImaginaryOptics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://marcom.github.io/ImaginaryOptics.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/marcom/ImaginaryOptics.jl",
    devbranch="master",
)
