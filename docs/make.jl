using FreePoleHeom
using Documenter

DocMeta.setdocmeta!(FreePoleHeom, :DocTestSetup, :(using FreePoleHeom); recursive=true)

makedocs(;
    modules=[FreePoleHeom],
    authors="Malte Krug <malte.krug@proton.me> and contributors",
    sitename="FreePoleHeom.jl",
    format=Documenter.HTML(;
        canonical="https://krugx.github.io/FreePoleHeom.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/krugx/FreePoleHeom.jl",
    devbranch="main",
)
