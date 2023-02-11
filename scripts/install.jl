using Pkg
using DrWatson

@quickactivate "metaplex-experiments"

packages = [
    "BenchmarkTools",
    "DataFrames",
    "DifferentialEquations",
    "DrWatson",
    "Graphs",
    "JLD2",
    "Plots",
    "Random",
    "Statistics",
]

for pkg in packages
    #Pkg.add(pkg)
end

Pkg.add(url="https://github.com/csimal/NetworkEpidemics.jl.git")
