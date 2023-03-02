using DrWatson
@quickactivate "metaplex-experiments"

using Random
using Graphs
using JLD2

include(srcdir("HeterogeneousMetaplexExperiment.jl"))
include(srcdir("metaplex_sweep.jl"))

if !isempty(ARGS)
    seed = ARGS[1]
else
    seed = 2023
end

M = 100
N = 1000* M

ba_k = 5

β = 1.0
γ = 1.0
D = 1.0

k_base = 20

β_factor = 1
k_factor = 100

tmax = 1000.0

params = Dict(
    
)