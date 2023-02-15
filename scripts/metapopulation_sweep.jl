using DrWatson
@quickactivate "metaplex-experiments"

using JLD2
using BenchmarkTools

include(srcdir("HeterogeneousMetaplexExperiment.jl"))
include(srcdir("metapopulation_sweep.jl"))

if !isempty(ARGS)
    seed = ARGS[1]
else
    seed = 2023
end


M = 100
N = 1000 * M

ba_k = 5

βs = 1.0#LinRange(0.0,20.0,200)
γs = 1.0
Ds = 1.0#LinRange(0.0,10.0,50)

k_base = 20

β_factor = 10
k_factor = 10

tmax = 1000.0

params = Dict(
    "rngseed" => 1:100,
    "β" => βs,
    "γ" => γs,
    "D" => Ds,
    "M" => M,
    "N" => N,
    "beta_factor" => β_factor,
    "k_factor" => 1,
    "tmax" => tmax,
    "k" => k_base,
    "ba_k" => ba_k,
    "u_seed" => "uniform",
)
dicts = dict_list(params)


@benchmark map(makesim, dicts)