using DrWatson
@quickactivate "metaplex-experiments"

using Plots
using Random
using BenchmarkTools
using JLD2

include(srcdir("HeterogeneousMetaplexExperiment.jl"))

Random.seed!(2023)

M = 10
N = 10 * M

ba_k = 5

h = barabasi_albert(M, ba_k)

β_base = 1.0
γ = 1.05
D = 1.0

k_base = 20

β_factor = 10
k_factor = 10

tmax = 1000.0

gs = [erdos_renyi(N,k_base/N) for _ in 1:M]
dynamics = [SIS(β_base,γ) for _ in 1:M]

mpx = HeterogeneousMetaplex(
    h,
    gs,
    [D, D],
    dynamics,
)

u0, u0_mp = initial_condition(mpx, 1)

#sol = final_infection(mpx, u0, 100.0)

b = @benchmark full_experiment(mpx, u0, β_factor, k_base, k_factor, tmax)

jldsave("benchmark.jld2", b)
