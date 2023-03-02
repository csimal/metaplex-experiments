using DrWatson
@quickactivate "metaplex-experiments"

using Random
using BenchmarkTools
using JLD2

include(srcdir("HeterogeneousMetaplexExperiment.jl"))

Random.seed!(2023)

M = 10
N = 100 * M

ba_k = 5

h = barabasi_albert(M, ba_k)

β_base = 0.5
γ = 0.5
D = 0.1

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

du = copy(u0)

f!, p = meanfield_fun(mpx)

@benchmark f!(du, u0, p, 0.0)

mp = HeterogeneousMetapopulation(mpx)
du_mp = copy(u0_mp)

f_mp!, p_mp = meanfield_fun(mp)
@benchmark f_mp!(du_mp, u0_mp, p_mp, 0.0)

final_infection(mp, u0_mp, 10.0, solver=Rodas5())

final_infection(mpx, u0, 10.0, solver=Rosenbrock32())

@benchmark final_infection(mpx, u0, 10.0)

df1 = @time full_experiment(mpx, u0, β_factor, k_base, k_factor, tmax)
df2 = @time full_experiment(mpx, u0, β_factor, k_base, k_factor, tmax)


jldsave("benchmark.jld2", df1, df2)
