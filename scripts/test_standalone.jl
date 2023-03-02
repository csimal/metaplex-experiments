using DrWatson
@quickactivate "metaplex-experiments"

using Random
using JLD2
using DifferentialEquations
using NetworkEpidemics
using LSODA

include(srcdir("HeterogeneousMetaplexExperiment.jl"))

Random.seed!(2023)

M = 10
N = 100

ba_k = 5

h = barabasi_albert(M, ba_k)

β_base = 1.0
γ = 1.05
D = 10.0

k_base = 100

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

f!, p = meanfield_fun(mpx)

u0, u0_mp = initial_condition(mpx, 1)

du = similar(u0)

f!(du, u0, p, 0.0)

function foo(f!, p, u0, solver)
    prob = ODEProblem(f!, u0, (0.0,10.0), p)
    return solve(prob, solver;
        #save_everystep = false,
        #save_start = false
    )
end

println("Starting simulation")
#@time fi = final_infection(mpx, u0, 10.0, solver=Rodas5())
@time sol = foo(f!, p, u0, Tsit5())
@time foo(f!, p, u0, lsoda())
println("Simulation complete!")
