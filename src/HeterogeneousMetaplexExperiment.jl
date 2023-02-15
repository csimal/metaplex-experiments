using NetworkEpidemics
using Graphs
using DifferentialEquations
#using LSODA
using Parameters: @with_kw
using Statistics
using DataFrames


@with_kw struct HeterogeneousMetaplexExperiment
    seed::Int = 2023
    N::Int
    M::Int
    h::SimpleGraph
    gs::Vector{SimpleGraph}
    g_mod::SimpleGraph
    β::Float64 = 0.1
    γ::Float64 = 0.11
    D::Float64 = 0.1
    tmax::Float64 = 1000.0
    nmax::Int = 100_000_000
    nsims::Int = 1000
    nbins::Int = 200
end

function HeterogeneousMetaplexExperiment(h::SimpleGraph, gs, g_mod, dyn, D; kws...)
    return HeterogeneousMetaplexExperiment(
        M=nv(h),
        N=nv(g_mod),
        h=h,
        gs=gs,
        g_mod=g_mod,
        dynamics=dyn,
        D=D;
        kws...
    )
end

struct HeterogeneousMetaplexResult
    ode_us::Vector{Float64}
    ts::Vector{Float64}
end


function single_experiment(exp::HeterogeneousMetaplexExperiment, μ, u_0; kws...)
    gs = copy(exp.gs)
    gs[μ] = exp.g_mod
    dynamics =[SIS(exp.β, exp.γ) for _ in vertices(exp.h)]
    mpx = HeterogeneousMetaplex(exp.h, gs, exp.D, dynamics)
    fi = final_infection(mpx, u_0, exp.tmax; kws...)
    us = vector_array_flatten(fi.u)
    return us, fi.t
end

"""
    single_experiment(mpx, f!, p, u_0, u_0, tmax, μ; kws...)
"""
function single_experiment(mpx::HeterogeneousMetaplex, f!, p, u_0, tmax, μ; kws...)
    fi = final_infection(f!, p, u_0, tmax; kws...)
    M = nv(mpx.h)
    df = DataFrame(
        :h_hash => hash(mpx.h),
        :g_hash => hash(mpx.g),
        :h => mpx.h,
        :modified_node => μ, # 0 denotes no modified node
        :g_ks => (mean ∘ degree).(mpx.g),
        :beta => [d.β for d in mpx.dynamics],
        :gamma => mpx.dynamics[1].γ,
        :D => mpx.D[1],
        :final_infection_pattern => reshape(sum(fi[2,:,:,end], dims=1),M),
        :final_infection => sum(fi[2,:,:,end]),
        :t_final => fi.t[end],
    )
    return df
end

function single_experiment(mp::HeterogeneousMetapopulation, f!, p, u_0, tmax, μ; kws...)
    fi = final_infection(f!, p, u_0, tmax; kws...)
    
end

function full_experiment(mpx::HeterogeneousMetaplex, u_0, β_factor, k_base, k_factor, tmax; kws...)
    #mp = HeterogeneousMetapopulation(mpx)
    β = mpx.dynamics[1].β
    N = nv(mpx.g[1])
    f_mpx!, p_mpx = meanfield_fun(mpx)
    #f_mp!, p_mp = meanfield_fun(mp)
    p_μ = deepcopy(p_mpx)
    mpx_μ = deepcopy(mpx)
    df_mpx = DataFrame()
    #df_mp = DataFrame()
    df_ = single_experiment(mpx, f_mpx!, p_mpx, u_0, tmax, 0; kws...)
    append!(df_mpx, df_)
    #df_ = single_experiment(mp, f_mp!, p_mp, u0_mp, tmax, 0; kws...)
    for μ in vertices(mpx.h)
        g = erdos_renyi(N, (k_base*k_factor)/N)
        modify_system!(mpx_μ, p_μ, μ, β_factor*β, g)
        df_ = single_experiment(mpx_μ, f_mpx!, p_μ, u_0, tmax, μ; kws...)
        reset_system!(mpx_μ, p_μ, mpx, p_mpx, μ)
        append!(df, df_)
    end
    return df_mpx
end

"""
    modify_system(mp::HeterogeneousMetapopulation, μ, β, k)

Modify the system `mp` by changing the mean contact degree and transmission rate of metanode `μ`.

Return a new system with the updated parameters.
"""
function modify_system(mp::HeterogeneousMetapopulation{SIS}, μ, β, k)
    ks = copy(mp.ks)
    dynamics = copy(mp.dynamics)
    γ = dynamics[μ].γ
    ks[μ] = k
    dynamics[μ] = SIS(β,γ)
    return HeterogeneousMetapopulation(mp.h, mp.N, ks, mp.D, dynamics)
end

function modify_system(mpx::HeterogeneousMetaplex{SIS}, μ, β, g)
    gs = copy(mpx.g)
    dynamics = copy(mpx.dynamics)
    γ = dynamics[μ].γ
    gs[μ] = g
    dynamics[μ] = SIS(β,γ)
    return HeterogeneousMetaplex(mpx.h, gs, mpx.D, dynamics)
end

"""
    modify_system(p, mpx, μ, β, g))

Modify the ODE parameters of the mpx system
"""
function modify_system!(mpx::HeterogeneousMetaplex{SIS}, p, μ, β_new, g)
    mpx.g[μ] = g
    p[1][μ] = β_new
    p[4][μ] = float.(adjacency_matrix(g))
end

function reset_system!(mpx_μ, p_μ, mpx, p, μ)
    mpx_μ.g[μ] = mpx.g[μ]
    p_μ[1][μ] = p[1][μ]
    p_μ[4][μ] = p[4][μ]
end

"""
    final_infection(mpx, u_0, tmax; kws...)

Compute the solution of the system `mpx` with initial condition `u_0` up to time `tmax`.
"""
function final_infection(mpx, u_0, tmax; solver=Tsit5(), kws...)
    f!, p = meanfield_fun(mpx)
    return final_infection(f!, p, u_0, tmax; solver=solver, kws...)
end

function final_infection(f!, p, u_0, tmax; solver=Tsit5(), kws...)
    prob = ODEProblem(f!, u_0, (0.0,tmax), p)
    return solve(prob, solver; 
                 callback = TerminateSteadyState(), # Stop when near equilibrium
                 save_everystep = false, # only save last step
                 save_start = false, # don't save initial condition
                )
end

function solve_infection(mpx, u_0, tmax; solver=Tsit5(), kws...)
    f!, p = meanfield_fun(mpx)
    return solve_infection(f!, p, u_0, tmax; solver=solver, kws...)
end

function solve_infection(f!, p, u_0, tmax; solver=Tsit5(), kws...)
    prob = ODEProblem(f!, u_0, (0.0,tmax), p)
    return solve(prob, solver)
end

"""
    initial_condition(mpx, μ; uniform_seed = true, fraction_infected = 0.05)

Construct an initial condition for the system `mpx` that is a small perturbation around the disease-free state.

By default, the infected seed is uniformly spread around the metanodes. If `uniform_seed` is equal to `false`, then the infected seed is only located in metanode `μ`. `fraction_infected` controls the expected fraction of infected individuals in the seed.

Returns the initial condition for the Metaplex system, as well as the reduced metapopulation initial condition.
"""
function initial_condition(mpx::HeterogeneousMetaplex, μ; uniform_seed = true, fraction_infected = 0.05)
    M = nv(mpx.h)
    N = nv(mpx.g[1])
    u0 = zeros(num_states(mpx.dynamics), N, M)
    u0[1,:,:] .= 1/M
    if uniform_seed == false
        u0[2,:,μ] .= rand(N) .* fraction_infected
    else
        u0[2,:,:] .= rand(N,M) .* fraction_infected
    end
    u0[1,:,:] .-= u0[2,:,:]
    u0_mp = reshape(sum(u0, dims=2), (num_states(mpx.dynamics), M))
    return u0, u0_mp
end

function initial_condition(mp::HeterogeneousMetapopulation, μ; uniform_seed = true, fraction_infected = 0.05)
    M = nv(mp.h)
    N = mp.N
    u0 = zeros(num_states(mp.dynamics), M)
    u0[1,:] .= N/M
    if uniform_seed == false
        u0[2,μ] = rand() * (N/M) * fraction_infected
    else
        u0[2,:] .= rand(M) .* (N/M) * fraction_infected
    end
    u0[1,:] .-= u0[2,:]
    return u0 
end
