using NetworkEpidemics
using Graphs
using DifferentialEquations
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

function compute_experiment(exp::HeterogeneousMetaplexExperiment)
    ts = LinRange(0.0, exp.tmax, exp.nbins)
    mpx = HeterogeneousMetaplex(exp.gs, exp.h, exp.D, exp.dynamics)

    for μ in vertices(exp.h)
        single_experiment(exp, μ, u_0)
    end
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

function single_experiment(mpx::HeterogeneousMetaplex, u_0, tmax, μ; kws...)
    fi = final_infection(mpx, u_0, tmax; kws...)
    M = nv(mpx.h)
    df = DataFrame(
        :h_hash => hash(mpx.h),
        :g_hash => hash(mpx.g),
        :h => mpx.h,
        :g_ks => (mean ∘ degree).(mpx.g),
        :beta => [d.β for d in mpx.dynamics],
        :gamma => mpx.dynamics[1].γ,
        :D => mpx.D[1],
        :final_infection => reshape(sum(fi[2,:,:,end], dims=1),M),
        :t_final => fi.t[end],
        :modified_node => μ # 0 denotes no modified node
    )
    return df
end

function full_experiment(mpx::HeterogeneousMetaplex, u_0, β_factor, k_base, k_factor, tmax; kws...)
    β = mpx.dynamics[1].β
    N = nv(mpx.g[1])
    df = DataFrame()
    df_ = single_experiment(mpx, u_0, tmax, 0; kws...)
    append!(df, df_)
    for μ in vertices(mpx.h)
        g = erdos_renyi(N, (k_base*k_factor)/N)
        mpx_ = modify_system(mpx, μ, β_factor*β, g)
        df_ = single_experiment(mpx_, u_0, tmax; kws...)
        append!(df, df_)
    end
    return df
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
    modify_system!(p, mpx, μ, β, g))

Modify the ODE parameters of the mpx system
"""
function modify_system!(p, mpx::HeterogeneousMetaplex{SIS}, μ, β, g)
    
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
    u0 = zeros(M, num_states(mp.dynamics))
    u0[:,1] .= N/M
    if uniform_seed == false
        u0[μ,2] = rand() * (N/M) * fraction_infected
    else
        u0[:,2] .= rand(M) .* (N/M) * fraction_infected
    end
    u0[:,1] .-= u0[:,2]
    return u0 
end
