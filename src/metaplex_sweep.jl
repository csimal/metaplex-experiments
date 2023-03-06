using Random
using Graphs
using NetworkEpidemics

function get_system_and_u0(rngseed, M, N, ba_k, β, γ, D, k, u_seed)
    rng = Xoshiro(rngseed)
    h = barabasi_albert(M, ba_k, rng=rng)
    g = [erdos_renyi(N, k/N, rng=rng) for _ in 1:M]
    mpx = HeterogeneousMetaplex(h, g, D, SIS(β,γ))
    u0, u0_mp = initial_condition(mpx, 0; uniform_seed=u_seed, rng=rng)
    return mpx, u0, u0_mp
end

function makesim(d::Dict)
    @unpack rngseed, β, γ, D, M, N, beta_factor, k_factor, tmax, k, ba_k, u_seed = d
    mpx, u0 = get_system_and_u0(rngseed, M, N, ba_k, β, γ, D, k, u_seed)
    rng = Xoshiro(rngseed)
    fi = final_infection(mpx, u0, tmax; solver=Tsit5())
    u_base = sum(fi[2,:,:,end]) / N
    t_base = fi.t[end]
    us = zeros(M)
    ts = zeros(M)
    modified_degrees = zeros(M)
    for μ in 1:M
        g_μ = erdos_renyi(N, k*k_factor/N, rng=rng)
        modified_degrees[μ] = mean(degree(g_μ))
        mpx_ = modify_system(mpx, μ, beta_factor*β, g_μ)
        fi = final_infection(mpx_, u0, tmax; solver=Tsit5())
        us[μ] = sum(fi[2,:,:,end]) / N
        ts[μ] = fi.t[end]
    end
    d["u_base"] = u_base
    d["t_base"] = t_base
    d["base_degrees"] = (mean ∘ degree).(mpx.g)
    d["us"] = us
    d["ts"] = ts
    d["modified_degrees"] = modified_degrees
    d["h"] = h
    wsave(datadir("simulations/metaplex", savename(d, "jld2")), d)
end