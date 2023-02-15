using Random
using Graphs

function get_system_and_u0(rngseed,M,N,ba_k,β,γ,D,k,uniform_seed)
    rng = Xoshiro(rngseed)
    h = barabasi_albert(M,ba_k, rng=rng)
    mp = HeterogeneousMetapopulation(h,N,k,D,SIS(β,γ))
    u0 = initial_condition(mp, 0; rng=rng)
    return mp, u0
end

function makesim(d::Dict)
    @unpack rngseed, β, γ, D, beta_factor, k_factor, tmax, k, ba_k, u_seed = d
    mp, u0 = get_system_and_u0(rngseed, M, N, ba_k, β, γ, D, k, u_seed)
    fi = final_infection(mp, u0, tmax; solver=Rodas5())
    u_base = sum(fi[2,:,end]) / N
    t_base = fi.t[end]
    us = zeros(M)
    ts = zeros(M)
    for μ in 1:M
        mp_ = modify_system(mp, μ, β_factor*β, k_factor*k)
        fi = final_infection(mp_, u0, tmax; solver=Rodas5())
        us[μ] = sum(fi[2,:,end]) / N
        ts[μ] = fi.t[end]
    end
    d["u_base"] = u_base
    d["t_base"] = t_base
    d["us"] = us
    d["ts"] = ts
    d["h"] = mp.h
    wsave(datadir("simulations/metapopulation", savename(d, "jld2")), d)
    return true
end