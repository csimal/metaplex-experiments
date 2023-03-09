using Distributions

function random_geometric_network(N::Integer, d = 2, r)
    xs = rand(N,d)
    edges = [(i,j) for i in 1:N, j in i+1:N if norm(xs[i,:]-xs[j,:]) < r]
    return SimpleGraphFromIterator(edges)
end

function random_spatial_network(N::Integer; # number of nodes
    d::Integer = 2, # dimension
    Pk = Poisson(10), # Expected degree distribution
    Pf = Normal(0.0, 0.03) # Distance kernel
    )
    ks = rand(Pk, N)
    xs = rand(N,d)
    dists = Array{Float64,2}(undef, N, N)
    for i in 1:N, j in i:N
        dists[i,j] = norm(xs[i,:]-xs[j,:])
        dists[j,i] = dists[i,j]
    end
    f = x -> pdf(Pf, x)
    return random_spatial_network(dists, ks, f, 1), xs
end

function random_spatial_network(distances::AbstractMatrix, ks, f, ρ)
    N = size(distances, 1)
    k = mean(ks)
    ps = (ks * ks') .* f.(distances) / (ρ * k)
    return edge_probability_graph(ps)
end

function edge_probability_graph(ps::AbstractMatrix)
    N = size(ps, 1)
    edges = [(i,j) for i in 1:N-1, j in (i+1):N if rand() < ps[i,j]]
    return SimpleGraphFromIterator(edges)
end