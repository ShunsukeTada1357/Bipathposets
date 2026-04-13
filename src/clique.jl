#Contributed by Toshitaka Aoki
import Combinatorics as Comb
import SimpleGraphs as SG

#Get a bipath from R^2.# clique.jl
# Contributed by Toshitaka Aoki
import Combinatorics as Comb
import SimpleGraphs as SG

# Get a bipath from R^2.
function get_rectangular_paths(init, ending, partition::Int) # init, ending ∈ R^2
    underpath = []
    uppath = []

    lx = (ending[1] - init[1]) / partition
    ly = (ending[2] - init[2]) / partition

    for i in 1:partition
        push!(underpath, init + [lx * (i - 1), 0])
        push!(uppath,   init + [0, ly * (i - 1)])
    end

    push!(underpath, [ending[1], init[2]])
    push!(uppath,    [init[1], ending[2]])

    for i in 1:partition-1
        push!(underpath, [ending[1], init[2] + ly * i])
        push!(uppath,    [init[1] + lx * i, ending[2]])
    end

    push!(underpath, ending)
    push!(uppath, ending)

    return uppath, underpath
end

function _clique_random_core(G, faces, path1, path2)
    edges = collect(G.E)
    nvtx = length(G.V)
    nedges = binomial(nvtx, 2)

    # vertices + higher-dimensional faces
    total_len = nvtx + length(faces)

    w1 = rand(nedges)
    sort!(w1)

    w2 = Dict(zip(edges, rand(nedges)))
    edges_sorted_by_w2 = sort(edges, by = e -> w2[e])

    d1 = Vector{Any}(undef, total_len)
    d2 = Vector{Any}(undef, total_len)

    # vertices
    for i in 1:nvtx
        d1[i] = [[i], 1]
        d2[i] = [[i], 1]
    end

    for (i, f) in enumerate(faces)
        idx1 = findlast(e -> collect(e) ⊆ f, edges)
        idx2 = findlast(e -> collect(e) ⊆ f, edges_sorted_by_w2)

        birth_f = [
            w1[idx1],
            w2[edges_sorted_by_w2[idx2]]
        ]

        d1[nvtx + i] = [f, findfirst(p -> (birth_f[1] <= p[1]) && (birth_f[2] <= p[2]), path1)]
        d2[nvtx + i] = [f, findfirst(p -> (birth_f[1] <= p[1]) && (birth_f[2] <= p[2]), path2)]

        if i % 10000 == 0
            println(i)
        end
    end

    valid_idx1 = findall(x -> x[2] !== nothing, d1)
    valid_idx2 = findall(x -> x[2] !== nothing, d2)

    d1 = d1[valid_idx1]
    d2 = d2[valid_idx2]

    return [d1, length(path1)], [d2, length(path2)]
end

function clique_random(G, faces, path1, path2)
    return _clique_random_core(G, faces, path1, path2)
end

function clique_random_ith(n, path1, path2, ith)
    G = SG.Complete(n)
    faces = collect(Comb.powerset(collect(G.V), 2, ith + 2))
    return _clique_random_core(G, faces, path1, path2)
end

