#FSCUtils.jl
function vectorizationofFSC_index(FSC)
    sortedFSC = sort([s[1] for s in FSC[1]])
    return Dict(sortedFSC[i] => i for i in 1:length(sortedFSC))
end

function vectorizationofSC_support(FSC_index::Dict, sumofcomplex)
    supp = Int[]
    scratch = Int[]

    for s in sumofcomplex
        idx = FSC_index[s]
        xor_sorted_vectors!(scratch, supp, [idx])
        empty!(supp)
        append!(supp, scratch)
    end

    return supp
end

function boundary_operation(simplex::AbstractVector)
    n = length(simplex)
    if n == 1
        return typeof(simplex)[]
    end
    return [deleteat!(copy(simplex), i) for i in 1:n]
end

function boundary_operation(cell::Tuple)
    tag = cell[1]

    if tag == :v
        return Tuple[]
    elseif tag == :h
        _, x, y = cell
        return [(:v, x, y), (:v, x + 1, y)]
    elseif tag == :w
        _, x, y = cell
        return [(:v, x, y), (:v, x, y + 1)]
    elseif tag == :s
        _, x, y = cell
        return [
            (:h, x, y),
            (:h, x, y + 1),
            (:w, x, y),
            (:w, x + 1, y)
        ]
    else
        error("Unknown cell type: $cell")
    end
end

function cell_dimension(cell)
    tag = cell[1]
    if tag == :v
        return 0
    elseif tag == :h || tag == :w
        return 1
    elseif tag == :s
        return 2
    else
        error("Unknown cell type: $cell")
    end
end

simplicial_basis_dimension(basis) = length(basis[1]) - 1
cubical_basis_dimension(basis) = cell_dimension(basis[1])
