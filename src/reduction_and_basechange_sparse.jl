# reduction_and_basechange_sparse
function xor_sorted_vectors!(out::Vector{Int}, a::Vector{Int}, b::Vector{Int})
    empty!(out)
    i, j = 1, 1

    while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            i += 1
            j += 1
        elseif a[i] < b[j]
            push!(out, a[i])
            i += 1
        else
            push!(out, b[j])
            j += 1
        end
    end

    while i <= length(a)
        push!(out, a[i])
        i += 1
    end

    while j <= length(b)
        push!(out, b[j])
        j += 1
    end

    return out
end

"""
out = Int[]
println(xor_sorted_vectors!(out,[2,4,7], [4,5,7]))  # [2,5]
println(xor_sorted_vectors!(out,[1,3], [2,4]))      # [1,2,3,4]
println(xor_sorted_vectors!(out, [1,2,3], [1,2,3]))  # Int[]
"""

function low_of_sparse_column(col::Vector{Int})
    return isempty(col) ? -1 : col[end]
end

function get_boundary_columns_sparse(FSC)
    n = length(FSC[1])

    # cell/simplex -> index
    inj = Dict(zip([s[1] for s in FSC[1]], 1:n))

    cols = Vector{Vector{Int}}(undef, n)

    for j in 1:n
        cell = FSC[1][j][1]
        bd = boundary_operation(cell)

        if isempty(bd)
            cols[j] = Int[]
        else
            rows = [inj[s] for s in bd]
            sort!(rows)
            cols[j] = rows
        end
    end

    return cols
end

function reduction_and_basechange_sparse_from_columns(cols0)
    cols = [copy(c) for c in cols0]
    n = length(cols)
    m = maximum(vcat([0], [isempty(c) ? 0 : c[end] for c in cols]))

    col_change = Vector{Tuple{Int,Int}}()
    pivot_to_col = fill(0, m)
    scratch = Int[]
    sizehint!(scratch, m)

    for j in 1:n
        lj = low_of_sparse_column(cols[j])

        while lj != -1 && pivot_to_col[lj] != 0
            i = pivot_to_col[lj]

            xor_sorted_vectors!(scratch, cols[j], cols[i])
            empty!(cols[j])
            append!(cols[j], scratch)

            push!(col_change, (i, j))
            lj = low_of_sparse_column(cols[j])
        end

        if lj != -1
            pivot_to_col[lj] = j
        end
    end

    return cols, col_change
end

function reduce_basis_with_coords(colsA::Vector{Vector{Int}})
    n = length(colsA)
    Ared = [copy(c) for c in colsA]

    # coords[j] = Ared[j] が元の A のどの列の xor か
    coords = [Int[j] for j in 1:n]

    m = maximum(vcat([0], [isempty(c) ? 0 : c[end] for c in Ared]))
    pivot_to_col = fill(0, m)

    scratch_col = Int[]
    scratch_coord = Int[]

    for j in 1:n
        lj = low_of_sparse_column(Ared[j])

        while lj != -1 && pivot_to_col[lj] != 0
            i = pivot_to_col[lj]

            xor_sorted_vectors!(scratch_col, Ared[j], Ared[i])
            empty!(Ared[j])
            append!(Ared[j], scratch_col)

            xor_sorted_vectors!(scratch_coord, coords[j], coords[i])
            empty!(coords[j])
            append!(coords[j], scratch_coord)

            lj = low_of_sparse_column(Ared[j])
        end

        if lj == -1
            error("spaceDown columns are linearly dependent.")
        end

        pivot_to_col[lj] = j
    end

    return Ared, coords, pivot_to_col
end

