#AtypeMethod.jl
#Using this code, we obtain intervals and i-th basis from FSC(filtered).
#This code generates a matrix encoding boundary operator of simplices. 
using OrderedCollections


# cols_sparse is a sparse reduced boundary representation (vector of support indices)
function interval_lt(I::Vector{Int}, J::Vector{Int})
    return I[1] < J[1] || (I[1] == J[1] && I[2] < J[2])
end

function get_intervals_with_representatives(FSC, cols_sparse, info_col_change)
    info = info_col_change
    n = length(FSC[1])

    basis = [[s[1]] for s in FSC[1]]
    newbasis = [[s[1]] for s in FSC[1]]

    # finite intervals
    imagebasis = Vector{Tuple{Any,Vector{Int}}}()

    # 全 intervals を (basis, interval) で持つ
    pair_list = Vector{Tuple{Any,Vector{Int}}}()

    # old E = collect(1:n) / setdiff! の代わり
    alive = trues(n)

    # finite intervals
    for i in 1:n
        j = low_of_sparse_column(cols_sparse[i])

        if j >= 0
            supp = cols_sparse[i]
            basisI = [basis[k][1] for k in supp]
            interval = [j, i - 1]

            push!(pair_list, (basisI, interval))
            push!(imagebasis, (basisI, interval))

            alive[i] = false
            alive[j] = false
        end
    end

    # 列変形情報を representative に反映
    for (src, dst) in info
        append!(newbasis[dst], newbasis[src])
    end

    # infinite intervals
    for i in 1:n
        if alive[i]
            push!(pair_list, (newbasis[i], [i, n]))
        end
    end

    # 区間で sort
    sort!(pair_list; by = last, lt = interval_lt)

    pairs = OrderedCollections.OrderedDict{Any,Any}()
    for (b, intv) in pair_list
        pairs[b] = intv
    end

    # imagebasis は旧版に合わせて [basis, interval] 形式
    imagebasis_oldstyle = [[b, intv] for (b, intv) in imagebasis]

    return [pairs, imagebasis_oldstyle]
end

##############################################################
#We contract birth and death of topological features [b,d] to some [b',d']. 
function contractbirth(birth,FSC)# I is an interval [b,d]
    return FSC[1][birth][2]
end
function contractdeath(death,FSC)
    if death == length(FSC[1])
        return FSC[2]
    else
        return FSC[1][death+1][2]-1
    end
end
function contractinterval(I, FSC)
    return [contractbirth(I[1],FSC), contractdeath(I[2],FSC)]
end
#################################################################

function baseswithintervals(FSC)
    fsc = [[FSC[1][i][1], i] for i in 1:length(FSC[1])]

    cols0 = get_boundary_columns_sparse([fsc, "any"])
    cols_red, info_col = reduction_and_basechange_sparse_from_columns(cols0)

    pairs, imagebasis = get_intervals_with_representatives(FSC, cols_red, info_col)

    newpairs = Dict()
    for k in keys(pairs)
        Interval = contractinterval(pairs[k], FSC)
        if Interval[1] <= Interval[2]
            newpairs[k] = Interval
        end
    end

pairvec = collect(newpairs)
sort!(pairvec; by = last, lt = interval_lt)

newpairs_sorted = OrderedCollections.OrderedDict{Any,Any}()
for p in pairvec
    newpairs_sorted[first(p)] = last(p)
end

newpairs = newpairs_sorted
    imagebasisleft = [a[1] for a in imagebasis if contractinterval(a[2], FSC)[2] == 0]
    imagebasis_only = [a[1] for a in imagebasis]

    return [newpairs, imagebasisleft, imagebasis_only]
end