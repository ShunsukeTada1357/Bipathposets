#BipathMatMethod.jl
#This file is for the function interval_decomposition(FSCa,FSCb), where the input FSCa and FSCb are FSC.

#This code separate intervals into four types.
function separateintervals(pairs, k) # k is FSC[2]
    basis = collect(keys(pairs))
    intervals = [pairs[b] for b in basis]    
    separation = Dict(
        :left => ([], []),
        :center => ([], []),
        :right => ([], []),
        :others => ([], [])
    )
    
    for (i, interval) in enumerate(intervals)
        key = if interval == [1, k] # [1,"∞"]
            :center
        elseif interval[1] == 1
            :left
        elseif interval[2] == k
            :right
        else
            :others
        end
        push!(separation[key][1], interval)
        push!(separation[key][2], basis[i])
    end
    return separation[:left], separation[:center], separation[:right], separation[:others]
end


function make_space_sparse(FSC_index::Dict, int_basis)
    cols = Vector{Vector{Int}}(undef, length(int_basis))
    for j in 1:length(int_basis)
        cols[j] = vectorizationofSC_support(FSC_index, int_basis[j])
    end
    return cols
end

function get_repmat_sparse_from_cols(
    colsUp::Vector{Vector{Int}},
    colsDown::Vector{Vector{Int}},
    mkeep::Int,
)
    Ared, coords, pivot_to_col = reduce_basis_with_coords(colsDown)

    coeff_cols = Vector{Vector{Int}}(undef, length(colsUp))

    scratch_b = Int[]
    scratch_x = Int[]
    scratch_trunc = Int[]

    for j in 1:length(colsUp)
        b = copy(colsUp[j])
        x = Int[]

        lb = low_of_sparse_column(b)
        while lb != -1
            i = (lb <= length(pivot_to_col)) ? pivot_to_col[lb] : 0
            if i == 0
                error("Column $j of spaceUp is not in the span of spaceDown.")
            end

            xor_sorted_vectors!(scratch_b, b, Ared[i])
            empty!(b)
            append!(b, scratch_b)

            xor_sorted_vectors!(scratch_x, x, coords[i])
            empty!(x)
            append!(x, scratch_x)

            lb = low_of_sparse_column(b)
        end

        # 旧 get_repmat(...)[1:m,1:m] に対応:
        # down 側の最初の mkeep 本だけ残す
        empty!(scratch_trunc)
        for t in x
            if t <= mkeep
                push!(scratch_trunc, t)
            end
        end
        coeff_cols[j] = copy(scratch_trunc)
    end

    return coeff_cols
end


function connect_updown_sparse(upintervalswithB, downintervalswithB, coeff_cols)
    conn = []
    n = length(upintervalswithB[1])

    cols_red, _ = reduction_and_basechange_sparse_from_columns(coeff_cols)

    for i in 1:n
        l = low_of_sparse_column(cols_red[i])
        if l == -1
            error("Column $i became zero in connect_updown_sparse.")
        end

        push!(conn, [
            [upintervalswithB[1][i], upintervalswithB[2][i]],
            [downintervalswithB[1][l], downintervalswithB[2][l]]
        ])
    end

    return conn
end

########
# Function to print intervals for intL
function print_intL(interval)
    s, t = interval[2][2] - 1, interval[1][2] - 1
    return "<" * (s == 0 ? "̂0" : string(s)) * "' ," * (t == 0 ? "̂0" : string(t)) * "> "
end

# Function to print intervals for intR
function print_intR(interval)
    s, t = interval[1][1] - 1, interval[2][1] - 1
    return "<" * (interval[1][1] == interval[1][2] ? "̂1" : string(s) )* "," * (interval[2][1] == interval[2][2] ?  "̂1" : string(t)*"'") * "> "
end
# Function to print intervals for up
function print_up(interval)
    return "<" * string(interval[1] - 1) * "," * string(interval[2] - 1) * "> "
end
# Function to print intervals for down
function print_down(interval)
    return "<" * string(interval[1] - 1) * "', " * string(interval[2] - 1) * "'> "
end
function print_intervals(header, intervals, printter)
    print(header)
    for interval in intervals
        print(printter(interval))
    end
    println(" ")
end
################################


function _bipath_decomposition_raw(FSCa, FSCb, basis_dim_fn)
    pairsa = baseswithintervals(FSCa)[1]
    pairsb, imb_left, imb = baseswithintervals(FSCb)

    dima = Set([basis_dim_fn(k) for k in keys(pairsa)])
    dimb = Set([basis_dim_fn(k) for k in keys(pairsb)])
    dims = sort(collect(union(dima, dimb)))
   
    sepa = separateintervals(pairsa, FSCa[2])
    sepb = separateintervals(pairsb, FSCb[2])
    FSCa_index = vectorizationofFSC_index(FSCa)
    FSCb_index = vectorizationofFSC_index(FSCb)
    right_basis_down = append!(copy(sepb[3][2]), imb, sepb[2][2])
    left_basis_down  = append!(copy(sepb[1][2]), imb_left)

    RcolsUp   = make_space_sparse(FSCa_index, sepa[3][2])
    LcolsUp   = make_space_sparse(FSCa_index, sepa[1][2])
    RcolsDown = make_space_sparse(FSCb_index, right_basis_down)
    LcolsDown = make_space_sparse(FSCb_index, left_basis_down)

    Lcoeff = get_repmat_sparse_from_cols(LcolsUp, LcolsDown, length(sepb[1][2]))
    Rcoeff = get_repmat_sparse_from_cols(RcolsUp, RcolsDown, length(sepb[3][2]))

    intLwithB = connect_updown_sparse(sepa[1], sepb[1], Lcoeff)
    intRwithB = connect_updown_sparse(sepa[3], sepb[3], Rcoeff)
    return sepa, sepb ,dims,intLwithB ,intRwithB     
end

function _interval_decomposition_common(FSCa, FSCb, basis_dim_fn; verbose=true)
    sepa, sepb ,dims,intLwithB ,intRwithB = _bipath_decomposition_raw(FSCa, FSCb, basis_dim_fn)
    
    i_thhomology = Dict()
    for i in dims
        intL = [[int[1][1], int[2][1]] for int in intLwithB if basis_dim_fn(int[1][2]) == i]
        intR = [[int[1][1], int[2][1]] for int in intRwithB if basis_dim_fn(int[1][2]) == i]
        up = [sepa[4][1][a] for a in 1:length(sepa[4][1]) if basis_dim_fn(sepa[4][2][a]) == i]
        center = [sepa[2][1][x] for x in 1:length(sepa[2][1]) if basis_dim_fn(sepa[2][2][x]) == i]
        down = [sepb[4][1][x] for x in 1:length(sepb[4][1]) if basis_dim_fn(sepb[4][2][x]) == i]

        i_thhomology[i] = [intL, intR, up, center, down]
        if verbose == true
            if isempty(intL) && isempty(intR) && isempty(up) && isempty(center) && isempty(down) #iszero(i_thhomology[i])
                println("∄ ", i, "_th homology")
                println("................")
            else
                println(" ∃ ", i, "_th homology, ", "#[̂0,̂1] is ", length(center))
                print_intervals("intervals with ̂0: ", intL, print_intL)
                print_intervals("intervals with ̂1: ", intR, print_intR)
                print_intervals("intervals up: ", up, print_up)
                print_intervals("intervals down: ", down, print_down)
                println("................")            
            end
        
        end
    end

    return i_thhomology, FSCa[2], FSCb[2]
end


"""
    interval_decomposition(FSCa, FSCb; verbose=true)

Compute the interval decomposition of bipath PH obtained from two filtered simplicial complexes

It returns the decomposition grouped by homological degree.

"""

function interval_decomposition(FSCa, FSCb; verbose=true)
    if is_cubical_fsc(FSCa) || is_cubical_fsc(FSCb)
        error("interval_decomposition is not for cubical filtrations. Use interval_decomposition_cubical instead.")
    end

    return _interval_decomposition_common(FSCa, FSCb, simplicial_basis_dimension; verbose=verbose)
end



"""
    interval_decomposition_cubical(FSCa, FSCb; verbose=true)

Compute the interval decomposition of bipath PH obtained from two cubical simplicial complexes

This is the main entry point for bipath interval decomposition in the cubical setting,
especially for data coming from binary images or image filtrations.

"""

function interval_decomposition_cubical(FSCa, FSCb; verbose=true)
    if !is_cubical_fsc(FSCa) || !is_cubical_fsc(FSCb)
        error("interval_decomposition_cubical is for cubical filtrations. Use interval_decomposition instead.")
    end

    return _interval_decomposition_common(FSCa, FSCb, cubical_basis_dimension; verbose=verbose)
end
