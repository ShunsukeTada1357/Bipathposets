
function left_records(intLwithB, dim, n, m)
    #RR = n + m + 2
    records = NamedTuple[]

    for int in intLwithB
        up_data = int[1]      # [up_interval, up_basis]
        down_data = int[2]    # [down_interval, down_basis]

        if cubical_basis_dimension(up_data[2]) == dim

            p = bipath_left_point(up_data[1], down_data[1], n, m)

            push!(records, (
                kind = :left,
                dim = dim,
                point = p,
                up_interval = up_data[1],
                down_interval = down_data[1],
                up_basis = up_data[2],
                down_basis = down_data[2]
            ))
        end
    end

    return records
end

function right_records(intRwithB, dim, n, m)
    records = NamedTuple[]

    for int in intRwithB
        up_data = int[1]
        down_data = int[2]

        if cubical_basis_dimension(up_data[2]) == dim
  
            p = bipath_right_point(up_data[1], down_data[1], n, m)
            push!(records, (
                kind = :right,
                dim = dim,
                point = p,
                up_interval = up_data[1],
                down_interval = down_data[1],
                up_basis = up_data[2],
                down_basis = down_data[2]
            ))
        end
    end

    return records
end

function up_records(sepa, dim, n, m)
    records = NamedTuple[]
    for a in 1:length(sepa[4][1])
        interval = sepa[4][1][a]
        basis = sepa[4][2][a]
        if cubical_basis_dimension(basis) == dim
 
            p = bipath_up_point(interval, n, m)
            push!(records, (
                kind = :up,
                dim = dim,
                point = p,
                up_interval = interval,
                down_interval = nothing,
                up_basis = basis,
                down_basis = nothing
            ))
        end
    end
    return records
end

function center_records(sepa, dim, n, m)
    RR = n + m + 2
    records = NamedTuple[]
    for a in 1:length(sepa[2][1])
        interval = sepa[2][1][a]
        basis = sepa[2][2][a]
        if cubical_basis_dimension(basis) == dim
            push!(records, (
                kind = :center,
                dim = dim,
                point = bipath_center_point(n, m),
                up_interval = interval,
                down_interval = interval,
                up_basis = basis,
                down_basis = basis
            ))
        end
    end
    return records
end

function down_records(sepb, dim, n, m)
    records = NamedTuple[]
    for a in 1:length(sepb[4][1])
        interval = sepb[4][1][a]
        basis = sepb[4][2][a]
        if cubical_basis_dimension(basis) == dim

            p = bipath_down_point(interval, n, m)
            push!(records, (
                kind = :down,
                dim = dim,
                point = p,
                up_interval = nothing,
                down_interval = interval,
                up_basis = nothing,
                down_basis = basis
            ))
        end
    end
    return records
end



function bipath_records_cubical(FSCa, FSCb)
    sepa, sepb ,dims,intLwithB ,intRwithB = _bipath_decomposition_raw(FSCa, FSCb, cubical_basis_dimension)
    n, m = FSCa[2] - 2, FSCb[2] - 2
    records = NamedTuple[]
    for dim in dims
        append!(records, left_records(intLwithB, dim, n, m))
        append!(records, right_records(intRwithB, dim, n, m))
        append!(records, up_records(sepa, dim, n, m))
        append!(records, center_records(sepa, dim, n, m))
        append!(records, down_records(sepb, dim, n, m))
    end
    return records
end
