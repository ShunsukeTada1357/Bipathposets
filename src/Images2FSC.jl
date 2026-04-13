#images2FSC.jl
function make_or_filtration(A_list)
    n = length(A_list)
    B_list = Vector{Matrix{Int}}(undef, n)
    current = copy(Int.(A_list[1] .!= 0))
    B_list[1] = copy(current)

    for i in 2:n
        current = current .| Int.(A_list[i] .!= 0)
        B_list[i] = copy(current)
    end
    return B_list
end


function matrix_to_cubical_cells(B::AbstractMatrix{<:Integer})
    m, n = size(B)

    vertices = Set{Tuple}()
    hedges   = Set{Tuple}()
    vedges   = Set{Tuple}()
    squares  = Set{Tuple}()

    for x in 1:m, y in 1:n
        if B[x, y] == 1
            # 2-cell
            push!(squares, (:s, x, y))

            # その境界の4辺
            push!(hedges, (:h, x, y))
            push!(hedges, (:h, x, y + 1))
            push!(vedges, (:w, x, y))
            push!(vedges, (:w, x + 1, y))

            # その境界の4頂点
            push!(vertices, (:v, x, y))
            push!(vertices, (:v, x + 1, y))
            push!(vertices, (:v, x, y + 1))
            push!(vertices, (:v, x + 1, y + 1))
        end
    end

    return (
        vertices = collect(vertices),#リスト化
        hedges = collect(hedges),
        vedges = collect(vedges),
        squares = collect(squares)
    )
end

function cubical_filtration_from_matrices(B_list)
    cell_birth = Dict{Tuple, Int}()

    for i in 1:length(B_list)
        cells = matrix_to_cubical_cells(B_list[i])

        for c in cells.vertices
            if !haskey(cell_birth, c)#辞書cell_birth がkey c を持つ(has)たないか
                cell_birth[c] = i
            end
        end
        for c in cells.hedges
            if !haskey(cell_birth, c)
                cell_birth[c] = i
            end
        end
        for c in cells.vedges
            if !haskey(cell_birth, c)
                cell_birth[c] = i
            end
        end
        for c in cells.squares
            if !haskey(cell_birth, c)
                cell_birth[c] = i
            end
        end
    end

    return cell_birth
end


function cubical_filtration_to_FSC(B_list)
    cell_birth = cubical_filtration_from_matrices(B_list)

    cells = collect(keys(cell_birth))
    sort!(cells, by = c -> (cell_birth[c], cell_dimension(c)))

    fsc_list = [[c, cell_birth[c]] for c in cells]
    max_birth = isempty(cells) ? 0 : maximum(values(cell_birth))

    return [fsc_list, max_birth]
end
