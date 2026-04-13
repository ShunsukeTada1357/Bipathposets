#FSC (filtered simplicial complex)
using Combinatorics

function is_sc(SC)
    n = length(SC)
    for i in 1:n
        sub=combinations(SC[i])
        for p in sub
            if sort(p) in SC  
                #do nothing
            else
                println(p, "is not in SC")
                return false
            end    
        end
    end
    return true
end

function is_fsc(FSC)
    n = length(FSC[1])
    k =0
    for s in FSC[1] 
        if  s[2] < k
            println("order of birth is wrong")
            return false, s
        else
            k=s[2]
        end
    end

    println("order OK")

    if length([s[1] for s in FSC[1]]) != length(Set([s[1] for s in FSC[1]]))
        println("∃duble")
        return false 
    end

    println("¬ ∃ duble")

    
    SC=[FSC[1][j][1] for j in 1:n]
    if is_sc(SC) 
    else
        return is_sc(SC) 
    end
    
    return true
end


function is_cubical_fsc(FSC)
    isempty(FSC[1]) && return false
    cell = FSC[1][1][1]

    return cell isa Tuple && !isempty(cell) && cell[1] in (:v, :h, :w, :s)
end


function check_bipath_endpoint_condition(B_list_a, B_list_b)
    isempty(B_list_a) && error("B_list_a is empty.")
    isempty(B_list_b) && error("B_list_b is empty.")

    if B_list_a[1] != B_list_b[1]
        error("Invalid input for bipath API: B_list_a[1] must equal B_list_b[1].")
    end

    if B_list_a[end] != B_list_b[end]
        error("Invalid input for bipath API: B_list_a[end] must equal B_list_b[end].")
    end

    return true
end

