#FSC (filtered simplicial complex)
using Combinatorics

function IsSC(SC)
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

function IsFSC(FSC)
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

    for i in 1:n
        SC=[FSC[1][j][1] for j in 1:n]
        if IsSC(SC) 
        else
            return IsSC(SC) 
        end
    end
    return true
end
