using AbstractAlgebra
R = AbstractAlgebra.GF(2)

#############################
#Linear algebra, matrix op
function MakeElementary(n,i,j)
    M = identity_matrix(R, n)
    M[i,j] = 1
        return M
    end
    
function standardbasis(dimension,field)
    R = field
    n=dimension 
    Id = identity_matrix(R, n)
    listmat= [Id[:,i] for i in 1:n  ]
    return listmat
end

function _vectorizationofFSC(FSC)
    sortedFSC= sort([s[1] for s in FSC[1]])
    n = length(sortedFSC)
    stdbasis=standardbasis(n,R)
    dic = Dict()
    for i in 1:n
        dic[sortedFSC[i]] = stdbasis[i]
    end
    return dic
end


function vectorizationofSC(FSC,sumofcomplex) #sumofcomplex=[[1,2],[2,3],[1,3]]
    FS = copy(FSC)
    dic=_vectorizationofFSC(FS)
    vec=0
    n = length(sumofcomplex)
    for i in 1:n
        vec = vec +dic[sumofcomplex[i]]
    end
    return vec
end

function find_lowest_nonzero(vect) #changed name of func. lows=>find_lowest_nonzero
    n = length(vect)
    j = 1
     if iszero(vect)
        return -1
     else
        for i in 1:n
            if vect[i] != 0
                j = i
            end
        end
    end
    return j
end

function reducematfromLtoR(matrix) #we can assume the matrix is regular. We use it to connect up and down intervals.
    mat = copy(matrix)
    n = size(mat)[2]
    
    for i in 1:n-1
        l = find_lowest_nonzero(mat[:,i]) # l is not -1.
        for j in i+1:n
            if mat[l,j] != 0  ##use a field F_2
                mat[:,j]=mat[:,i]+mat[:,j]
            end
        end
    end 
    
    return mat
end

#####################################################
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


function subfiltFSC(FSC,i)
    sub =[s for s in FSC[1] if s[2]<= i]
    return [sub,i]
end

function IsFSC_weak(FSC)
    n = FSC[2]
    for i in 1:n
        sub = subfiltFSC(FSC,i)
       SC = [s[1] for s in sub[1]]
       if IsSC(SC) 
       else
           return false
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
        #for j in 1:i
        #    push!(SC,FSC[1][j][1])
        #end

        if IsSC(SC) 
        else
            return IsSC(SC) 
        end
    end
    return true
end

function FSCtoDict(FSC)
    FS = [s[1] for s in FSC[1]]
    index =[i for i in 1:length(FS)]
    return Dict(zip(FS,index))
end

function boundary_operation(complex)
    C=complex
    n = length(C)
    lis=[]
    if n == 1
        return 0
    end
    for i in 1:n
      fix = copy(C)  
       push!(lis, deleteat!(fix,i))
    end
    return lis
end
