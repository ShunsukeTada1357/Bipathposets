using AbstractAlgebra
R = AbstractAlgebra.GF(2)
#############################
function standardbasis(dimension,field)
    #R = field
    Id = identity_matrix(R, dimension )
    listmat= [Id[:,i] for i in 1:dimension ]
    return listmat
end

function _vectorizationofFSC(FSC)
    sortedFSC= sort([s[1] for s in FSC[1]])
    n = length(sortedFSC)
    stdbasis=standardbasis(n,R)
    dic = Dict(zip(sortedFSC,stdbasis))
    return dic
end

function vectorizationofSC(FSC,sumofcomplex) #sumofcomplex=[[1,2],[2,3],[1,3]]
    FS = copy(FSC)
    dic=_vectorizationofFSC(FS)
    vec = [R(0) for i in 1:length(FSC[1])] 
    for s in sumofcomplex
        vec = vec +dic[s]
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

function reducematfromLtoR(matrix) #we assume the matrix is regular. We use it to connect up and down intervals.
    mat = copy(matrix)
    n = size(mat)[2]
    for i in 1:n-1
        l = find_lowest_nonzero(mat[:,i]) # l is not -1 since we assume mat is regular.
        for j in i+1:n
            if mat[l,j] != 0  ##use a field F_2
                mat[:,j] = mat[:,i] + mat[:,j]
            end
        end
    end 
    return mat
end

function boundary_operation(complex)#complex may be [1,2,3], or [1], and so on.
    n = length(complex)
    boundary=[]
    if n == 1
        return 0
    end
    for i in 1:n
      fix = copy(complex)  
       push!(boundary, deleteat!(fix,i))
    end
    return boundary
end
