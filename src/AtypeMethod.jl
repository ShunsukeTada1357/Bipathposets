#Using this code, we obtain intervals and i-th basis from FSC(filtered).
#This code generates a matrix encoding boundary operator of simplices. 
using AbstractAlgebra
R = GF(2);
function getB(FSC)
    n = length(FSC[1])  
    B = zero_matrix(R,n,n)
    inj = Dict(zip([s[1] for s in FSC[1]], Vector(1:n) )) #injection: simplex -> index
    for i in 1:n
        if length(FSC[1][i][1]) != 1 # boundary_operation(FS[1][i][1])!=0
            S = boundary_operation(FSC[1][i][1])
            for s in S
                B[inj[s],i] = 1
            end
        end
    end
    return B
end

function getBhat_and_basechange(BB) #BB is a boundary matrix from a filtered simplicial complex over a field with two elements.
    Bhat = copy(BB)
    col_change=[]
    for j in 2 : size(BB)[2]
        i =1
        while i<j
            if  (find_lowest_nonzero(Bhat[:,i]) == find_lowest_nonzero(Bhat[:,j])) && find_lowest_nonzero(Bhat[:,j]) != -1
                Bhat[:,j] = Bhat[:,j]+Bhat[:,i]
                push!(col_change,[i,j])
                i= 1
            else
                i =i+1
            end 
        end
    end  
    return [Bhat,col_change]
end

function setF_Fdash_E(Bhat)
    n = size(Bhat)[2]
    F = []
    bijFdashtoF = Dict()
    Fdash =[]
    E =Vector(1:n)
    for i in 1:n
        a = find_lowest_nonzero(Bhat[:,i])
        if a >= 0
            push!(F,i)
            push!(Fdash,a)
            bijFdashtoF[a]=i
        end 
    end
    setdiff!(E,F)
    setdiff!(E,Fdash)
    return F, Fdash, E, bijFdashtoF
end

function get_intervals_with_basis(FSC,Bhat,info_col_change)
    FS = copy(FSC)
    info = info_col_change
    d = Dict()
    newbasis = [[s[1]] for s in FS[1] ]
    basis =  [[s[1]] for s in FS[1]]

    for vec in info
        for i in 1 : length(newbasis[vec[1]])
             push!(newbasis[vec[2]], newbasis[vec[1]][i])
        end
    end
 
    FFdashE = setF_Fdash_E(Bhat)
    Fdash = FFdashE[2]
    bijFdashtoF = FFdashE[4]
    E = FFdashE[3]
    imagebasis = []

    for k in Fdash
        j = bijFdashtoF[k]
        basisI = [basis[i][1] for i in 1:length(FS[1]) if Bhat[i,j] !=0 ]
        int = [k,j-1]
        d[basisI] = int
        push!(imagebasis,[basisI,int])
    end

    for i in E
        d[newbasis[i]] =[i,length(FS[1])]
    end
    pairs = sort(d;byvalue=true)
    return  [pairs,imagebasis] #It returns pairs of intervals and bases. 
end

function get_intervals_with_basis_contract(FSC,pairs)
    d = Dict()
    for k in keys(pairs)
        Int = contractinterval(pairs[k],FSC)
        if Int[1]<= Int[2]
            d[k] = contractinterval(pairs[k],FSC)
        end
    end    
    return d
end

function baseswithintervals(FSC)
    fsc = [[FSC[1][i][1],i] for i in 1: length(FSC[1])]
    Bhats = getBhat_and_basechange( getB([fsc,"any"]))
    info_col = Bhats[2]
    Bhat = Bhats[1]
    X = get_intervals_with_basis(FSC,Bhat,info_col)
    pairs = sort( get_intervals_with_basis_contract(FSC,X[1]);byvalue=true)
    imagebasisleft = [a[1] for a in X[2] if contractinterval(a[2],FSC)[2]==0] ## The left image of the vector space.  
    return [pairs,imagebasisleft, X[2]]                                       ## It returns pairs (intervals and bases) and basis of image of boundary operation.        
end
