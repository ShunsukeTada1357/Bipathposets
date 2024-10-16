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

function get_intervals_with_basis(FSC, Bhat, info_col_change)
    info = info_col_change
    d = Dict()
        
    basis =  [[s[1]] for s in FSC[1]]    #equal to FS[1]
    imagebasis = []
    E = Vector(1: length(FSC[1]))
    for i in 1 : length(FSC[1])
        j = find_lowest_nonzero(Bhat[:,i])
        if j >= 0
            basisI = [basis[k][1] for k in 1:length(FSC[1]) if Bhat[k,i] !=0 ]
            interval = [j,i-1]
            d[basisI] = interval
            push!(imagebasis,[basisI,interval])
            setdiff!(E,[i,j])
        end 
    end

    newbasis = [[s[1]] for s in FSC[1] ] #equal to FSC[1]
    for vec in info
        for i in 1 : length(newbasis[vec[1]])
             push!(newbasis[vec[2]], newbasis[vec[1]][i])
        end
    end

    for i in E
        d[newbasis[i]] =[i,length(FSC[1])]
    end
    pairs = sort(d;byvalue=true)
    return  [pairs,imagebasis] #It returns pairs of intervals and bases. 
end

function baseswithintervals(FSC)
    fsc = [[FSC[1][i][1],i] for i in 1: length(FSC[1])]
    Bhat, info_col = getBhat_and_basechange( getB([fsc,"any"]))
    pairs,imagebasis = get_intervals_with_basis(FSC,Bhat,info_col)

    newpairs = Dict()
    for k in keys(pairs)
        Interval = contractinterval(pairs[k],FSC)
        if Interval[1]<= Interval[2]
            newpairs[k] = Interval
        end
    end   
    newpairs = sort(newpairs;byvalue=true)
    imagebasisleft = [a[1] for a in imagebasis if contractinterval(a[2],FSC)[2]==0] ## The left image of the vector space.  
    return [newpairs,imagebasisleft, imagebasis]                                       ## It returns pairs (intervals and bases) and basis of image of boundary operation.        
end