#This file is for the function bipathpersistence(FSCa,FSCb).
using AbstractAlgebra
export bipathpersistence
R= AbstractAlgebra.GF(2)
function separateintervals(pairs,k)#k is FSC[2], #This code separate intervals into four types.  
    basis = collect(keys(pairs))
    intervals =[pairs[b] for b in basis ]
    n = length(intervals)
    leftintervals = []
    rightintervals = []
    center = []
    centerbasis = []
    others = []
    othersbasis = []
    leftbasis = []
    rightbasis = []
    for i in 1:n
        if intervals[i] ==[1,k] # [1,"∞"]
            push!(center,intervals[i]) 
            push!(centerbasis, basis[i])
        elseif intervals[i][1] == 1
            push!(leftintervals,intervals[i])
            push!(leftbasis, basis[i])
        elseif intervals[i][2] == k #"∞"
            push!(rightintervals,intervals[i])
            push!(rightbasis, basis[i])
        else 
            push!(others, intervals[i])
            push!(othersbasis, basis[i])
        end
    end
    return [leftintervals,leftbasis], [center,centerbasis], [rightintervals,rightbasis], [others,othersbasis] 
end

function get_basis_intervals_left(FSC,left_int_basis)   #The input "left_int_basis" is the result of the function separateintervals(-,-)[1].
    leftB = [vectorizationofSC(FSC,i) for i in left_int_basis[2]]
    if length(leftB) == 0
        return  matrix_space(R,length(FSC[1]),length(leftB))()#R=GF(2). Call the result Lsapce. 
    else
        return hcat(leftB...)
    end
end

function get_basis_intervals_right(FSC, right_int_basis)#The input "right_int_basis" is the result of the function separateintervals(-,-)[3].
    rightB = [vectorizationofSC(FSC,i) for i in right_int_basis[2]]
    if length(rightB) == 0
        return  matrix_space(R,length(FSC[1]),length(rightB))()#R=GF(2), call the result Lsapce 
    else
        return hcat(rightB...)
    end
end

function get_basis_intervals_left(FSC,left_int_basis,imagebasis_left) 
    leftB = [vectorizationofSC(FSC,i) for i in left_int_basis[2]]
    if length(leftB) == 0
         Lspace = matrix_space(R,length(FSC[1]), length(leftB))()#R=GF(2)
    else
        Lspace = hcat(leftB...)
    end 
    leftBB = [vectorizationofSC(FSC,i) for i in imagebasis_left]
    if length(leftBB) == 0
        LLspace = matrix_space(R,length(FSC[1]), length(leftBB))()#R=GF(2)
    else
        LLspace = hcat(leftBB...)
    end
    return [Lspace LLspace]
end

function get_basis_intervals_right(FSC,right_int_basis,imagebasis,centerwithB) 
    rightB = [vectorizationofSC(FSC,i) for i in right_int_basis[2]]
    if length(rightB) == 0
        Rspace=matrix_space(R,length(FSC[1]),length(rightB))()#R=GF(2)
    else 
        Rspace = hcat(rightB...)
    end

    rrightB=[vectorizationofSC(FSC,i[1]) for i in imagebasis]
    if length(rrightB) == 0
        RRspace=matrix_space(R,length(FSC[1]),length(rrightB) )()#R=GF(2)
    else
        RRspace= hcat(rrightB...)
    end

    rrrightB=[vectorizationofSC(FSC,i) for i in centerwithB[2]]
    if length(rrrightB)  == 0
        aspace=matrix_space(R,l,length(rrrightB) )()#R=GF(2)
    else
        aspace = hcat(rrrightB...)
    end 
    return [Rspace RRspace aspace]
end

function get_two_repmat(LspaceUp,RspaceUp,LspaceDown,RspaceDown)
    m = size(LspaceUp)[2]
    leftrepmat =solve(LspaceDown,LspaceUp)[1:m,1:m]
    n = size(RspaceUp)[2]
    rightrepmat =solve(RspaceDown,RspaceUp)[1:n,1:n]
    return leftrepmat, rightrepmat
end

function connect_updown(upintervalswithB,downintervalswithB,matrix) #matrix method!!
    conn = []
    n = length(upintervalswithB[1])# number of intervals
    redmat = reducematfromLtoR(matrix)
    for i in 1:n
        l = find_lowest_nonzero(redmat[:,i]) 
        push!(conn,[ [upintervalswithB[1][i],upintervalswithB[2][i] ],[downintervalswithB[1][l],downintervalswithB[2][l] ]])
    end
    return conn
end

function _bipathpersistence(FSCa,FSCb)
    pairsa = baseswithintervals(FSCa)[1]
    pairsb, imb_left, imb = baseswithintervals(FSCb)
    dima=Set([length(k[1])-1 for k in keys(pairsa)]) #k=[[4, 5], [5, 6], [4, 6]] can be a basis of 1-dimensional homology
    dimb=Set([length(k[1])-1 for k in keys(pairsb)]) 
    dims = sort(collect(union(dima,dimb)))

    a = separateintervals(pairsa,FSCa[2])
    b = separateintervals(pairsb,FSCb[2])

    RspaceDown = get_basis_intervals_right(FSCb,b[3],imb,b[2]) 
    LspaceDown = get_basis_intervals_left(FSCb,b[1],imb_left) 
    RspaceUp = get_basis_intervals_right(FSCa,a[3])
    LsapceUp = get_basis_intervals_left(FSCa,a[1])  
    mats= get_two_repmat(LsapceUp,RspaceUp,LspaceDown,RspaceDown)
    Lrepmat = mats[1]
    Rrepmat = mats[2]
    intL= connect_updown(a[1],b[1],Lrepmat)
    intR = connect_updown(a[3],b[3],Rrepmat)
    return intL,intR,a[4],a[2],b[4],dims
end

function bipathpersistence(FSCa,FSCb)
    X = _bipathpersistence(FSCa,FSCb)
    i_thhomology = Dict()
    for i in X[6]
       intL = [[int[1][1], int[2][1]] for int in X[1] if length(int[1][2][1])== i+1]
       intR = [[int[1][1], int[2][1]] for int in X[2] if length(int[1][2][1])== i+1]
       up= [X[3][1][a] for a in 1:length(X[3][1]) if length(X[3][2][a][1])==i+1 ]
       center= [X[4][1][a] for a in 1:length(X[4][1]) if length(X[4][2][a][1])==i+1 ]
       down = [X[5][1][a] for a in 1:length(X[5][1]) if length(X[5][2][a][1])==i+1]

       i_thhomology[i] =  [intL, intR, up,center,down]
       if iszero(i_thhomology[i])
          println( "¬  ∃ ",i,"_th homology" )
          println("................")
       else
       println( " ∃ ",i,"_th homology, ", "#[̂0,̂1] is ", length(center) )
       println("intervals with ̂0 : ",intL)
       println("intervals with ̂1 : ",intR)
       println("intervals up: ",up)
       println("intervals down: ",down)
       println("................")
       end
    end
    return i_thhomology, FSCa[2] ,FSCb[2]
end
