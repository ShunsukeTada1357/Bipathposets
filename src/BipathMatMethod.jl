#This file is for the function bipathpersistence(FSCa,FSCb), where the input FSCa and FSCb are FSC.
using AbstractAlgebra
R= AbstractAlgebra.GF(2)

export bipathpersistence
#This code separate intervals into four types.
function separateintervals(pairs,k)#k is FSC[2]   
    basis = collect(keys(pairs))
    intervals =[pairs[b] for b in basis]
    
    n = length(intervals)
    leftintervals = []
    rightintervals = []
    center = []
    centerbasis = []
    others = []
    othersbasis =[]
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

# the input "left_int_basis" is the result of the function separateintervals(-,-)[1].
function get_basis_intervals_left(FSC,left_int_basis)   
    l=length(FSC[1])
    leftB = [vectorizationofSC(FSC,i) for i in left_int_basis[2]]
    cols = length(leftB)
    Lspace=matrix_space(R,l,cols)()#R=GF(2)
    for i in 1:cols
        Lspace[:,i]= leftB[i]
    end 
    return Lspace
end

# the input "right_int_basis" is the result of the function separateintervals(-,-)[3].
function get_basis_intervals_right(FSC, right_int_basis)# l=length(FSC[1])   
    l=length(FSC[1])
    rightB = [vectorizationofSC(FSC,i) for i in right_int_basis[2]]
    cols = length(rightB)
    Rspace=matrix_space(R,l,cols)()#R=GF(2)
    for i in 1:cols
        Rspace[:,i]= rightB[i]
    end 
    return Rspace
end

function get_basis_intervals_left(FSC,left_int_basis,imagebasis_left) 
    l=length(FSC[1])
    leftB = [vectorizationofSC(FSC,i) for i in left_int_basis[2]]
    cols1 = length(leftB)
    Lspace=matrix_space(R,l,cols1)()#R=GF(2)
    for i in 1:cols1
        Lspace[:,i]= leftB[i]
    end 
    leftBB=[vectorizationofSC(FSC,i) for i in imagebasis_left]

    cols2 = length(leftBB) 
    LLspace=matrix_space(R,l,cols2)()#R=GF(2)
    for i in 1:cols2
        LLspace[:,i]= leftBB[i]
    end 
    LLLspace=[Lspace LLspace]
    return LLLspace
end

function get_basis_intervals_right(FSC,right_int_basis,imagebasis,centerwithB) 
    l = length(FSC[1])
    rightB = [vectorizationofSC(FSC,i) for i in right_int_basis[2]]
    cols1 = length(rightB)
    Rspace=matrix_space(R,l,cols1)()#R=GF(2)
    for i in 1:cols1
        Rspace[:,i]= rightB[i]
    end 

    rrightB=[vectorizationofSC(FSC,i[1]) for i in imagebasis]
    cols2 = length(rrightB) 
    RRspace=matrix_space(R,l,cols2)()#R=GF(2)
    for i in 1:cols2
        RRspace[:,i]= rrightB[i]
    end 
    RRRspace=[Rspace RRspace]
    
    rrrightB=[vectorizationofSC(FSC,i) for i in centerwithB[2]]
    cols3 = length(rrrightB) 
    aspace=matrix_space(R,l,cols3)()#R=GF(2)
    for i in 1:cols3
        aspace[:,i]= rrrightB[i]
    end 
    RRRRspace=[RRRspace aspace]
    return RRRRspace
end

function get_two_repmat(LspaceUp,RspaceUp,LspaceDown,RspaceDown)
    m = size(LspaceUp)[2]
    leftrepmat =solve(LspaceDown,LspaceUp)[1:m,1:m]
    n = size(RspaceUp)[2]
    rightrepmat =solve(RspaceDown,RspaceUp)[1:n,1:n]
    return leftrepmat, rightrepmat
end
#matrixmethod!!
function connect_updown(upintervalswithB,downintervalswithB,matrix)
    conn = []
    n = length(upintervalswithB[1])# number of intervals
    redmat = copy(reducematfromLtoR(matrix))
    for i in 1:n
        l = find_lowest_nonzero(redmat[:,i]) 
        push!(conn,[ [upintervalswithB[1][i],upintervalswithB[2][i] ],[downintervalswithB[1][l],downintervalswithB[2][l] ]])
    end
    return conn
end

function _bipathpersistence(FSCa,FSCb)
    pairsa = baseswithintervals(FSCa)[1]
    pairsb, imb_left , imb = baseswithintervals(FSCb)
    dima=dim_int_list(pairsa)
    dimb=dim_int_list(pairsb)
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
####

function bipathpersistence(FSCa,FSCb)
    X = _bipathpersistence(FSCa,FSCb)
    dims =X[6]

    i_thhomology =Dict()
    for i in dims
       intL_old = X[1]
       intL = [[int[1][1], int[2][1]] for int in intL_old if length(int[1][2][1])== i+1]
       
       intR_old = X[2]
       intR = [[int[1][1], int[2][1]] for int in intR_old if length(int[1][2][1])== i+1]
    
       int_up_old=X[3][1]
       b_up_old=X[3][2]
       up= [int_up_old[x] for x in 1:length(int_up_old) if length(b_up_old[x][1])==i+1 ]
       up = [up[j] for j in findall(upper_than_diag,up)]
    
       int_center_old = X[4][1]
       b_center_old = X[4][2]
       center=[int_center_old[x] for x in 1:length(int_center_old) if length(b_center_old[x][1])==i+1 ]
    
       int_down_old=X[5][1]
       b_down_old=X[5][2]
       down =[ int_down_old[x] for x in 1:length(int_down_old) if length(b_down_old[x][1])==i+1]
       down = [down[j] for j in findall(upper_than_diag,down)]

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
