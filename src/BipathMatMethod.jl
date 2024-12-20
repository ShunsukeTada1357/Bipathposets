#This file is for the function interval_decomposition(FSCa,FSCb), where the input FSCa and FSCb are FSC.
using AbstractAlgebra
R= AbstractAlgebra.GF(2)

#This code separate intervals into four types.
function separateintervals(pairs, k) # k is FSC[2]
    basis = collect(keys(pairs))
    intervals = [pairs[b] for b in basis]    
    result = Dict(
        :left => ([], []),
        :center => ([], []),
        :right => ([], []),
        :others => ([], [])
    )
    
    for (i, interval) in enumerate(intervals)
        key = if interval == [1, k] # [1,"∞"]
            :center
        elseif interval[1] == 1
            :left
        elseif interval[2] == k
            :right
        else
            :others
        end
        push!(result[key][1], interval)
        push!(result[key][2], basis[i])
    end
    return result[:left], result[:center], result[:right], result[:others]
end

# the input "left_int_basis" is the result of the function separateintervals(-,-)[1].
function get_basis_intervals_left(FSC,left_int_basis)   
    leftB = [vectorizationofSC(FSC,i) for i in left_int_basis[2]]
    Lspace=matrix_space(R,length(FSC[1]),length(leftB))()#R=GF(2)
    for i in 1:length(leftB)
        Lspace[:,i]= leftB[i]
    end 
    return Lspace
end

# the input "right_int_basis" is the result of the function separateintervals(-,-)[3].
function get_basis_intervals_right(FSC, right_int_basis)# l=length(FSC[1])   
    rightB = [vectorizationofSC(FSC,i) for i in right_int_basis[2]]
    Rspace=matrix_space(R,length(FSC[1]),length(rightB))()#R=GF(2)
    for i in 1:length(rightB)
        Rspace[:,i]= rightB[i]
    end 
    return Rspace
end

function get_basis_intervals_left(FSC,left_int_basis,imagebasis_left) 
    l = length(FSC[1])
    leftB = [vectorizationofSC(FSC,i) for i in left_int_basis[2]]
    Lspace = matrix_space(R,l,length(leftB))()#R=GF(2)
    for i in 1:length(leftB)
        Lspace[:,i] = leftB[i]
    end 

    leftBB = [vectorizationofSC(FSC,i) for i in imagebasis_left]
    LLspace = matrix_space(R,l,length(leftBB) )()#R=GF(2)
    for i in 1:length(leftBB) 
        LLspace[:,i] = leftBB[i]
    end 
    return  [Lspace LLspace]
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
    
    rrrightB=[vectorizationofSC(FSC,i) for i in centerwithB[2]]
    cols3 = length(rrrightB) 
    aspace=matrix_space(R,l,cols3)()#R=GF(2)
    for i in 1:cols3
        aspace[:,i]= rrrightB[i]
    end 
    return [Rspace RRspace aspace]
end

function get_two_repmat(LspaceUp,RspaceUp,LspaceDown,RspaceDown)
    m = size(LspaceUp)[2]
    leftrepmat =solve(LspaceDown,LspaceUp,side = :right)[1:m,1:m]
    n = size(RspaceUp)[2]
    rightrepmat =solve(RspaceDown,RspaceUp,side = :right)[1:n,1:n]
    return leftrepmat, rightrepmat
end
#matrixmethod!!
function connect_updown(upintervalswithB,downintervalswithB,matrix)
    conn = []
    n = length(upintervalswithB[1])# number of intervals
    redmat = reducematfromLtoR(matrix)
    for i in 1:n
        l = find_lowest_nonzero(redmat[:,i]) 
        push!(conn,[ [upintervalswithB[1][i],upintervalswithB[2][i] ],[downintervalswithB[1][l],downintervalswithB[2][l] ]])
    end
    return conn
end

function interval_decomposition(FSCa,FSCb)
    pairsa = baseswithintervals(FSCa)[1]
    pairsb, imb_left , imb = baseswithintervals(FSCb)
    dima=Set([length(k[1])-1 for k in keys(pairsa)]) #k=[[4, 5], [5, 6], [4, 6]] can be a basis of 1-dimensional homology
    dimb=Set([length(k[1])-1 for k in keys(pairsb)]) 
    dims = sort(collect(union(dima,dimb)))           #We want to know the existence of non-zero qth homology modules. dims tells us this information.  

    sepa = separateintervals(pairsa,FSCa[2])
    sepb = separateintervals(pairsb,FSCb[2])

    RspaceDown = get_basis_intervals_right(FSCb,sepb[3],imb,sepb[2]) 
    LspaceDown = get_basis_intervals_left(FSCb,sepb[1],imb_left) 
    RspaceUp = get_basis_intervals_right(FSCa,sepa[3])
    LsapceUp = get_basis_intervals_left(FSCa,sepa[1])  
    mats= get_two_repmat(LsapceUp,RspaceUp,LspaceDown,RspaceDown)

    intLwithB= connect_updown(sepa[1],sepb[1],mats[1])
    intRwithB = connect_updown(sepa[3],sepb[3],mats[2])
    #return intLwithB, intRwithB, sepa[4], sepa[2], b[4], dims

    i_thhomology = Dict()
    for i in dims
       intL = [[int[1][1], int[2][1]] for int in intLwithB if length(int[1][2][1])== i+1]
       intR = [[int[1][1], int[2][1]] for int in intRwithB if length(int[1][2][1])== i+1]
       up= [sepa[4][1][a] for a in 1:length(sepa[4][1]) if length(sepa[4][2][a][1])==i+1 ]
       center=[sepa[2][1][x] for x in 1:length(sepa[2][1]) if length(sepa[2][2][x][1])==i+1 ]
       down =[ sepb[4][1][x] for x in 1:length(sepb[4][1]) if length(sepb[4][2][x][1])==i+1]

       i_thhomology[i] =  [intL, intR, up,center,down]

       if iszero(i_thhomology[i])
          println( "¬  ∃ ",i,"_th homology" )
          println("................")
       else
          println( " ∃ ",i,"_th homology, ", "#[̂0,̂1] is ", length(center) )
          print_intervals("intervals with ̂0: ", intL, format_intL)
          print_intervals("intervals with ̂1: ", intR, format_intR)
          print_intervals("intervals up: ", up, format_up)
          print_intervals("intervals down: ", down, format_down)
          println("................")
       end       
    end
    return i_thhomology, FSCa[2] ,FSCb[2]
end