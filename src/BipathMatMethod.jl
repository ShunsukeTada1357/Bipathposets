#This file is for the function bipathpersistence(FSCa,FSCb), where the input FSCa and FSCb are FSC.
using AbstractAlgebra
R= AbstractAlgebra.GF(2)

#This code separate intervals into four types.
function separateintervals(pairs,k)#k is FSC[2] 
    basis = collect(keys(pairs))
    intervals =[pairs[b] for b in basis ]
    leftintervals = []
    rightintervals = []
    center = []
    centerbasis = []
    others = []
    othersbasis = []
    leftbasis = []
    rightbasis = []
    for i in 1:length(intervals)
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

function _bipathpersistence(FSCa,FSCb)
    pairsa = baseswithintervals(FSCa)[1]
    pairsb, imb_left , imb = baseswithintervals(FSCb)
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
    intL= connect_updown(a[1],b[1],mats[1])
    intR = connect_updown(a[3],b[3],mats[2])
    return intL,intR,a[4],a[2],b[4],dims
end

function bipathpersistence(FSCa,FSCb)
    X = _bipathpersistence(FSCa,FSCb)
    dims =X[6]
    i_thhomology = Dict()
    for i in dims
       intL = [[int[1][1], int[2][1]] for int in X[1] if length(int[1][2][1])== i+1]
       intR = [[int[1][1], int[2][1]] for int in X[2] if length(int[1][2][1])== i+1]
       up= [X[3][1][a] for a in 1:length(X[3][1]) if length(X[3][2][a][1])==i+1 ]
       center=[X[4][1][x] for x in 1:length(X[4][1]) if length(X[4][2][x][1])==i+1 ]
       down =[ X[5][1][x] for x in 1:length(X[5][1]) if length(X[5][2][x][1])==i+1]
       i_thhomology[i] =  [intL, intR, up,center,down]
       if iszero(i_thhomology[i])
          println( "¬  ∃ ",i,"_th homology" )
          println("................")
       else
       println( " ∃ ",i,"_th homology, ", "#[̂0,̂1] is ", length(center) )
       print("intervals with ̂0: ")        #println("intervals with ̂0 : ",intL) 
       for int in intL
        s,t = int[2][2]-1, int[1][2]-1
          if s !=0 &  t != 0 
             print("<"*string(s)*","*string(t)*"> ")
          elseif (t == 0) & (s != 0 )
             print("<"*string(s)*"', ̂0> ")
          elseif (t!=0) & (s == 0) 
            print("<̂0,"*string(t)*"> ")
          else
            print("<̂0, ̂0>  ")
          end
       end
       println(" ")
       print("intervals with ̂1: ")        #println("intervals with ̂1 : ",intR) 
       for int in intR
        s,t = int[1][1]-1, int[2][1]-1
          if (int[1][1] != int[1][2]) &  (int[2][1] != int[2][2]) 
             print("<"*string(s)*","*string(t)*"'> ")
          elseif (int[1][1] == int[1][2])  & (int[2][1] != int[2][2]) 
             print("<̂1,"*string(t)*"'> ")
          elseif (int[1][1] != int[1][2])  & (int[2][1] == int[2][2]) 
            print("<̂"*string(s)*",̂1> ")
          else
            print("<̂1, ̂1> ")
          end
       end
       println(" ")
       print("intervals up: ")        #println("intervals up: ",up)
       for int in up
        print("<"*string(int[1]-1)*","*string(int[2]-1)*"> ")
       end
       println(" ")
       print("intervals down: ")    #println("intervals down: ",down)
       for int in down
        print("<"*string(int[1]-1)*"',"*string(int[2]-1)*"'> ")
       end
       println(" ")
       println("................")
       end       
    end
    return i_thhomology, FSCa[2] ,FSCb[2]
end
