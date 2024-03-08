using SimpleGraphs
using Random
using Combinatorics
#using .bipathposets


function comp(a, b)
    if (length(a[1]) < length(b[1])) 
        return true
    end
    return false
end

function reorderingFSC(FSC)
    days = FSC[2]
    want= []
    new=[]
    for i in 1:days
        new=[s for s in FSC[1] if s[2]== i ]
        new = sort(new, lt = comp)
        append!(want,new)
    end
    return [want,days]
end


function AttachProb(alledges)
    p = Dict()
    q = Dict()
    for edge in alledges
        p[edge] = rand()
       q[edge] = rand()   
    end
    w = p,q
    return w
end

##########################################

###########################################
function weightofface(face,w1,w2)
    edgesofface = collect(combinations(face,2))
    weightlis1 = [(w1[(sort(f)[1],sort(f)[2])]) for f in edgesofface] 
    weightlis2 = [(w2[(sort(f)[1],sort(f)[2])]) for f in edgesofface] 
    return findmax(weightlis1)[1], findmax(weightlis2)[1]
end
###########################################

###########################################
function weightoffaces(faces,w1,w2)
    wfaces = Dict()
    for face in faces
        if length(face) == 1 #not needed
        else
            wfaces[face] = weightofface(face,w1,w2)
        end
    end
    return wfaces
end


function get_wfaces(n)
    G = Complete(n)
    alledges = G.E
    vertices = G.V
    Pvertices = collect(vertices) 
    faces = collect(powerset(Pvertices,1,length(Pvertices)))
    w = AttachProb(alledges)
    w1=w[1]
    w2=w[2]
    wfaces = weightoffaces(faces,w1,w2)
    return wfaces
end
############################
function makeFClique(numberofvertices,weightoffaces,path)
    n = numberofvertices
    max = length(path)
    FCG=[ [[i],1] for i in 1:n]

    for face in keys(weightoffaces)
        if (weightoffaces[face][1]  <= path[1][1]) && (weightoffaces[face][2]  <= path[1][2])
            face = sort(face)    
            push!(FCG,[face,1])       
        end
    end 

    if max == 1
        return  reorderingFSC([FCG, max])
    end

    for i in 2:max
        for face in keys(weightoffaces)
            if (weightoffaces[face][1]  <= path[i][1]) && (weightoffaces[face][2]  <= path[i][2])
                if (weightoffaces[face][1]  > path[i-1][1]) || (weightoffaces[face][2]  > path[i-1][2])
                    push!(FCG,[sort(face),i])
                end
            end
        end
    end    
    FCG = reorderingFSC([FCG, max])
    return FCG
end


############################
function get_rectangular_paths(init,ending,partition::Int) #init in R^2
    underpath =[]
    uppath=[]
    lx = (ending[1] - init[1])/partition
    ly = (ending[2] - init[2])/partition
    for i in 1:partition
        push!(underpath,init+[lx*(i-1),0])   
        push!(uppath,init+[0,ly*(i-1)])  
    end
    push!(underpath,[ending[1],init[2]])
    push!(uppath,[init[1],ending[2]])
    for i in 1:partition-1
        push!(underpath,[ending[1],init[2]+ly*(i)])
        push!(uppath,[init[1]+lx*(i),ending[2]])
    end  
    push!(underpath,ending)
    push!(uppath,ending)

    return uppath, underpath
 
end

############################################################
function roof_function(r2,path)
    n = length(path)
    for i in 1:n
        if (r2[1] <= path[i][1]) && (r2[2] <= path[i][2])
        return i
        end
    end
    return false
end


function makebipathFClique(numberofvertices,weightoffaces,path1,path2)
    n = numberofvertices
    max1 = length(path1)
    max2 = length(path2)
    FCG1=[ [[i],1] for i in 1:n]
    FCG2 = copy(FCG1)

    for face in keys(weightoffaces)
        if (weightoffaces[face][1]  <= path1[1][1]) && (weightoffaces[face][2]  <= path1[1][2])    
            push!(FCG1,[sort(face),1])       
        end
        if (weightoffaces[face][1]  <= path2[1][1]) && (weightoffaces[face][2]  <= path2[1][2])
            push!(FCG2,[sort(face),1])       
        end
    end 

    M = maximum([max1,max2])
    for i in 2:M
        if max1 <= M
            for face in keys(weightoffaces)
               if (weightoffaces[face][1]  <= path1[i][1]) && (weightoffaces[face][2]  <= path1[i][2])
                   if (weightoffaces[face][1]  > path1[i-1][1]) || (weightoffaces[face][2]  > path1[i-1][2])
                       push!(FCG1,[sort(face),i])
                   end
               end
            end
        end

        if max2 <= M
            for face in keys(weightoffaces)
               if (weightoffaces[face][1]  <= path2[i][1]) && (weightoffaces[face][2]  <= path2[i][2])
                   if (weightoffaces[face][1]  > path2[i-1][1]) || (weightoffaces[face][2]  > path2[i-1][2])
                       push!(FCG2,[sort(face),i])
                   end
               end
            end
        end
    end    
    FCG1 = reorderingFSC([FCG1, max1])
    FCG2 = reorderingFSC([FCG2, max2])
    return FCG1, FCG2
end

##############################################

function comp(a, b)
    if (length(a[1]) < length(b[1])) 
        return true
    end
    return false
end


function make_clique_random(G,path1,path2)
    E = collect(G.E)
    SE = Set.(E)
    n= length(G.V)
    m = binomial(n,2)
    w1= rand(m)
    w2= rand(m)
    Δ1 = [] #Vector(undef, 2^n -1)
    Δ2 = [] #Vector(undef, 2^n -1)
    faces = collect(powerset(collect(G.V),2, n))
    faces = sort.(faces)
    #SE = Set.(E)
    t=1
    #E= sort.(E)
    for f in faces

        edges_in_f = Set.(collect(combinations(f,2)))
        birth_f = maximum.([[w1[findfirst(x-> x == y, SE)], w2[findfirst(x-> x == y, SE)]]  for y in edges_in_f] )
        
        if roof_function(birth_f,path1) != false
            push!(Δ1, [f, roof_function(birth_f,path1)])
        end
        if  roof_function(birth_f,path2) != false
            push!(Δ2,  [f, roof_function(birth_f,path2)])
        end
        if rem(t,200) == 0
           println(t,"/",length(faces))
        end
        t =t+1
    end 
    append!(Δ1, [[[i],1] for i in  1:n])
    append!(Δ2, [[[i],1] for i in  1:n])
    return reorderingFSC([Δ1,length(path1)]), reorderingFSC([Δ2,length(path2)]) 
end


