#Contributed by Toshitaka Aoki
import Combinatorics as Comb
import SimpleGraphs as SG
using BenchmarkTools

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

function roof_function(r2,path)
    n = length(path)
    for i in 1:n
        if (r2[1] <= path[i][1]) && (r2[2] <= path[i][2])
        return i
        end
    end
    return false
end

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

function clique_random(G,faces,path1,path2)##writen by Toshitaka Aoki
    SE1 = collect(G.E)
    n = length(G.V)
    m = binomial(n,2)
    l = 2^n - 1 
    w1 = rand(m);
    sort!(w1) 
    w2 = Dict(zip(SE1,rand(m))); 
    SE2 = sort(SE1, by = x -> w2[x])
    d1 = Vector(undef, l);
    d2 = Vector(undef, l);
    for i in 1:n 
        d1[i] = [[i], 1] 
        d2[i] = [[i], 1] 
    end
 
    for i in 1:length(faces)
        f = faces[i]
        birth_f = [w1[findlast(e -> collect(e) ⊆ f, SE1)], w2[SE2[findlast(e -> collect(e) ⊆ f, SE2)]]]
        d1[i+n], d2[i+n] = [f, roof_function(birth_f,path1)], [f,roof_function(birth_f,path2)] 
        if i % 10000 == 0 
            println(i)
        end
    end

    d1 = d1[findall(x -> d1[x][2] != false, 1:l)]
    d2 = d2[findall(x -> d2[x][2] != false, 1:l)]
    return reorderingFSC([d1,length(path1)]), reorderingFSC([d2,length(path2)]) 
end  
