#Contributed by Toshitaka Aoki
import Combinatorics as Comb
import SimpleGraphs as SG

#Get a bipath from R^2.
function get_rectangular_paths(init,ending,partition::Int) #init, ending ∈R^2
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

 

function clique_random(G,faces,path1,path2)##Contributed by Toshitaka Aoki
    SE1 = collect(G.E)  # G のエッジ（辺）のリスト
    n = length(G.V)     # 頂点数
    m = binomial(n,2)   # 完全グラフにおける最大エッジ数
    l = 2^n - 1         # 全てのfaceの数
    w1 = rand(m);
    sort!(w1)           #ランダムな m 個の数値を持つリストを作成し、ソート。SE1の重みに使う。
    w2 = Dict(zip(SE1,rand(m)));  #各エッジにランダムな重みを割り当てる辞書を作成。
    SE2 = sort(SE1, by = x -> w2[x]) #w2の値でエッジSE1をソート
    d1 = Vector(undef, l);
    d2 = Vector(undef, l);
    
    #初期状態（頂点単体のクリーク）
    for i in 1:n 
        d1[i] = [[i], 1] 
        d2[i] = [[i], 1] 
    end
 
    for i in 1:length(faces)
        f = faces[i]
        birth_f = [w1[findlast(e -> collect(e) ⊆ f, SE1)], w2[SE2[findlast(e -> collect(e) ⊆ f, SE2)] ]]#∈R^2
        d1[i+n] = [f, findfirst(p -> (birth_f[1] <= p[1]) && (birth_f[2] <= p[2]), path1)]
        d2[i+n] = [f, findfirst(p -> (birth_f[1] <= p[1]) && (birth_f[2] <= p[2]), path2)] 
        
        if i % 10000 == 0 
            println(i)
        end
    end

    d1 = d1[findall(x -> d1[x][2] != nothing, 1:l)]
    d2 = d2[findall(x -> d2[x][2] != nothing, 1:l)]
    return [d1,length(path1)], [d2,length(path2)] 
end  

function clique_random_ith(n,path1,path2,ith)
    G = SG.Complete(n) 
    faces = collect(Combinatorics.powerset(collect(G.V), 2, ith+2))
    l = length(faces)  
    SE1 = collect(G.E)  # G のエッジ（辺）のリスト
    n = length(G.V)     # 頂点数
    m = binomial(n,2)   # 完全グラフにおける最大エッジ数       
    w1 = rand(m);
    sort!(w1)           #ランダムな m 個の数値を持つリストを作成し、ソート。SE1の重みに使う。
    w2 = Dict(zip(SE1,rand(m)));  #各エッジにランダムな重みを割り当てる辞書を作成。
    SE2 = sort(SE1, by = x -> w2[x]) #w2の値でエッジSE1をソート
    d1 = Vector(undef, l+n);
    d2 = Vector(undef, l+n);
    
    #初期状態（頂点単体のクリーク）
    for i in 1:n 
        d1[i] = [[i], 1] 
        d2[i] = [[i], 1] 
    end
 
    for i in 1:length(faces)
        f = faces[i]
        birth_f = [w1[findlast(e -> collect(e) ⊆ f, SE1)], w2[SE2[findlast(e -> collect(e) ⊆ f, SE2)] ]]#∈R^2
        d1[i+n] = [f, findfirst(p -> (birth_f[1] <= p[1]) && (birth_f[2] <= p[2]), path1)]
        d2[i+n] = [f, findfirst(p -> (birth_f[1] <= p[1]) && (birth_f[2] <= p[2]), path2)] 
        
        if i % 10000 == 0 
            println(i)
        end
    end

    d1 = d1[findall(x -> d1[x][2] != nothing, 1:l)]
    d2 = d2[findall(x -> d2[x][2] != nothing, 1:l)]
    return [d1,length(path1)], [d2,length(path2)] 
end  
