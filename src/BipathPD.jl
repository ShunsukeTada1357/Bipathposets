##This file is for a function plotintlist(intlist,a,b)
##That visualize bipath intervals.     
using Plots
using LaTeXStrings
function intervalstocartesiancoordinate(intL,intR,up,down,center,n,m)
    intlist =[]
    for I in intL
        a = n+m+3-I[2][2] # an end point of I
        b = (I[1][2] + I[2][2]-2) # persistence of the interval I.
        push!(intlist,[a,(a+b)%(n+m+2)])
    end
    for I in intR
        a = I[1][1]-1
        b = (I[1][2]-I[1][1] ) + ( I[2][2] -I[2][1]) # persistence of the interval I.
        push!(intlist, [a, a+b])
    end
    for I in up
        a = I[1]-1
        b = I[2]-I[1]
        push!(intlist,[a, a+b])
    end
    for I in down
        a = n+m+3 -I[2]
        b = I[2]-I[1]
        push!(intlist,[a, a+b])
    end
    for I in center
        push!(intlist,[0,n+m+2])
    end
    return intlist
end

function countfunction(intlist)#[I,I,J,J,J]=> [[I,2],[J,3]]
    lis = copy(intlist)
    int_num = []
    while iszero(lis) != true 
        push!(int_num, [lis[1],length([lis[1] for J in lis if J == lis[1]]) ] ) 
        lis = [J for J in lis if J != lis[1]]
    end
    return int_num
end

#(5) color of points means the number of intervals C_{n,m}
function colorfunc(col,max::Int,mult::Int)
    n = length(col)
    d = Int64(floor(n/max))
    if d <= 1
        return 1
    elseif max == 1
        return 1
    else
        return min(1+(d)*(mult-1),n)
    end
end

function plotpoints(points,n,m)#points=[[a,b],[c,d],[a,b],...]
    plot()
    RR = n + m +2
    ints = countfunction(points)
    maxim = maximum([i[2] for i in ints])
    col =  cgrad(:turbo )
    num = length(col)
    d = div(num,maxim)

vline!([(n+1)%RR], linestyle=:dash, color=:black, label="1-hat ="*string(n+1))
hline!([(n+1)%RR], linestyle=:dash, color=:black, label="1-hat ="*string(m+1))
hline!([0], linestyle=:dash, color=:black, label="0-hat ="*string(m+1))
vline!([RR], linestyle=:dash, color=:black, label="0-hat ="*string(m+n+1))

lll=[]
for i in 1:n 
    push!(lll,string(i)) 
end
push!(lll,L"\widehat{1}")
for i in 1:m 
    push!(lll,string(m-i+1)*"'") 
end
uu = -1.5*ones(RR)
xaxe=copy(lll) 
push!(xaxe,L"\widehat{0}")
annotate!(1:RR, uu, text.(xaxe,10,"Computer Modern"))
yaxe= copy(lll)
insert!(yaxe,1,L"\widehat{0}")
annotate!(uu,0:RR-1, text.(yaxe,10,"Computer Modern"))

for i in 1:length(ints) ##plot points
    if ints[i][1][1] % RR  !=0 
       plot!([(ints[i][1][1])%RR], [(ints[i][1][2])%RR], 
       proj=:identity, ms=6,msw=:0.3,
       m=:o,
       c = col[colorfunc(col,maxim,ints[i][2])],#c = col[colorfunc(col,multiplicity[i]*d)],
       lt=:scatter,
       markerstrokealpha=1, legend = false,markeralpha=1,left_margin = 10Plots.mm, 
       bottom_margin = 10Plots.mm, aspect_ratio=:equal, label="",xlims=(-1,RR+0.5), ylims=(-1,RR+0.5),  
       fa=1, xticks=[], yticks=[])
    elseif ints[i][1][2] == RR
        plot!([0], [RR], 
        proj=:identity, ms=6,msw=:0.3,
        m=:o,
        c = col[colorfunc(col,maxim,ints[i][2])],
        lt=:scatter,
        markerstrokealpha=1,markeralpha=1,legend = false,left_margin = 10Plots.mm, 
        bottom_margin = 10Plots.mm, aspect_ratio=:equal,  label="",xlims=(-1,RR+0.5), ylims=(-1,RR+0.5),  
        fa=1, xticks=[], yticks=[])
    else
        plot!([RR], [(ints[i][1][2])], 
        proj=:identity, ms=6,msw=:0.3,
        m=:o,
        c = col[colorfunc(col,maxim,ints[i][2])],
        lt=:scatter,markeralpha=1, legend = false,left_margin = 10Plots.mm, 
        bottom_margin = 10Plots.mm  , label="",xlims=(-1,RR+0.5), ylims=(-1,RR+0.5),  
        fa=1, xticks=[], yticks=[],aspect_ratio=:equal)
    end
end ## End "plots points"
current()
end

#The first input is the result of a function "bipathpersistence" in BipathmatrixMthod.jl
function plotintlist(dict_resultof_bipathpersistence,i_th::Int) 
    X =dict_resultof_bipathpersistence
    if (i_th in keys(X[1])) == false #check that we have points to plot. 
        println("no ",i_th," homology.")
        return false
    end
    dict_int=X[1][i_th]
    n,m = X[2]-2, X[3]-2
    intL, intR, up, center, down = dict_int[1], dict_int[2], dict_int[3], dict_int[4], dict_int[5]
    intlist = intervalstocartesiancoordinate(intL,intR,up,down,center,n,m)
    plotpoints(intlist,n,m)
end
