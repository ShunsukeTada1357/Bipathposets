##This file is for a function plotintlist(intlist,a,b)
##That visualize bipath intervals.     

using Plots
using LaTeXStrings

function BfindendpointL(I,n,m)
    return n+m+3-I[2][2]
end

function BfindendpointR(I)
    return I[1][1]-1
end

function BfindendpointUp(I)
    return I[1]-1
end

function BfindendpointDown(I,n,m)
    return n+m+3 -I[2]
end

function persistenceL(I)
    return (I[1][2] + I[2][2]-2)
end

function persistenceR(I)
    return (I[1][2]-I[1][1] ) + ( I[2][2] -I[2][1])
end

function persistenceUp(I)
    return I[2]-I[1]
end

function persistenceDown(I)
    return I[2]-I[1]
end

#(3) intervals to cartesian coordinate
function intervaltocartesiancoordinateL(I,n,m)
    return [BfindendpointL(I,n,m), (persistenceL(I)+BfindendpointL(I,n,m))%(n+m+2)]
end

function intervaltocartesiancoordinateR(I)
    return [BfindendpointR(I), persistenceR(I)+BfindendpointR(I)]
end

function intervaltocartesiancoordinateUp(I)
    return [BfindendpointUp(I), persistenceUp(I)+BfindendpointUp(I)]
end

function intervaltocartesiancoordinateDown(I,n,m)
    return [BfindendpointDown(I,n,m), persistenceDown(I)+BfindendpointDown(I,n,m)]
end

function intervalstocartesiancoordinate(intL,intR,up,down,center,n,m)
    intlist =[]
    for I in intL
        push!(intlist,intervaltocartesiancoordinateL(I,n,m))
    end

    for I in intR
        push!(intlist,intervaltocartesiancoordinateR(I))
    end

    for I in up
        push!(intlist,intervaltocartesiancoordinateUp(I))
    end

    for I in down
        push!(intlist,intervaltocartesiancoordinateDown(I,n,m))
    end

    for I in center
        push!(intlist,[0,n+m+2])
    end

    return intlist
end

function intlisttointlistwithmultiplicity(intlist)#[I,I,J,J,J]=> [[I,2],[J,3]]
    lis=copy(intlist)
    int_num=[]
    while iszero(lis) != true 
        push!(int_num, [lis[1],length([lis[1] for J in lis if J ==lis[1]]) ] ) 
        lis = [J for J in lis if J!= lis[1]]
    end
    return int_num
end

function intervalstocartesiancoordinatewithmulti(intlist)
    newintlistwith = intlisttointlistwithmultiplicity(intlist)
    xlist = []
    ylist = []
    multi = []
    for I in newintlistwith
        push!(xlist, I[1][1])
        push!(ylist, I[1][2])
        push!(multi, I[2])
    end
    return xlist,ylist,multi
end
#(5) color of points means the number of intervals C_{n,m}
function colorfunc(col,r::Int)
    n = length(col)
    if r <= 1
        return 1
    elseif r>= n
        return n
    else 
        return r
    end
end

function colorfunc2(col,max::Int,mult::Int)
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
#The first input is the result of a function "bipathpersistence" in BipathmatrixMthod.jl


function plotintlist(dict_resultof_bipathpersistence,i_th::Int)#####a=0,b=0
    ###check that we have points to plot. 
    X =dict_resultof_bipathpersistence
    dict_int=X[1][i_th]
    if (i_th in keys(X[1])) == false
        println("no ",i_th," homology.")
        return false
    else
        if iszero(dict_int[1]) &&  iszero(dict_int[2]) && iszero(dict_int[3]) && iszero(dict_int[5]) 
            println("all the intervals are bipath posets")
            return false 
        end
    end
###End "check that we have points to plot" 
###plots points 
    n = X[2]-2
    m = X[3]-2
    RR = n + m +2
        
    intL = dict_int[1]
    intR= dict_int[2]
    up =  dict_int[3]
    center = dict_int[4]
    down= dict_int[5]
    intlist = intervalstocartesiancoordinate(intL,intR,up,down,center,n,m)
    plot()
    x_ax, y_ax, multiplicity = intervalstocartesiancoordinatewithmulti(intlist)
    maxim = maximum(multiplicity)
    col =  cgrad(:turbo )
    num = length(col)
    d = div(num,maxim)

## draw labels on a plane
a=0
b=0
vline!([(n+1+a)%RR], linestyle=:dash, color=:black, label="1-hat ="*string(n+1))
hline!([(n+1+b)%RR], linestyle=:dash, color=:black, label="1-hat ="*string(m+1))
hline!([b%RR], linestyle=:dash, color=:black, label="0-hat ="*string(m+1))
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
## End "draw labels on a plane"

##plot points
for i in 1:length(x_ax)
    if x_ax[i] % RR  !=0 
       plot!([(x_ax[i]+a)%RR], [(y_ax[i]+b)%RR], 
       proj=:identity, ms=6,msw=:0.3,
       m=:o,
       c = col[colorfunc2(col,maxim,multiplicity[i])],#c = col[colorfunc(col,multiplicity[i]*d)],
       lt=:scatter,
       markerstrokealpha=1, legend = false,markeralpha=1,left_margin = 10Plots.mm, 
       bottom_margin = 10Plots.mm, aspect_ratio=:equal, label="",xlims=(-1,RR+0.5), ylims=(-1,RR+0.5),  
       fa=1, xticks=[], yticks=[])
    elseif y_ax[i] == RR
        plot!([0], [RR], 
        proj=:identity, ms=6,msw=:0.3,
        m=:o,
        c = col[colorfunc2(col,maxim,multiplicity[i])],
        lt=:scatter,
        markerstrokealpha=1,markeralpha=1,legend = false,left_margin = 10Plots.mm, 
        bottom_margin = 10Plots.mm, aspect_ratio=:equal,  label="",xlims=(-1,RR+0.5), ylims=(-1,RR+0.5),  
        fa=1, xticks=[], yticks=[])
    else
        plot!([RR], [y_ax[i]], 
        proj=:identity, ms=6,msw=:0.3,
        m=:o,
        c = col[colorfunc2(col,maxim,multiplicity[i])],
        lt=:scatter,markeralpha=1, legend = false,left_margin = 10Plots.mm, 
        bottom_margin = 10Plots.mm  , label="",xlims=(-1,RR+0.5), ylims=(-1,RR+0.5),  
        fa=1, xticks=[], yticks=[],aspect_ratio=:equal)
    end
end
## End "plots points"
current()
end
