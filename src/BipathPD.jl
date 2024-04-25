##This file is for a function plotintlist(intlist,a,b)
##That visualize bipath intervals. Input: intervals (represented by dimension vector), a,b: integers for rotating PD.    

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


function intmultiplicity(I, intlist)
    return length([J for J in intlist if J ==I])
end

function deletespecific(I,intlist)
    return [J for J in intlist if J!= I]
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

#intervalstocartesiancoordinate(intL,intR,others,4,4)
function intervalstocartesiancoordinate(intL,intR,up,down,n,m)
    intlist =[]
    up 
    down
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
    return intlist#, center
end

function intlisttointlistwithmultiplicity(intlist)#[I,I,J,J,J]=> [[I,2],[J,3]]
    a=copy(intlist)
    want=[]
    while iszero(a) != true 
        push!(want,[a[1],intmultiplicity(a[1], a) ])
        a=deletespecific(a[1],a)
    end
    return want
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

####################################################
##########################################################################################
#the first input is the result of a function "bipathpersistence" in BipathmatrixMthod.jl
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
    down= dict_int[5]
    intlist = intervalstocartesiancoordinate(intL,intR,up,down,n,m)
    plot()
    x_ax, y_ax, multiplicity = intervalstocartesiancoordinatewithmulti(intlist)
    maxim = maximum(multiplicity)
    #col = cgrad(:lajolla50,rev =true )
    col =  cgrad(:turbo )

    num = length(col)
    d = div(num,maxim)

## draw labels on a plane
a=0
b=0
vline!([(n+1+a)%RR], linestyle=:dash, color=:black, label="1-hat ="*string(n+1))
hline!([(n+1+b)%RR], linestyle=:dash, color=:black, label="1-hat ="*string(m+1))

vline!([a%RR], linestyle=:dash, color=:black, label="0-hat ="*string(n+1))

hline!([b%RR], linestyle=:dash, color=:black, label="0-hat ="*string(m+1))
vline!([RR], linestyle=:dash, color=:black, label="0-hat ="*string(m+n+1))

#hline!([RR+1], linestyle=:dash, color=:black, label="∞ ="*string(RR+1))

#hat = L"\widehat{0}"
hat = " "
lll=[hat]
for i in 1:n 
    push!(lll,string(i)) 
end
hat = L"\widehat{1}"

push!(lll,hat)
for i in 1:m 
    push!(lll,string(m-i+1)*"'") 
end
uu = -1.5*ones(RR)
annotate!(0:RR-1, uu, text.(lll,16,"Computer Modern"))
uu = -1.5*ones(RR)
annotate!(uu,0:RR-1, text.(lll,16,"Computer Modern"))
#annotate!([-1.5],RR:RR, text.(uuuu,16,"Computer Modern"))
#plot!(annotation = (RR, -1.5,L"\widehat{0}"))

## End "draw labels on a plane (or torus)"

##plot points
for i in 1:length(x_ax)
    plot!([(x_ax[i]+a)%RR], [(y_ax[i]+b)%RR], 
    proj=:identity, ms=6,msw=:0.3,
    m=:o,
    c = col[colorfunc2(col,maxim,multiplicity[i])],#c = col[colorfunc(col,multiplicity[i]*d)],
    lt=:scatter,
    markerstrokealpha=1, legend = false,markeralpha=1,
    xlims=(-1,RR+0.5), ylims=(-1,RR+0.5),  
    fa=1, xticks=[], yticks=[],
    grid=true, left_margin = 10Plots.mm, 
    bottom_margin = 10Plots.mm, 
    aspect_ratio=:equal, label="")
end
## End "plots points "
current()
end

#savefig("bipathPD3.png") 
#plotintlist(intlist,a,b)
