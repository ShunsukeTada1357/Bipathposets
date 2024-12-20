# This file is for the function plotintlist(FSCa,FSCb,i_th::Int) which visualizes intervals in the plane.     
using Plots
using LaTeXStrings
using StatsBase
function intervalstoplane(intL,intR,up,down,center,n,m)
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

#Color of points means the number of intervals in the bipath posets B_{n,m}
function colorfunc(col, max::Int, mult::Int)
    d = div(length(col), max)
    return (d <= 1 || max == 1) ? 1 : min(1 + d * (mult - 1), length(col))
end

function plotpoints(points,n,m)#points=[[a,b],[c,d],[a,b],...]
    plot()
    RR = n + m +2
    ints = StatsBase.countmap(points) # 点の頻度を辞書に 
    maxim= maximum(values(ints))
    col =  cgrad(:turbo)
    num = length(col)
    d = div(num,maxim)
   
    # Axis annotation
    lll=[string(i) for i in 1:n]
    push!(lll,L"\widehat{1}")
    append!(lll, [string(m-i+1)*"'" for i in 1:m])
    uu = -1.5*ones(RR)
    xaxis = vcat(lll,L"\widehat{0}")
    yaxis = vcat(L"\widehat{0}",lll)

    annotate!(1:RR, uu, text.(xaxis,10,"Computer Modern"))
    annotate!(uu,0:RR-1, text.(yaxis,10,"Computer Modern"))

     # Draw axis lines
     vline!([(n+1)%RR], linestyle=:dash, color=:black, label="1-hat ="*string(n+1))
     hline!([(n+1)%RR], linestyle=:dash, color=:black, label="1-hat ="*string(m+1))
     hline!([0], linestyle=:dash, color=:black, label="0-hat ="*string(m+1))
     vline!([RR], linestyle=:dash, color=:black, label="0-hat ="*string(m+n+1)) 

    # Plot points
    for key in keys(ints) 
        (p,q) = key[1] % RR  !=0 ? ((key[1])%RR, (key[2])%RR) : key[2] == RR ? (0, RR) : (RR, key[2])
        plot!([p], [q], 
        proj=:identity, ms=6,msw=:0.3, m=:o,
        c = col[colorfunc(col,maxim,ints[key])],#c = col[colorfunc(col,multiplicity[i]*d)],
        lt=:scatter, markerstrokealpha=1, markeralpha=1, legend = false,
        left_margin = 10Plots.mm, bottom_margin = 10Plots.mm,
        label="",xlims=(-1,RR+0.5), ylims=(-1,RR+0.5),
        fa=1, xticks=[], yticks=[], aspect_ratio=:equal)
    end     
current()
end
#This is the desired function.
function plotintlist(FSCa,FSCb,i_th::Int) 
    X = interval_decomposition(FSCa,FSCb)
    if (i_th in keys(X[1])) == false #check that we have points to plot. 
        println("no ",i_th," homology.")
        return false
    end
    dict_int = X[1][i_th]
    n, m = X[2] - 2, X[3] - 2
    intL, intR, up, center, down = dict_int[1], dict_int[2], dict_int[3], dict_int[4], dict_int[5]
    intlist = intervalstoplane(intL,intR,up,down,center,n,m)
    plotpoints(intlist,n,m)
end