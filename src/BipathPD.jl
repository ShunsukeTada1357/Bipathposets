##This file is for a function plotintlist(intlist,a,b)
## that visualize bipath intervals. Input: intervals (represented by dimension vector), a,b: integers for rotating PD.    

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


function shaffle(list,a)##change name!! want good name
    r = length(list)
    a = a%(r)
    if a == 0
        return list
    else
        want=list[a+1:r]
        for i in 1:a
            push!(want,list[i])
        end
        return want
    end
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

    #center = others[2]
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

function intervalstocartesiancoordinatewithmulti(intlist,n,m)
    newintlistwith = intlisttointlistwithmultiplicity(intlist)
    #n = length(Intlist[1][1])-2
    #m = length(Intlist[1][2])-2
    #R =n+m+2
    xlist = []
    ylist = []
    multi = []
    for I in newintlistwith
        push!(xlist, I[1][1])
        push!(ylist, I[1][2])
        push!(multi, I[2])
    end
    return xlist,ylist,multi,n,m
end



#(5) color of points means the number of intervals C_{n,m}
function plotintlist(intL,intR,up,down,a,b,n,m)#####a=0,b=0
    intlist = intervalstocartesiancoordinate(intL,intR,up,down,n,m)
    plot()
    XX = intervalstocartesiancoordinatewithmulti(intlist,n,m)
    #a,b = 0,0
    maxim = maximum(XX[3])
    col =  cgrad(:lajolla50,rev =true )
    num = length(col)
    d = div(num,maxim)
    n = XX[4]
    m = XX[5]
    R = n + m +2

#color area where points are not plotted
    #plot!([0, n+1], [0, n+1], fillrange=0, fillalpha=0.4, fillcolor=:lightblue,linewidth=0, msw=:0,label=false)
    #plot!([5,R],[n+1,R],fillrange = [n+1,n+1],fillalpha = 0.4,fillcolor =:lightblue, msw=:0,linewidth=0, label=false)
    
for i in 1:length(XX[1])
    #sizeofpoint= 3 + I[2]/all
    plot!([(XX[1][i]+a)%R], [(XX[2][i]+b)%R], 
    proj=:identity, ms=6,msw=:0,
    m=:o, c =col[XX[3][i]*d], lt=:scatter,
    markerstrokealpha=1, legend = false,
    xlims=(-1,R), ylims=(-1,R),  
    fa=0.4, xticks=[], yticks=[],
    grid=true, left_margin = 10Plots.mm, 
    bottom_margin = 10Plots.mm, 
    aspect_ratio=:equal, label="")
end

vline!([(n+1+a)%R], linestyle=:dash, color=:red, label="1-hat ="*string(XX[4]+1))
hline!([(n+1+b)%R], linestyle=:dash, color=:red, label="1-hat ="*string(XX[5]+1))
vline!([a%R], linestyle=:dash, color=:blue, label="0-hat ="*string(XX[4]+1))
hline!([b%R], linestyle=:dash, color=:blue, label="0-hat ="*string(XX[5]+1))
# draw labels

hat = L"\widehat{0}"
lll=[hat]
for i in 1:n 
    push!(lll,string(i)) 
end
hat = L"\widehat{1}"

push!(lll,hat)
for i in 1:m 
    push!(lll,string(m-i+1)*"'") 
end
uu = -1.5*ones(R)
uuuu=shaffle(lll,(R-a)%R)
annotate!(0:R-1, uu, text.(uuuu,18,"Computer Modern"))
uu = -1.5*ones(R)
uuuu=shaffle(lll,(R-b)%R)
annotate!(uu,0:R-1, text.(uuuu,18,"Computer Modern"))
current()
end

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

function upper_than_diag(x)
    if x[1]<= x[2]
      return true
    else
      return false
    end
  end

function plotintlist(FSCa,FSCb,a,b,i_th)#####a=0,b=0
    X =bipathpersistence(FSCa,FSCb,i_th)[i_th+1]
    n = FSCa[2]-2
    m = FSCb[2]-2
    intL = X[1]
    intR=X[2]
    up = X[3]
    down=X[5]
    intlist = intervalstocartesiancoordinate(intL,intR,up,down,n,m)
    plot()
    XX = intervalstocartesiancoordinatewithmulti(intlist,n,m)
    #a,b = 0,0
    maxim = maximum(XX[3])
    col =  cgrad(:lajolla50,rev =true )
    num = length(col)
    d = div(num,maxim)
    n = XX[4]
    m = XX[5]
    RR = n + m +2

#color area where points are not plotted
    #plot!([0, n+1], [0, n+1], fillrange=0, fillalpha=0.4, fillcolor=:lightblue,linewidth=0, msw=:0,label=false)
    #plot!([5,R],[n+1,R],fillrange = [n+1,n+1],fillalpha = 0.4,fillcolor =:lightblue, msw=:0,linewidth=0, label=false)
    
for i in 1:length(XX[1])
    #sizeofpoint= 3 + I[2]/all
    plot!([(XX[1][i]+a)%RR], [(XX[2][i]+b)%RR], 
    proj=:identity, ms=6,msw=:0,
    m=:o#,c =col[min(length(col),XX[3][i]*d)]
    ,
    c = col[colorfunc(col,XX[3][i]*d)],
    lt=:scatter,
    markerstrokealpha=1, legend = false,
    xlims=(-1,RR), ylims=(-1,RR),  
    fa=0.4, xticks=[], yticks=[],
    grid=true, left_margin = 10Plots.mm, 
    bottom_margin = 10Plots.mm, 
    aspect_ratio=:equal, label="")
end

vline!([(n+1+a)%RR], linestyle=:dash, color=:red, label="1-hat ="*string(XX[4]+1))
hline!([(n+1+b)%RR], linestyle=:dash, color=:red, label="1-hat ="*string(XX[5]+1))
vline!([a%RR], linestyle=:dash, color=:blue, label="0-hat ="*string(XX[4]+1))
hline!([b%RR], linestyle=:dash, color=:blue, label="0-hat ="*string(XX[5]+1))
# draw labels

hat = L"\widehat{0}"
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
uuuu=shaffle(lll,(RR-a)%RR)
annotate!(0:RR-1, uu, text.(uuuu,18,"Computer Modern"))
uu = -1.5*ones(RR)
uuuu=shaffle(lll,(RR-b)%RR)
annotate!(uu,0:RR-1, text.(uuuu,18,"Computer Modern"))
current()
end

#the first input is the result of a function "bipathpersistence" in BipathmatrixMthod.jl
function plotintlist(dict_resultof_bipathpersistence,a::Int,b::Int,i_th::Int)#####a=0,b=0
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
    n = X[2]-2
    m = X[3]-2
    RR = n + m +2
    
    intL = dict_int[1]
    #intL =  [X[1][1][i] for i in findall(upper_than_diag,X[1][1])]
    intR= dict_int[2]
    #intR =  [X[2][1][i] for i in findall(upper_than_diag,X[2][1])]
    infinity = dict_int[4]
    #println("#[̂0,̂1] is ", length(infinity))
    
    #up =  [X[3][i] for i in findall(upper_than_diag,X[3])]
    #down=  [X[5][i] for i in findall(upper_than_diag,X[5])]
    up =  dict_int[3]
    down= dict_int[5]
    intlist = intervalstocartesiancoordinate(intL,intR,up,down,n,m)
    plot()
    XX = intervalstocartesiancoordinatewithmulti(intlist,n,m)
    #a,b = 0,0
    maxim = maximum(XX[3])
    col =  cgrad(:lajolla50,rev =true )
    num = length(col)
    d = div(num,maxim)


#color area where points are not plotted
    #plot!([0, n+1], [0, n+1], fillrange=0, fillalpha=0.4, fillcolor=:lightblue,linewidth=0, msw=:0,label=false)
    #plot!([5,R],[n+1,R],fillrange = [n+1,n+1],fillalpha = 0.4,fillcolor =:lightblue, msw=:0,linewidth=0, label=false)
    
for i in 1:length(XX[1])
    #sizeofpoint= 3 + I[2]/all
    plot!([(XX[1][i]+a)%RR], [(XX[2][i]+b)%RR], 
    proj=:identity, ms=6,msw=:0,
    m=:o#,c =col[min(length(col),XX[3][i]*d)]
    ,
    c = col[colorfunc(col,XX[3][i]*d)],
    lt=:scatter,
    markerstrokealpha=1, legend = false,
    xlims=(-1,RR), ylims=(-1,RR),  
    fa=0.4, xticks=[], yticks=[],
    grid=true, left_margin = 10Plots.mm, 
    bottom_margin = 10Plots.mm, 
    aspect_ratio=:equal, label="")
end

vline!([(n+1+a)%RR], linestyle=:dash, color=:red, label="1-hat ="*string(XX[4]+1))
hline!([(n+1+b)%RR], linestyle=:dash, color=:red, label="1-hat ="*string(XX[5]+1))
vline!([a%RR], linestyle=:dash, color=:blue, label="0-hat ="*string(XX[4]+1))
hline!([b%RR], linestyle=:dash, color=:blue, label="0-hat ="*string(XX[5]+1))
# draw labels

hat = L"\widehat{0}"
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
uuuu=shaffle(lll,(RR-a)%RR)
annotate!(0:RR-1, uu, text.(uuuu,18,"Computer Modern"))
uu = -1.5*ones(RR)
uuuu=shaffle(lll,(RR-b)%RR)
annotate!(uu,0:RR-1, text.(uuuu,18,"Computer Modern"))
current()
end

dd = Dict()
dd[1]= "a"
1 in keys(dd)
###(6)Exampleintlist

#savefig("bipathPD3.png") 
#plotintlist(intlist,a,b)

