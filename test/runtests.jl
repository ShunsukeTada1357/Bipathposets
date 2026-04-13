using Bipathposets
using Combinatorics
using BenchmarkTools
using AbstractAlgebra
R = GF(2)

function random_GF2_matrix(n::Int, m::Int, p::Float64=0.5)
    return matrix(R, n, m, [R(rand() < p ? 1 : 0) for _ in 1:(n*m)])
end
A = random_GF2_matrix(20, 20, 0.5)

function dense_to_sparse_columns(BB)
    nrows = size(BB, 1)
    ncols = size(BB, 2)

    cols = Vector{Vector{Int}}(undef, ncols)

    for j in 1:ncols
        nz = Int[]
        for i in 1:nrows
            if !iszero(BB[i, j])
                push!(nz, i)
            end
        end
        cols[j] = nz
    end

    return cols
end

##########################################
#Example(1) 
Test = [[[1],1],[[2],2],[[1,2],3],[[4],4],[[5],5],[[4,5],6],[[2,4],7],[[1,5],8],[[2,5],9],[[2,4,5],10],[[1,2,5],11]]
Test=[Test,11]

aa = Bipathposets.baseswithintervals(Test)
aa =Bipathposets.interval_decomposition(Test, Test, verbose=true)

Bipathposets.plot_bipath_diagram(aa,0)
Bipathposets.plot_bipath_diagram(aa,1)

##########################################
#Example(2) 
FSCa = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,4],2],[[1,5],2],[[4,5],2],[[3,4],3],[[3,5],4],[[1,3,5],5]]
FSCa=[FSCa,5]
FSCb = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[3,5],2],[[4,5],3],[[1,4],3],[[1,3,5],3],[[3,4],4]]
FSCb =[FSCb,4]
Bipathposets.baseswithintervals(FSCa)[1]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plot_bipath_diagram(aa,0)
Bipathposets.plot_bipath_diagram(aa,1)
###################################

A1 = [
    0 0 0 0 0 0 0 0 0 0;
    0 1 1 1 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 1 1;
    0 0 0 0 0 0 0 1 0 1;
    1 1 1 0 0 0 0 1 1 1;
    1 0 1 0 0 0 0 0 0 0;
    1 0 1 0 0 0 0 0 0 0;
    1 1 1 0 0 0 0 0 0 0
]

A2 = [
    0 0 0 0 0 0 0 1 1 1;
    0 0 0 0 0 0 0 1 0 1;
    0 0 0 1 0 0 0 1 1 1;
    0 0 0 1 0 0 0 0 0 0;
    0 0 1 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0;
    0 0 1 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0
]

A3 = [
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 1 1 0 0;
    0 0 0 0 0 1 0 1 0 0;
    0 0 0 0 0 1 1 1 0 0;
    0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0
]

A4 = [
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0
]

A5 = [
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 1 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0
]

A_list = [A1, A2, A3, A4, A5]
B_list = Bipathposets.make_or_filtration(A_list)
FSCa = Bipathposets.cubical_filtration_to_FSC(B_list)

aa =Bipathposets.interval_decomposition_cubical(FSCa,FSCa)
Bipathposets.plot_bipath_diagram(aa,0)
Bipathposets.plot_bipath_diagram(aa,1)

cls=Bipathposets.bipath_records_cubical(FSCa, FSCa)
Bipathposets.interactive_bipath_viewer_all(cls, B_list, B_list, FSCa, FSCa; dim=1)
Bipathposets.interactive_bipath_viewer(cls, B_list, B_list, FSCa, FSCa; dim=1)
###################################
#Example(3)
FSCa = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[6],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[1,4],3],[[4,5],3],[[3,4],3],[[5,6],3],[[3,5],4],[[1,3,5],5]]
FSCa=[FSCa,5]
FSCb = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[6],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[3,5],2],[[4,5],3],[[1,4],3],[[5,6],3],[[1,3,5],3],[[3,4],4]]
FSCb =[FSCb,4]
Bipathposets.baseswithintervals(FSCa)[2]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plot_bipath_diagram(aa,0)
Bipathposets.plot_bipath_diagram(aa,1)
using Plots
savefig("bipathPD1th_empty.png") 
##########################################
#Example (4) 
FSCa = [[ [[1],1], [[2],1], [[1,2],2] ],  5]
FSCb = [[ [[1],1], [[2],1], [[1,2],3] ],  5]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plot_bipath_diagram(aa,0)

##########################################
##########################################
#Example (5) 
FSCa =[ [ [[1],1],[[2],1],[[3],1],[[4],1],[[5],1], [[1,2],1], [[1,3],1], [[2,3],1], [[3,4],1], [[3,5],1], [[4,5],1], [[1,5],2], [[1,4],3], [[3,4,5],4],[[1,3,5],5]   ] ,5]
FSCb =[[ [[1],1],[[2],1],[[3],1],[[4],1], [[5],1], [[1,2],1], [[1,3],1], [[2,3],1], [[3,4],1], [[3,5],1], [[4,5],1], [[1,4],2], [[1,5],3], [[1,3,5],3], [[3,4,5],4] ],  4]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plot_bipath_diagram(aa,1)
##########################################

#Example(6)
FSCa =[ [ [[1],1],[[2],1],[[3],1],[[4],1],[[5],1], [[1,2],1], [[1,3],1], [[2,3],1], [[3,4],1], [[3,5],1], [[4,5],1], [[1,5],2], [[1,4],3], [[3,4,5],4],[[1,3,5],5]   ] ,5]
FSCb =[[ [[1],1],[[2],1],[[3],1],[[4],1], [[5],1], [[1,2],1], [[1,3],1], [[2,3],1], [[3,4],1], [[3,5],1], [[4,5],1], [[1,4],2], [[1,5],3], [[1,3,5],3], [[3,4,5],4] ],  4] 
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plot_bipath_diagram(aa,1)
aa = Bipathposets.baseswithintervals(FSCa)

#Example(8)
FSCa = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[1,2],1],[[2,3],1],[[1,3],2],
[[1,4],2],[[1,5],2],[[4,5],2],[[3,4],3],[[3,5],4],[[1,3,5],5]]
FSCa=[FSCa,5]
FSCb = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[1,2],1],[[2,3],1],[[1,3],2],
[[1,5],2],[[3,5],2],[[4,5],3],[[1,4],3],[[3,4],4], [[1,3,5],4]]
FSCb =[FSCb,4]
Bipathposets.baseswithintervals(FSCa)[2]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plot_bipath_diagram(aa,0)
Bipathposets.plot_bipath_diagram(aa,1)

using Plots
savefig("aaa.png")
