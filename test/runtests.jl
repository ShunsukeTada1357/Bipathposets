using bipathposets
using AbstractAlgebra
using Combinatorics
using SimpleGraphs 
R =GF(2)
##########################################
##########################################
#Example(1) from red text book P.82.
Test = [[[1],1],[[2],2],[[1,2],3],[[4],4],[[5],5],[[4,5],6],[[2,4],7],[[1,5],8],[[2,5],9],[[2,4,5],10],[[1,2,5],11]]
Test=[Test,11]
aa = Bipathposets.bipathpersistence(Test,Test)
Bipathposets.plotintlist(aa,0)
Bipathposets.plotintlist(aa,1)
##########################################
#Example(3) from our paper's simple exampe.
FSCa = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,4],2],[[1,5],2],[[4,5],2],[[3,4],3],[[3,5],4],[[1,3,5],5]]
FSCa=[FSCa,5]
FSCb = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[3,5],2],[[4,5],3],[[1,4],3],[[1,3,5],3],[[3,4],4]]
FSCb =[FSCb,4]
Bipathposets.baseswithintervals(FSCa)[2]
aa =Bipathposets.bipathpersistence(FSCa,FSCb)
Bipathposets.plotintlist(aa,0)
Bipathposets.plotintlist(aa,1)
##########################################
#Example(3) from our paper's simple exampe.
FSCa = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[6],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[1,4],3],[[4,5],3],[[3,4],3],[[5,6],3],[[3,5],4],[[1,3,5],5]]
FSCa=[FSCa,5]
FSCb = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[6],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[3,5],2],[[4,5],3],[[1,4],3],[[5,6],3],[[1,3,5],3],[[3,4],4]]
FSCb =[FSCb,4]
Bipathposets.baseswithintervals(FSCa)[2]
aa =Bipathposets.bipathpersistence(FSCa,FSCb)
Bipathposets.plotintlist(aa,0)
Bipathposets.plotintlist(aa,1)
using Plots
savefig("bipathPD1th_empty.png") 
#Example (4) 
FSCa = [[ [[1],1], [[2],1], [[1,2],2] ],  5]
FSCb = [[ [[1],1], [[2],1], [[1,2],3] ],  5]
aa =Bipathposets.bipathpersistence(FSCa,FSCb)
Bipathposets.plotintlist(aa,0)
##########################################
#Example(6.1) random model 1
path1, path2 = Bipathposets.get_rectangular_paths([0.4,0.4],[0.7,0.7],5)
n =25
G = Complete(n) 
faces = collect(Combinatorics.powerset(collect(G.V), 2, n))
aa=Bipathposets.clique_random(G,faces,path1,path2)
FSCa = aa[1]
FSCb = aa[2]
aa =Bipathposets.bipathpersistence(FSCa,FSCb)
Bipathposets.plotintlist(aa,1)
