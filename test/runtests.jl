using Bipathposets
using Combinatorics
using SimpleGraphs 
##########################################
#Example(1) 
Test = [[[1],1],[[2],2],[[1,2],3],[[4],4],[[5],5],[[4,5],6],[[2,4],7],[[1,5],8],[[2,5],9],[[2,4,5],10],[[1,2,5],11]]
Test=[Test,11]
aa = Bipathposets.interval_decomposition(Test,Test)
Bipathposets.plotintlist(aa,0)
Bipathposets.plotintlist(aa,1)
##########################################
#Example(2) 
FSCa = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,4],2],[[1,5],2],[[4,5],2],[[3,4],3],[[3,5],4],[[1,3,5],5]]
FSCa=[FSCa,5]
FSCb = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[3,5],2],[[4,5],3],[[1,4],3],[[1,3,5],3],[[3,4],4]]
FSCb =[FSCb,4]
Bipathposets.baseswithintervals(FSCa)[2]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plotintlist(aa,0)
Bipathposets.plotintlist(aa,1)
##########################################
#Example(3)
FSCa = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[6],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[1,4],3],[[4,5],3],[[3,4],3],[[5,6],3],[[3,5],4],[[1,3,5],5]]
FSCa=[FSCa,5]
FSCb = [[[1],1],[[2],1],[[3],1],[[4],1],[[5],1],[[6],1],[[1,2],1],[[2,3],1],[[1,3],1],
[[1,5],2],[[3,5],2],[[4,5],3],[[1,4],3],[[5,6],3],[[1,3,5],3],[[3,4],4]]
FSCb =[FSCb,4]
Bipathposets.baseswithintervals(FSCa)[2]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plotintlist(aa,0)
Bipathposets.plotintlist(aa,1)
using Plots
savefig("bipathPD1th_empty.png") 
##########################################
#Example (4) 
FSCa = [[ [[1],1], [[2],1], [[1,2],2] ],  5]
FSCb = [[ [[1],1], [[2],1], [[1,2],3] ],  5]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plotintlist(aa,0)
##########################################
##########################################
#Example (5) 
FSCa =[ [ [[1],1],[[2],1],[[3],1],[[4],1],[[5],1], [[1,2],1], [[1,3],1], [[2,3],1], [[3,4],1], [[3,5],1], [[4,5],1], [[1,5],2], [[1,4],3], [[3,4,5],4],[[1,3,5],5]   ] ,5]
FSCb =[[ [[1],1],[[2],1],[[3],1],[[4],1], [[5],1], [[1,2],1], [[1,3],1], [[2,3],1], [[3,4],1], [[3,5],1], [[4,5],1], [[1,4],2], [[1,5],3], [[1,3,5],3], [[3,4,5],4] ],  4]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plotintlist(aa,0)
##########################################
#Example(6) random model
path1, path2 = Bipathposets.get_rectangular_paths([0.4,0.4],[0.7,0.7],5)
n =25
G = Complete(n) 
faces = collect(Combinatorics.powerset(collect(G.V), 2, n))
aa=Bipathposets.clique_random(G,faces,path1,path2)
FSCa = aa[1]
FSCb = aa[2]
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plotintlist(aa,1)
#Example(7)
FSCa =[ [ [[1],1],[[2],1],[[3],1],[[4],1],[[5],1], [[1,2],1], [[1,3],1], [[2,3],1], [[3,4],1], [[3,5],1], [[4,5],1], [[1,5],2], [[1,4],3], [[3,4,5],4],[[1,3,5],5]   ] ,5]
FSCb =[[ [[1],1],[[2],1],[[3],1],[[4],1], [[5],1], [[1,2],1], [[1,3],1], [[2,3],1], [[3,4],1], [[3,5],1], [[4,5],1], [[1,4],2], [[1,5],3], [[1,3,5],3], [[3,4,5],4] ],  4] 
aa =Bipathposets.interval_decomposition(FSCa,FSCb)
Bipathposets.plotintlist(aa,0)
