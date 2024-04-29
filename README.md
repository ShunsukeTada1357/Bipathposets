# Bipathposets
 A computation for bipath persistent homology using Julia. The computational method is given in the paper "Bipath Persistence" <a href="https://arxiv.org/abs/2404.02536"> arXiv:2404.02536 </a> by Toshitaka Aoki, Emerson G. Escolar, and Shunsuke Tada.
 
# Basic use

We treat bipath filtration of simplicial complexes, which is seen as a pair of filtration sharing same spaces at their ends. 
For example, 
```
julia> FSCa = [[ [[1],1], [[2],1], [[1,2],2] ],  5]
julia> FSCb = [[ [[1],1], [[2],1], [[1,2],3] ],  4]
```
are two filtrations. As for FSCa, simplicies [1] and [2] are born at 1, [1,2] born at 2. No simplicies are born at 3,4, and 5.  

Our main function is "Bipathposets.bipathpersistence" whose arguments are two filtartions of simplicial complexes sharing same spaces at their ends and its output is a list with three elements. the first is a dictionaly and the second, third are integers meaning the length of each filtration.
```
julia> bipath = Bipathposets.bipathpersistence(FSCa,FSCb)
(Dict{Any, Any}(0 => Vector{Any}[[[[1, 1], [1, 2]]], [], [], [[1, 5]], []]), 5, 4)
```
In calculatin the follwing will be printed.
```
 ∃ 0_th homology, #[̂0,̂1] is 1
intervals with ̂0: <1', ̂0>
intervals with ̂1:
intervals up:
intervals down:
```
The explanation of  the notation <-,->, ̂1, and ̂0 can be seen the above paper. 

If we want the persistence of i-th homology group in the bipath filtration, we compute
```
bipath[1][i]
```
If we want to visualize the persistence of i-th homology group in the bipath persistnce, we compute
```
Bipathposets.plotintlist(bipath,i)
```
For exampe, let i be 0, we obtain the following diagram.

  <img src="bipath.jpg" alt="bipath persistence diagram" width="200px" align="center">

