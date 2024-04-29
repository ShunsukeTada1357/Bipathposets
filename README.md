# Bipathposets
 A computation for bipath persistent homology using Julia. The computational method is given in the paper "Bipath Persistence" <a href="https://arxiv.org/abs/2404.02536"> arXiv:2404.02536 </a> by Toshitaka Aoki, Emerson G. Escolar, and Shunsuke Tada.
 
# Basic use

We treat bipath filtration of simplicial complexes, which is seen as a pair of filtration sharing same spaces at their ends. 
For example, 
```
FSCa = [[ [[1],1], [[2],1], [[1,2],2] ],  5]
FSCb = [[ [[1],1], [[2],1], [[1,2],3] ],  5]
```
are two filtrations. As for FSCa, simplicies [1] and [2] are born at 1, [1,2] born at 2. No simplicies are born at 3,4, and 5.  

Our main function is "Bipathposets.bipathpersistence" whose arguments are two filtartions of simplicial complexes sharing same spaces at their ends and its output is a dictionaly.
```
bipath = Bipathposets.bipathpersistence(FSCa,FSCb)
```
If we want the persistence of i-th homology group in the bipath filtration, we compute
```
bipath[i]
```
If we want to visualize the persistence of i-th homology group in the bipath persistnce, we compute
```
Bipathposets.plotintlist(aa,0)
```
