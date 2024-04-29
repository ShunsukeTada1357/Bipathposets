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




