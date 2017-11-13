Simple tree shape metrics: imbalance, branchiness
=================================================

The Yule process revisited
--------------------------
- We've previously seen the Yule process in the context of species delimitation, where
  certain algorithms (e.g. GMYC) attempt to find the inflection point between 
  diversification (Yule) and population-genetic (Coalescent) processes.
- The process assumes that every lineage is equally likely to speciate at any given time.
- Hence, the more lineages there are, the shorter the waiting time till the next 
  speciation, because there are more lineages playing the lottery.
- The average expected waiting time will the next speciation is therefore 1/_n_ (or any
  given waiting time, when simulating, is drawn from an exponential distribution).

![](lecture2/exponential-distribution.svg)

Birth/death processes
---------------------
- In the simplest birth/death processes, an additional parameter determines the 
  probability with which every lineage is to go extinct at any given time.
- The _net diversification rate_ is thus `speciation - extinction`.
- (In simulations this rate should obviously be positive, so the extinction rate should
  be smaller than the birth rate in order for there to 
  [grow a tree](http://naturalis.github.io/browbrow).)
  
![](lecture2/birth-death.png)

Diversification through time
----------------------------
But, is net diversification rate (speciation-extinction) constant through time? We might
expect 
[ecological opportunities to arise, triggering adaptive radiations](https://www.nature.com/scitable/knowledge/library/ecological-opportunity-trigger-of-adaptive-radiation-84160951),
processes that we might visualize, qualitatively, as lineage-through-time plots:

![](lecture2/ltt.png)

- (A) Even rates through time, the null hypothesis for patterns of diversification 
  (γ = 0.05). 
- (B) Early burst of cladogenesis and species accumulation, the expected pattern under 
  Ecological Opportunity (γ = -3.39). 
- (C) Late burst of speciation or early extinction (γ = 3.20).

Tree imbalance
--------------
![](lecture2/imbalance.jpg)

Colless' imbalance
------------------
**Colless, DH**, 1982. The theory and practice of phylogenetic systematics. 
_Systematic Zoology_ 31(1): 100-104

![](lecture2/ic.png)

Add up, for all (_n_-1) nodes in a tree with _n_ tips, the absolute difference between
the tips subtended by the child "on the left" and that of the child on the right
(i.e. | _T_<sub>R</sub> - _T_<sub>L</sub> |). Then, normalize this value by dividing
through the maximum value for a tree that size, which is ((_n_-1)*(_n_-2))/2

This value can be computed in R using `ape` and `apTreeShape` thusly:

```R
library(ape)
library(apTreeshape)
tree <- read.tree(text="((A,B),C);")
aptree <- as.treeshape(tree)
ic <- colless(aptree)
```

I2 imbalance
------------
**Mooers AO & Heard SB**, 1997. Inferring evolutionary process from phylogenetic tree 
shape. _Quarterly Review of Biology_ 72: 31–54.

![](lecture2/i2.png)

A perhaps reasonable critique of the _I_<sub>C</sub> index is that it weights "deep"
nodes heavier (consider how the diff between left and right may be much higher for deep
nodes than for shallow ones). An alternative index might therefore, as in this case,
normalize each node right away using _j_ = the number of tips subtended by the focal 
node.

Which one might compute, for example, thusly:

```perl
use Bio::Phylo::IO 'parse_tree';
$ic = parse_tree(
	'-format' => 'newick',
	'-string' => '((A,B),C);',
)->calc_i2;
```

Effect of extinction on imbalance metrics
-----------------------------------------

Empirical results for tree balance
----------------------------------
![](lecture2/phylogenetic-tree-balance-as-a-function-of-tree-size.gif)


Tutorial
--------
https://recology.info/2012/10/phylogenetic-tree-balance/