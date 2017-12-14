Teach the controversy!
======================

![](bamm/schluter.jpg)

> [_The Ecology of Adaptive Radiation_](https://global.oup.com/academic/product/the-ecology-of-adaptive-radiation-9780198505235)
> is an important synthesis on the interplay between the selective environment and diversification

Wouldn't it be nice if we could detect _diversification rate shifts_ (so, changes in the net difference
between speciation, `λ`, and extinction, `μ`) in relation to state shifts? This would be our ticket to
detecting _adaptive ratiations_ and _key innovations_.

What are some of the options?
-----------------------------

**Maddison WP, Midford PE, Otto SP**, 2007. Estimating a binary character's effect on speciation and extinction.
_Syst Biol_ **56**(5):701-10
doi:[10.1080/10635150701607033](https://doi.org/10.1080/10635150701607033)

- The **Bi**nary-**S**tate **S**peciation and **E**xtinction model (**BiSSE**)
- Assumes that: 
  - An accurate rooted phylogenetic tree with branch lengths is known (the “inferred tree”) 
  - The character state is known for each of the terminal taxa
  - The tree is assumed complete: all extant species in the group have been found and included 
  - The tree is ultrametric (i.e., the total root-to-tip distance is the same for all tips)
- Estimates six parameters: 
  - The instantaneous rates of speciation λ<sub>0</sub> and extinction μ<sub>0</sub>, 
    when the lineage is in state `0` (e.g., _herbivory_)
  - These rates (λ<sub>1</sub>, μ<sub>1</sub>) when the lineage is in state `1` (e.g., _carnivory_)
  - The instantaneous rates of character state change (`0 ⟶ 1`, _q_<sub>01</sub> and `1 ⟶ 0`, _q_<sub>10</sub>)

How BiSSE estimates its model parameters
----------------------------------------

![](bamm/bisse.gif)

- D<sub>N0</sub>(t) or D<sub>N1</sub>(t) are the probabilities that a lineage beginning at time _t_ with state 0 or 1 
  evolves into a clade as observed descending from node _N_
- Differential equations track the changes in the parameters as the method traverses from the tips to the root to 
  calculate the likelihood
- (This exploits the general pattern known as _Felsenstein's pruning algorithm_)

The BAMM situation
------------------
