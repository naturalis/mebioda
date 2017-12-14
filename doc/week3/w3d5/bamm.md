Methods for detecting trait-dependent diversification
=====================================================

![](bamm/schluter.jpg)

> [_The Ecology of Adaptive Radiation_](https://global.oup.com/academic/product/the-ecology-of-adaptive-radiation-9780198505235)
> is an important synthesis on the interplay between the selective environment and diversification

Wouldn't it be nice if we could detect _diversification rate shifts_ (so, changes in the net difference
between speciation, `λ`, and extinction, `μ`) in relation to state shifts? This would be our ticket to
detecting _adaptive ratiations_ and _key innovations_.

The Binary-State Speciation and Extinction model (BiSSE)
--------------------------------------------------------

**Maddison WP, Midford PE, Otto SP**, 2007. Estimating a binary character's effect on speciation and extinction.
_Syst Biol_ **56**(5):701-10
doi:[10.1080/10635150701607033](https://doi.org/10.1080/10635150701607033)

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

Follow up on the SSE methods
----------------------------

- BiSSE can be used for hypothesis testing (likelihood ratio tests)
- BiSSE has been expanded to multistate (MuSSE), continuous (QuaSSE), geography (GeoSSE), etc.
- However, several criticisms have been suggested:
  - Maybe not very powerful: need a lot of taxa in a complete, dated tree; trait values should probably
    have some homoplasy, i.e. not one of many synapomorphies for a megadiverse group
    ([Davis et al. (2013)](https://doi.org/10.1186/1471-2148-13-38))
  - False positives: simulating a character onto an existing tree that already has diversification
    rate variation rejects the null ([Rabosky & Goldberg, 2015](https://doi.org/10.1093/sysbio/syu131))
  

Bayesian Analysis of Macroevolutionary Mixtures (BAMM)
------------------------------------------------------
