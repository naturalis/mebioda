Topological analysis
====================

**G-D Yang, P-M Agapow & G Yedi**, 2017. The tree balance signature of mass extinction is 
erased by continued evolution in clades of constrained size with trait-dependent 
speciation. _PLoS ONE_ **12**(6): e0179553
doi:[10.1371/journal.pone.0179553](https://doi.org/10.1371/journal.pone.0179553)

What happens to the signatures of mass extinctions or radiations (such as imbalance or
"stemminess")? Wouldn't they eventually be swamped?

Experimental evolution using [Avida](http://avida.devosoft.org/)
----------------------------------------------------------------

![](lecture/avida.png)

_Avida is a free, open source scientific software platform for conducting and analyzing 
experiments with self-replicating and evolving computer programs. It provides detailed 
control over experimental settings and protocols, a large array of measurement tools, 
and sophisticated methods to analyze and post-process experimental data._

Effect on tree balance of mass extinction
-----------------------------------------

![](lecture/fig1.png)

**Change in tree balance at select time points after mass extinction episode in 
communities of avida digital organisms.**

- Mass extinction treatments were applied _randomly and instantaneously_ (pulse) or by 
  _massive environmental change_ over a period of time (press), at strong and weak 
  intensities.
- The y-axis is Aldous's β [β<sub>A</sub>] a measure of tree balance applicable to 
  non-dichotomous trees; a Yule expectation is around zero, while more negative values 
  indicate trees more imbalanced than this expectation.
- Data points are averages of 100 replicates ± 2 standard errors. Solid traces are 
  maximum likelihood estimates of β<sub>A</sub>, dashed traces are 95% confidence 
  intervals around the calculated β<sub>A</sub> estimates. β<sub>A</sub> values (with 
  confidence intervals) were determined using a customized version of the 
  `maxlik.betasplit` function in the R package apTreeshape (courtesy M. Blum).

1. Simulated phylogenies were generated using the software [MESA](http://datadryad.org/resource/doi:10.5061/dryad.sm379/15).
2. Three mass extinction treatments were employed: 
   - _Random_ = taxa were culled from the tree regardeless of trait value or phylogenetic position. 
   - _Selective-on-diversifiers_ = taxa culled from the tree had the lowes trait value and highest speciation rates.
   - _Selective-on-relicts_ = those taxa with highest trait value and lowest speciation rates were culled preferentally. 
3. Each of the treatments occurred at intensity 90%, 0.75% and 0.5% of all the extant taxa in the tree.
4. The nexus files produced by MESA were then manipulated in R. 
5. The Colless index of imbalance was calculated using the function `colless` for each tree in a time series. 
   [MESA_output files](http://datadryad.org/resource/doi:10.5061/dryad.sm379)
6. For the three treatments combined with the different intensities and the control, the Colless index of imbalance values 
   were collected at the key `CSR time`, `CSR end-mid 1quartile`, `CSR end-mid`, `CSR end-mid 3quartile` and at pre extinction 
   time (`300`), post extinction event (`305`), at the CRS end-mid time (`470`) and at the end of simulation (`600`).
7. the Colless index of imbalance values collected were later rearranged in EXCEL and exported in txt files to be used for the 
   statistical analyses.   
8. The treatments involving a combination of extinction types and intensity were analyzed with two-way ANOVA and Tukey-
   corrected multiple comparison testing, in order to see if there are any significant treatment-by-intensity interactions. 
9. To test whether the various extinction treatment outcomes differed systematically from the control and pre-treatment 
   reference points,  Dunnett's tests were perfomed, using either pre-treatment or the control treatment as the reference 
   standard.
