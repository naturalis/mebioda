Topological analysis
--------------------
1.Simulated phylogenies were generated using the software [MESA](http://datadryad.org/resource/doi:10.5061/dryad.sm379/15).

2.Three mass extintion treatments were employed: 
Random = taxa were culled from the tree regardeless of trait value or phylogenetic position.
Selective-on-diversifiers = taxa culled from the tree had the lowes trait value and highest speciation rates.
Selective-on-relicts = those taxa with highest trait value and lowest speciation rates were culled preferentally. 

   2. 2a Each of the treatments occurred at intensity 90%, 0.75% and 0.5% of all the extant taxa in the tree.
  
3. The nexus files produced by MESA were then manipulated in R. 

   3. 3a The Colless index of imbalance was calculated using the function *colless* for each tree in a time series. [MESA_output files](http://datadryad.org/resource/doi:10.5061/dryad.sm379)
  
4.  For the three treatments combined with the different intensities and the control, the Colless index of imbalance values were collected at the key CSR time, CSR end-mid 1quartile, CSR end-mid, CSR end-mid 3quartile and at pre extintion time (300), post extintion event (305), at the CRS end-mid time (470) and at the end of simulation (600).
   
    4. 4a the Colless index of imbalance values collected were later rearranged in EXCEL and exported in txt files to be used for the statistical analyses.
   
5. The treatments involving a combination of extintion types and intensity were analyzed with two-way ANOVA and Tukey-corrected multiple comparison testing, in order to see if there are any significant treatment-by-intensity interactions. 
6. To test whether the various extinction treatment outcomes differed systematically from the control and pre-treatment reference points,  Dunnett's tests were perfomed, using either pre-treatment or the control treatment as the reference standard.
