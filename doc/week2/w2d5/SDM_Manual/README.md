Species Distribution Modelling Manual
=====================================

Adapted from a previous version by Dr. Niels Raes
-------------------------------------------------

### Outline

This manual gives a run-through of SDM on Windows computers. As such, some of the steps cannot
be performed on other operating systems with the programs described. However, other programs
may be equally suitable. Note also that the exercises in this manual are not a requirement of
the course, they are only provided here for informational purposes.

1. **[GIS and spatial modelling](1_Pointdata)** - 
   Chapter 1 starts with data pre-processing and some of the configuration challenges in dealing 
   with lat/lon data on Windows computers. Subsequently, the construction of base maps with 
   administrative features (country borders) is demonstrated using the DIVA-GIS program.
2. **[Preparing GIS data layers for Species Distribution Modelling](2_Data_layers)** -
   Chapter 2 continues by introducing the usage of DIVA-GIS in converting and cropping GIS data
   layers, such as those containing (bio)climatic data.
3. **[Species Distribution Modelling with MaxEnt](3_Modelling)** - 
   Chapter 3 introduces SDM using the maximum entropy algorithm as implemented in the portable
   Java-based program Maxent. Considerable attention is given to validation of the results, which
   are then visualized and transformed in DIVA-GIS. The chapter concludes with a (non-compulsory)
   exercise.
4. **[Projecting species distribution models to future and past climate scenarios, and to different geographical regions](4_Projecting)** - 
   Chapter 4 demonstrates how SDM can be used to project species' responses to a changing environment,
   i.e. temperature changes.