Species distribution modelling (SDM)
====================================

In this practical you are going to do a species distribution modelling (SDM) analysis. Such an analysis
is done to establish correlations between the occurrences of a species (which you downloaded and cleaned
from GBIF) and the environmental conditions in the locations where it occurs (data for which we will 
download today, if you haven't already). As we have seen, there are different ways in which this can be done,
e.g. using deep learning and various other machine learning techniques, and using 
[maximum entropy](https://en.wikipedia.org/wiki/Principle_of_maximum_entropy) as implemented in the Java
program [maxent](https://biodiversityinformatics.amnh.org/open_source/maxent/), the functionality of which
we will use today.

## Instructions

- [Full Manual](Mebioda_PracticalManual_2019.pdf) - **The sections in here correspond with the code below**
- Powerpoint Nuno - Preparing the data and varaiable selection
- Powerpoint Rosaleen - Discussion of the MaxEnt Results
- [Powerpoint Maarten - Making Maps of SDM (MaxEnt) output in ArcGIS](https://surfdrive.surf.nl/files/index.php/s/rcBYszz1J1GN7PL)
- [Writing the report](reporting.md)
- [Submitting your report](https://github.com/naturalis/mebioda/blob/master/doc/week2/w2d3/lecture3.md#exercise-contributing-to-the-course-repository) -
  deadline: **Wednesday 11 December 2019, 17:00**

## Code

0. [Setting up working environment](00_SettingUpWorkEnviroment.R)
1. [Cropping environmental variables](01_CroppingEnvVariables.R)
2. [Variable selection](02_VariableSelection.R)
3. [Making change maps](03_Making_ChangeMaps.R)

<!--

The SDM practical
-----------------

- [Instructions for the practical](SDM_Workshop_MethodsBiodiversity_08_12_17.pdf)
- [R script for clipping and variable selection](RScript_SDM_Workshop_VariableSelection_Clipping.R)

Lecture slides
--------------

- [SDM and its suite of applications](Principles_SDM_Raes_2017.pptx) - Niels Raes
- [Practical uses of SDM: Forecasting](Presentation_SDM_Forecasting_Leon_Marshall_08_12_17.pptx) - Leon Marshall


- [examples of cropping and conversion using various approaches](https://github.com/naturalis/mebioda/blob/master/doc/week2/w2d3/Workflow.md)
- [SDM Manual](SDM_Manual) - A collection of informational tutorial pages on DIVA-GIS and Maxent
- [MAXENT Data](Maxent_data.zip) - A Zip archive with example input data and Maxent 
  analysis results for the species [_Macaranga auriculata_](http://www.asianplant.net/MacMalBorneo/Macaranga%20auriculata.htm) -->

Background reading
------------------

- a [book chapter](Raes_Aguirre_2018_ch21.pdf) on SDM by Niels Raes and Jesús Aguirre‐Gutiérrez


