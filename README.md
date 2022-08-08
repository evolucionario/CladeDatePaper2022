### Code used in Claramunt 2022 CladeDate: empirical calibration density generator for divergence time estimation. Methods in Ecology and Evolution.

Santiago Claramunt

Department of Natural History, Royal Ontario Museum, and
Department of Ecology and Evolutionary Biology, University of Toronto, Ontario, Canada.

E-mail: claramunt.bio@gmail.com

##Content:

# CladeDate package

CladeDate is an R package for the generation of empirical calibration information from the fossil record. CladeDate uses simple mathematical models to estimate the age of a clade and its uncertainty based on fossil ages. Using a Monte Carlo approach, CladeDate generates empirical densities representing the uncertainty associated with the age of the clade and fits standard probability density functions that can be used in time-tree inference software such as BEAST2, MrBayes, and MCMCtree.

Instalation from the R console:

````
library(devtools)
install_github("evolucionario/CladeDate")
library(CladeDate)
````

# Code used for simulations

  R code used for generating phylogenetic trees, DNA sequencies, and fossil records
  
  R code for analyzing the simulated data using CladeDate and the time-tree estimation function 'chronos' (ape package) 

  R code for the comparison between CladeDate + chronos and CladeAge + BEST
