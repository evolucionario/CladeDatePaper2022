# Code used in "CladeDate: calibration information generator for divergence time estimation" (Claramunt 2022, _Methods in Ecology and Evolution_)

**Santiago Claramunt**

_Department of Natural History, Royal Ontario Museum, and Department of Ecology and Evolutionary Biology, University of Toronto_

Functions and code used in the paper "CladeDate: empirical calibration density generator for divergence time estimation" (Claramunt 2022 Methods in Ecology and Evolution), including the first published version of the R package CladeDate (Version 1.0). The most updated version of CladeDate is in https://github.com/evolucionario/cladedate.


## Content:

### CladeDate package version 1.0

`CladeDate` is an `R` package for the generation of empirical calibration information from the fossil record. `CladeDate` uses simple mathematical models to estimate the age of a clade and its uncertainty based on fossil ages. Using a Monte Carlo approach, `CladeDate` generates empirical densities representing the uncertainty associated with the age of the clade and fits standard probability density functions that can be used in time-tree inference software such as `BEAST2`, `MrBayes`, and `MCMCtree`.

Instalation from the `R` console:

````
library(devtools)
install_github("evolucionario/CladeDatePaper2022/CladeDatePackage")
library(CladeDate)
````

### Simulations used for the validation of the CladeDate algorithm

  1. Simulation of Datasets: `R` code used for generating phylogenetic trees, DNA sequencies, and fossil records.
  
  2. Simulations CladeDate with Chronos: `R` code for analyzing the simulated data using `CladeDate` and the time-tree estimation function `chronos` in the `ape` package (Paradis & Schliep 2019).

  3. Simulation for CladeAge Analysis: `R` code for the comparison between `CladeDate` + `chronos` and `CladeAge` + `BEAST` (Matschiner _et al._ 2017).


### References

- Claramunt, S. 2022 CladeDate: calibration information generator for divergence time estimation. _Methods in Ecology and Evolution_.

- Matschiner, M., Musilová, Z., Barth, J. M. I., Starostová, Z., Salzburger, W., Steel, M., & Bouckaert, R. R. (2017). Bayesian phylogenetic estimation of clade ages supports trans-Atlantic dispersal of cichlid fishes. _Systematic Biology_ 66(1):3–22.

- Paradis, E. & Schliep, K. (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. _Bioinformatics_ 35:526-528.
