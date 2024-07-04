# pclogit
[R package] Penalized conditional logistic regression

## Overview

Penalized Conditional (Unconditional) Logistic Regression: Fit a regularization path of conditional (unconditional) logistic regression model for a matched (unmatched) case-control response at a grid of values for regularization parameter lambda. A network-based penalty function was implemented to induce both sparse and smoothing variable selection. It is useful for analysis of high-dimensional data that has a specified correlation structure. Selection probability of each variable is provided.  

## Installation

### Windows users
1. Download and install the latest version of [R](https://cran.r-project.org/bin/windows/base/)
2. Download and install the latest version of [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

### Mac users
1. Download and install the latest version of [R](https://cloud.r-project.org/bin/macosx/)
2. Download and install the Fortran compliler [gfortran](https://cran.r-project.org/bin/macosx/tools/)

### Both Windows and Mac users 
```
## "devtools" package is required if you don't have it.  
install.packages('devtools')

library(devtools)
install_github("statsun78/pclogit")
```

## References

* **Sun, H.** and Wang, S. (2012) Penalized Logistic Regression for High-dimensional DNA Methylation Data with Case-Control Studies, *Bioinformatics* 28(10), p.1368-1375.
* **Sun, H.** and Wang, S. (2013) Network-based Regularization for Matched Case-Control Analysis of High-dimensional DNA Methylation Data, *Statistics in Medicine* 32(12), p.2127-2139.
