LSMM
===

LSMM (Latent Sparse Mixed Model), is an efficient statistical approach to integrating functional annotations with genome-wide association studies. 'LSMM' package provides model parameter estimation as well as statistical inference.

Installation
===========

To install the development version of LSMM, it's easiest to use the 'devtools' package. Note that LSMM depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("mingjingsi/LSMM")
```

Usage
===========

[The 'LSMM' vignette](https://github.com/mingjingsi/LSMM/blob/master/inst/doc/LSMM_package.pdf?raw=true) will provide a good start point for the genetic analysis using LSMM package. The following help page will also provide quick references for LSMM package and the example command lines:

```
library(LSMM)
package?LSMM
```

References
==========

Jingsi Ming, Mingwei Dai, Mingxuan Cai, Xiang Wan, Jin Liu, Can Yang; LSMM: A statistical approach to integrating functional annotations with genome-wide association studies, Bioinformatics, , bty187, https://doi.org/10.1093/bioinformatics/bty187


Reproducibility
==========

All the simulation results can be reproduced by using the code at [sim-LSMM](https://github.com/mingjingsi/sim-LSMM).


Development
==========

This R package is developed by Jingsi Ming and Can Yang (macyang@ust.hk)
