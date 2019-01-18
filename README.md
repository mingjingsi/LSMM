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

Jingsi Ming, Mingwei Dai, Mingxuan Cai, Xiang Wan, Jin Liu, Can Yang; LSMM: A statistical approach to integrating functional annotations with genome-wide association studies, Bioinformatics, Volume 34, Issue 16, 15 August 2018, Pages 2788–2796, https://doi.org/10.1093/bioinformatics/bty187


Reproducibility
==========

All the simulation results can be reproduced by using the code at [sim-LSMM](https://github.com/mingjingsi/sim-LSMM). We also provide an example to reproduce the results of real data analysis at [sim-LSMM/realdata_example.R](https://github.com/mingjingsi/sim-LSMM/blob/master/realdata_example.R). The data sets analyzed in this paper, including nine genic category annotations, 127 cell-type specific functional annotations, and the summary statistics of 30 GWASs can be download [here](https://drive.google.com/drive/folders/1YrhncaNpo6doHfnIYqXPlhPdE7YpqEjO).


Development
==========

This R package is developed by Jingsi Ming and Can Yang (macyang@ust.hk)
