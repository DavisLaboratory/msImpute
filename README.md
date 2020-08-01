# msImpute [![DOI](https://zenodo.org/badge/239129382.svg)](https://zenodo.org/badge/latestdoi/239129382)

Methods for label-free mass spectrometry proteomics imputation

**Installation (R)**

```
install.packages("devtools") # devtools is required to download and install the package
devtools::install_github("DavisLaboratory/msImpute")
```

**Quick Start**

```
library(reticulate)
library(msImpute)

selectFeatures(xna)  # xna is a numeric matrix with NAs (for MAR/MNAR diagnosis only)
xna <- scaleData(xna) 
msImpute(xna, rank.max = 2) # rank 2 approximaiton
xcomplete <- msImpute(xna)  # optimal rank determined by msImpute
computeStructuralMetrics(xcomplete, group, xna$E) # "group" denotes experimental condition (e.g. control, treatment etc). Requires python. See Manual for more information.
```

See [user manual](https://github.com/DavisLaboratory/msImpute/blob/master/msImpute_1.2.0.pdf) for help. 


**Reference**

Manuscript in preparation

