# msImpute [![DOI](https://zenodo.org/badge/239129382.svg)](https://zenodo.org/badge/latestdoi/239129382)
<!-- badges: start -->
[![Travis build status](https://travis-ci.org/DavisLaboratory/msImpute.svg?branch=master)](https://travis-ci.org/DavisLaboratory/msImpute)
<!-- badges: end -->

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


# Requires python. See Manual for more information.
top.hvp <- findVariableFeatures(xna$E)
computeStructuralMetrics(xcomplete, 
                         # "group" denotes experimental condition (e.g. control, treatment etc).
                         group, 
                         xna$E[rownames(top.hvp)[1:50],], 
                         k = 2) 


```

Need more help to start? We have collected a number of case studies [here](https://github.com/soroorh/proteomicscasestudies/blob/master/msImpute_user_guide.pdf)


**Reference**

Cite our [preprint](https://www.biorxiv.org/content/10.1101/2020.08.12.248963v1)

