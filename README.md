# msImpute
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/239129382.svg)](https://zenodo.org/badge/latestdoi/239129382)
<!-- badges: end -->

Methods for label-free mass spectrometry proteomics imputation

**Installation (R)**

Install from Github:
```
install.packages("devtools") # devtools is required to download and install the package
devtools::install_github("DavisLaboratory/msImpute")
```

Install from Bioconductor:
```
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("msImpute")
```



**Quick Start**

```
library(reticulate)
library(msImpute)

# xna is a numeric matrix with NAs (for MAR/MNAR diagnosis only)
# "group" defines experimental condition (e.g. control, treatment etc).

# select peptides missing in at least one experimental group
selectFeatures(xna, method="ebm", group=group) 


# select peptides that can be informative for
# exploring missing value patterns at high abundance
selectFeatures(xna, method="hvp", n_features=500) 


# impute MAR data by low-rank models (v2 is enhanced version of v1 implementation)
xcomplete <- msImpute(xna, method="v2") 


# impute MNAR data by low-rank models (adaptation of low-rank models for MNAR data)
xcomplete <- msImpute(xna, method="v2-mnar", group=group)  


# Requires python. See Manual for more information.
top.hvp <- findVariableFeatures(xna)
computeStructuralMetrics(xcomplete, 
                         group, 
                         xna[rownames(top.hvp)[1:50],], 
                         k = 2) 


```

Need more help to start? We have collected a number of **case studies** [here](https://github.com/soroorh/proteomicscasestudies/blob/master/msImputeUsersGuide.pdf)


**Reference**

Please consider to cite our [preprint](https://www.biorxiv.org/content/10.1101/2020.08.12.248963v1)

