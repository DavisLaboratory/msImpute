
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/239129382.svg)](https://zenodo.org/badge/latestdoi/239129382)
<!-- badges: end -->

msImpute - Methods for label-free mass spectrometry proteomics imputation
========================================

MsImpute is a R package for imputation of peptide intensity in proteomics experiments.
It additionally contains tools for MAR/MNAR diagnosis and assessment of distortions to the probability 
distribution of the data post imputation.  

The missing values are imputed by low-rank approximation of the underlying data matrix if they are MAR (method = "v2"), by Barycenter approach if missingness is MNAR ("v2-mnar"), or by Peptide Identity Propagation (PIP). While "v2" approach is more appropriate for imputation of data acquired by DIA, "v2-mnar" is designed for imputation of DDA, TMT and time-series datasets. However, the true dynamic range can not be reliably recovered by imputation, particularly in datasets with small sample sizes (for example, 3-5 replicates per experimental condition). Our PIP approach infers the missing intensity values for an identification based on similarity of LC-MS features of peptide-like signals detected in MS1 (e.g. by a feature detector) and the identified peptides. We currently support MaxQuant outputs, including DDA-PASEF datasets. **We strongly recommend the PIP approach for imputation of time-series, or datasets that suffer from large (> 50%) missing values per run**. Our PIP enhances data completeness, while reporting *weights* that measure the confidence in propagation. These can be used as observation-level weights in *limma* linear models to improve differential expression testing, by incorporating the uncertainty in intensity values that are inferred by PIP into the model.


Installation
--------------

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



Quick Start
----------------

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

New feature : msPIP
---------------------

**PIP reduces the proportion of missing values substantially, almost to zero.**

Figure: The proportion of missing peptides per sample in PASS00589 DDA dataset before and after PIP.

<img src="https://user-images.githubusercontent.com/7257233/121839424-5c3a4380-cd1d-11eb-84fa-437a387c44f2.png" width="700px" align="center">


**PIP recovers the low abundance peptides and re-constructs the true dynamic range**

Low-abundance peptides not quantified by MaxQuant are recovered, and differential expression results are improved. Note down regulated peptides that are not present in the volcano plot of DE test on MQ-reported data (bottom left), that are recovered by PIP (bottom right volcano plot) for the same experimental contrast.

<img src="https://user-images.githubusercontent.com/7257233/121839859-55600080-cd1e-11eb-998e-f7e60896b1bf.png" width="700px" align="center">


Need more help to start? We have collected a number of **case studies** [here](https://github.com/soroorh/proteomicscasestudies/blob/master/msImputeUsersGuide.pdf)


**Reference**

Please consider to cite our [preprint](https://www.biorxiv.org/content/10.1101/2020.08.12.248963v1)

