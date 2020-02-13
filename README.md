# msImpute [![DOI](https://zenodo.org/badge/239129382.svg)](https://zenodo.org/badge/latestdoi/239129382)

Methods for mass spectrometry label-free proteomics imputation

**Install (R)**

```
devtools::install_github("DavisLaboratory/msImpute")

```

**Run**

```
library(msImpute)

msImpute(xna) # xna is a numeric matrix with NAs
msImpute(xna, rank.max = 2) # rank 2 approximaiton
```

See `?msImpute` for help. 


**Reference**

Manuscript under preparation

