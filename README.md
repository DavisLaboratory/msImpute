
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/239129382.svg)](https://zenodo.org/badge/latestdoi/239129382)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-msimpute/badges/downloads.svg)](https://anaconda.org/bioconda/bioconductor-msimpute)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-msimpute/badges/license.svg)](https://anaconda.org/bioconda/bioconductor-msimpute)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-msimpute/badges/version.svg)](https://anaconda.org/bioconda/bioconductor-msimpute)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/bioconductor-msimpute/badges/latest_release_date.svg)](https://anaconda.org/bioconda/bioconductor-msimpute)
<!-- badges: end -->


msImpute - Methods for label-free mass spectrometry proteomics imputation
========================================

MsImpute is a R package for imputation of peptide intensity in proteomics experiments.
It additionally contains tools for MAR/MNAR diagnosis and assessment of distortions to the probability 
distribution of the data post imputation.  

The missing values are imputed by low-rank approximation of the underlying data matrix if they are MAR (method = "v2"), by Barycenter approach if missingness is MNAR ("v2-mnar"), or by Peptide Identity Propagation (PIP). While "v2" approach is more appropriate for imputation of data acquired by DIA, "v2-mnar" is designed for imputation of DDA, TMT and time-series datasets. However, the true dynamic range can not be reliably recovered by imputation, particularly in datasets with small sample sizes (for example, 3-5 replicates per experimental condition). 

Our PIP approach infers the missing intensity values for an identification based on similarity of LC-MS features of peptide-like signals detected in MS1 (e.g. by a feature detector) and the identified peptides. We currently support MaxQuant outputs, including DDA-PASEF datasets. **We strongly recommend the PIP approach for imputation of time-series, or datasets which suffer from large (> 50%) missing values per run**. Our PIP enhances data completeness, while reporting *weights* that measure the confidence in propagation. These can be used as observation-level weights in *limma* linear models to improve differential abundance testing, by incorporating the uncertainty in intensity values that are inferred by PIP into the model. **We have given a demo of PIP approach on a published DDA dataset below.**


Installation
--------------
**Please note R version 4.1.1 or later is required**

Install from Github:

```{r}
install.packages("devtools") # devtools is required to download and install the package
devtools::install_github("DavisLaboratory/msImpute")
```

Install from Bioconductor:
```{r}
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("msImpute")
```



Quick Start
----------------

```{r}

library(msImpute)

# Let xna be a numeric matrix of (unormalised) log-intensity with NAs 
# Let "group" defines a single experimental condition (e.g. control, treatment etc).
# Let "design" defines the experimental design (e.g. model.matrix(~0+group+batch)).

# select peptides missing in at least one experimental group
group <- factor(c('control','control','conditionA','conditionA'))
selectFeatures(xna, method="ebm", group=group) 


# select peptides that can be informative for
# exploring missing value patterns at high abundance
selectFeatures(xna, method="hvp", n_features=500) 


# Impute MAR data by low-rank approximation (v2 is enhanced version of v1 implementation tailored to small data)
xcomplete <- msImpute(xna, method="v2") 


# Impute complex MV mechanims (MNAR and MAR) as mixture of two normal distributions (known as the Barycenter approach) 
design <- model.matrix(~0+group+batch)
xcomplete <- msImpute(xna, method="v2-mnar", design=design)  


# Allow for features with very few (less than 4) measurements
xcomplete <- msImpute(xna, method="v2-mnar", design=design, relax_min_obs = TRUE)

# Rank-2 approximation for the modeling MAR MVs in small sample regimes
xcomplete <- msImpute(xna, method="v2-mnar", design=design, relax_min_obs = TRUE, rank.max = 2)


# Disable seed generator such that the lower component of the mixture corresponding to MNAR is stochastic and returns a different results with each call (Note this is not recommended for reproducibility)
xcomplete <- msImpute(xna, method="v2-mnar", design=design, relax_min_obs = TRUE, rank.max = 2, use_seed = FALSE)

```

News
---------------------
**22.03.2025**

The following changes have been made to function calls:
- The use of 'group' is now deprecated. msImpute now allows specifying a design matrix (which has to have zero intercept) to accommodate more complex missing value (MV) data generation processes such as LC batch.
- The new version models log-intensity as a mixture of two normal distributions, one for the MAR and one for the MNAR component. The weights of the mixture (equivalent to `a` or `alpha` in the old API) are determined according to a Dirichlet distribution learned from mv patterns, so you no longer need to specify the weights of the two distributions manually.
- The new version also allows for retaining peptides/proteins with very few measurements (e.g. less than 4) via `relax_min_obs`.
- In the old API, imputation was set to be deterministic for reproducibility purposes. If you wish to keep it stochastic for the lower component of the mixture that corresponds to MNAR distribution (sampling from down-shifted distribution) please set the use_seed argument.

The following dependencies were removed:
- reticulate
- scran

The following functions are deprecated:
- computeStructuralMetrics()

Tutorials 
---------------------
Example workflows can be found under `figures/` in the [reproducibility repository](https://github.com/DavisLaboratory/msImpute-reproducibility) associated with the manuscript.


New feature : msPIP
---------------------

We applied the PIP framework to a DDA dataset.The dataset consists of eight experimental condition, each with three replicates (total of 24 runs). Twelve non-human proteins were spiked at known concentrations into constant HEK-293 background. 
We examined proportion of missing peptides per run before and after PIP. The volcano plots represent data for comparing group 8 vs group 1. 

**PIP reduces the proportion of missing values substantially, almost to zero.**

Figure: The proportion of missing peptides per sample in PASS00589 DDA dataset before and after PIP. 

<img src="https://user-images.githubusercontent.com/7257233/121839424-5c3a4380-cd1d-11eb-84fa-437a387c44f2.png" width="700px" align="center">


**PIP recovers the low abundance peptides and re-constructs the true dynamic range**

Low-abundance peptides not quantified by MaxQuant are recovered, and differential abundance results are improved. Note down regulated peptides that are not present in the volcano plot of DE test on MQ-reported data (bottom left), that are recovered by PIP (bottom right volcano plot) for the same experimental contrast.

<img src="https://user-images.githubusercontent.com/7257233/121839859-55600080-cd1e-11eb-998e-f7e60896b1bf.png" width="700px" align="center">


The PIP workflow involves the following two function calls:

```{r}
dda_pip <- mspip("/path/to/combined/txt", k=3, thresh = 0.0, tims_ms = FALSE, skip_weights = FALSE)
y_pip <- evidenceToMatrix(dda_pip, return_EList = TRUE)
```
Test for differential abundance in *limma*:

```{r}
y_pip <- normalizeBetweenArrays(y_pip, method = "quantile")
design <- model.matrix(~ group)
fit <- lmFit(y_pip, design)
fit <- eBayes(fit)
summary(decideTests(fit))
```
*limma* automatically recognizes the `EListRaw` object created by `evidenceToMatrix`, applies log2 transformation to intensity
values, and passes the PIP confidence scores as observation-level weights to `lmFit`. 


Need more help to start? Please see documentation. We have also collected a number of **case studies** [here]()

**Questions?** Please consider openning an issue.


Reference
-----------
```
@article{hediyeh2023msimpute,
  title={MsImpute: Estimation of missing peptide intensity data in label-free quantitative mass spectrometry},
  author={Hediyeh-Zadeh, Soroor and Webb, Andrew I and Davis, Melissa J},
  journal={Molecular \& Cellular Proteomics},
  pages={100558},
  year={2023},
  publisher={Elsevier}
}
```

