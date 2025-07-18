<h1 align="center"><project-name></h1>

<p align="center"><project-description></p>

## Introduction

This repository is intended as an online supplement to the manuscript,
_Estimating the time-to-event distribution for
loan-level data within a consumer auto loan asset-backed security_.
Please attribute any citations of this repository to the original
manuscript.  For reference, the current version of the manuscript may be
found [here](https://img1.wsimg.com/blobby/go/e126e6bc-09bb-4685-8200-323fe6a91322/downloads/a7266967-1ce1-410d-acd3-d83b26b41db9/aoas-main-07182025.pdf?ver=1752781971609) with
its supplement [here](https://img1.wsimg.com/blobby/go/e126e6bc-09bb-4685-8200-323fe6a91322/downloads/832237ae-9742-40d9-b60f-fb8a07effca2/aoas-supplement-07182025.pdf?ver=1752781971609).


This repository includes:

- **raw-data** Scraped loan demographic and performance data from the ABS bonds
AART 2017-3 and AART 2019-3.

- **data-clean** Cleaned raw data into files used within the manuscript.  These
files are identical to the files created by `data-processing.R' in the **code**
folder.

- **code** Replication code files.  First run `data-processing.R` to create the
clean data files in a new folder, **processed-data** (alternatively, rename the
**data-clean** folder as **processed-data**).  Second, all results in the
manuscript can be replicated with `data-analysis.R`.
All results will either print in the R console or be
stored in a new folder, **results**.


## Screenshots

![Asymptotic Normality](/illustrative-figures/sim_comps.pdf)

![AART Application](/illustrative-figures/aart_comp.pdf)

## Lead, Corresponding Author

**Jackson P. Lautier**

- [Website](https://jacksonlautier.com/)

## Complete Authors

**Vladimir Pozdnyakov**

- [Website](https://vladimir-pozdnyakov.github.io/)

**Jun Yan**

- [Website](http://merlot.stat.uconn.edu/~jyan/)