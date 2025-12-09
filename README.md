<h1 align="center"><project-name></h1>

<p align="center"><project-description></p>

## Introduction

This repository is intended as an online supplement to the manuscript:

Jackson P. Lautier. Vladimir Pozdnyakov. Jun Yan. "Estimating the time-to-event distribution for loan-level data within a consumer auto loan asset-backed security." Ann. Appl. Stat. 19 (4) 2831 - 2851, December 2025. [https://doi.org/10.1214/25-AOAS2103](https://doi.org/10.1214/25-AOAS2103)

Please attribute any citations of this repository to the original
manuscript.

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