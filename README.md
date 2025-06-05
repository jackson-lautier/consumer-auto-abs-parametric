<h1 align="center"><project-name></h1>

<p align="center"><project-description></p>

## NOTE: REVISIONS IN PROGRESS

As of Spring 2025, this repository is in a revision state.  Once complete, this
note will be removed.

## Introduction

This repository is intended as an online supplement to the manuscript,
_Estimating the time-to-event distribution for
loan-level data within a consumer auto loan asset-backed security_.
Please attribute any citations of this repository to the original
manuscript.


This repository includes:

- **data** Loan demographic and performance data from the ABS bonds AART 2017-3
and AART 2019-3.

- **code** Replication code files.  To replicate the numeric results and
simulation studies, use 's5-simulation-studies.R'.  To replicate the application
results, use 's6-application.R'.  Please note these replication files rely on
functions also stored within the **code** folder.

## Workflow

There are two workstreams: (1) the numeric and simulation
studies and (2) the data application.  The latter will require
calls to **data** and the former is self-contained.  All output,
including figures, will be stored in a new folder 'results'.

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