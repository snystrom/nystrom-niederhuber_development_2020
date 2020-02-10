# Expression of E93 provides an instructive cue to control dynamic enhancer activity and chromatin accessibility during development
## Analysis Code
## Spencer Nystrom

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3660090.svg)](https://doi.org/10.5281/zenodo.3660090)

## To rerun code
First run `R_Pipeline/00_Runall.R`

Next run `Published_Analysis.Rmd` which will source the output of the pipeline.

## Motif Analysis Results
All motif analysis was done using custom-built scripts for running MEME suite
tools from within R. These can be found in the `Scripts/Nystlib/R/MEME_Utils.R`
script. This functionality is currently being ported to an R package
[DremeR](https://github.com/snystrom/dremeR).

I used meme 4.12.0 for DREME analysis, and AMEv5.1.0 for motif scanning. I have
saved the R data results from the AME and DREME runs in `Data/` which are by
default loaded in the analysis notebook. To repeat the meme-suite runs
yourself, you'll have to manually remove the `eval=F` flags in the
corresponding notebook chunks. Note that using "shuffled" input for dreme runs
in this version breaks when setting a random seed, such that repeat iterations
using the same random seed will produce different results. This does not affect
overall conclusions of the manuscript. See notebook for details.

## Raw Data availability
Raw data will be available at [GSE141738](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141738) soon. This repo will also provide further description of the raw & processed data files for increased clarity.

## Session Info
Session info for my analysis can be found [here](sessionInfo.txt).
