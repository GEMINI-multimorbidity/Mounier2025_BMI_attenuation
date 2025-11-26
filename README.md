# GEMINI - partial genetic correlation when adjusting for BMI genetics

<!-- badges: start -->
[![](https://img.shields.io/badge/version-1.0-informational.svg)](https://github.com/GEMINI-multimorbidity/partialLDSC)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17725332.svg)](https://doi.org/10.5281/zenodo.17725332)
<!-- badges: end -->

This repository contains the scripts for the paper:

 - Mounier et al. (2025) "Genetics identifies obesity as a shared risk factor for co-occurring multiple long-term conditions" Communications Medicine

If you use the scripts/data/results please cite the paper.

## Required data 

Required data not included in the repo are:

 - the GWAS results from GEMINI: https://doi.org/10.5281/zenodo.14284046
 - the BMI GWAS summary statistics from 2015 and 2018: https://giant-consortium.web.broadinstitute.org/GIANT_consortium_data_files

## Scripts

`01_run_partialLDSC_71LTCs.R` is used to get all the pairwise genetic correlations and create the [raw results file](_results/71LTCs_070823_200b_difference.tsv). 
Note that this script can be used with different number of blocks (for comparison), or other BMI data to generate the other raw results files.
It uses the `partial_ldsc.R` function from the [`LDSCpartial`](https://github.com/GEMINI-multimorbidity/partialLDSC) R-package.    

`02_process_partialLDSC_71LTCs.R` is used to create the [main table](_results/71LTCs_070823_200b.xlsx)  

`03_run_MR.R` is used to get the causal effect of BMI on the conditions using MR and create [MR_BMI2015_070723.tsv](_results/MR_BMI2015_070723.tsv) and [MR_BMI2018_070723.tsv](_results/MR_VMI2018_070723.tsv).)      

`04_run_bGWAS.R` is used to obtain direct effects using the [`bGWAS`](https://github.com/n-mounier/bGWAS) R-package for 23 conditions and re-estimate pairwise genetic correlations using these effects. 
Partial genetic correlation results are in [23LTCs_direct_111023_200b_difference.tsv](_results/23LTCs_direct_111023_200b_difference.tsv).    

`05_run_GWASpairs.R` is used to perform GWAS on 15 pairs of condition and relies on functions from `REGENIE.R` to perform GWASs using REGENIE. 

`06_run_MRpairs.R` is used to get the causal effect of BMI on the pairs using MR, and estimate the effect of intervening on BMI and create [BMIintervention2015_pairs_201123.tsv](_results/BMIintervention2015_pairs_201123.tsv)  

