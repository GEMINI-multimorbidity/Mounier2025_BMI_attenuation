#############################################
####                                     ####
####      GEMINI, run MR                 ####
####      to describe BMI role           ####
####                                     ####
####                                     ####
####      Ninon Mounier                  ####
####      31/05/2023                     ####
####                                     ####
#############################################

## to launch on linux:
# Rscript 03_run_MR.R 2015 070723

# assumed you are: 
#  1. running in the repo directory
#  2. have downloaded the "GEMINI GWAS" to directory GWASs_GEMINI/ 
#     get from Zenodo https://doi.org/10.5281/zenodo.14284046

# get arguments
BMI <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
my_date <- commandArgs(trailingOnly = TRUE)[2]

# load packages
library(TwoSampleMR)
library(tidyr)
library(dplyr)

## exposure data
if(BMI == 2018){
  # Yengo et al. (2018)
  exp_BMI <- TwoSampleMR::extract_instruments("ieu-b-40")
  nrow(exp_BMI) # 507 IVs
  mean(exp_BMI$beta.exposure^2/exp_BMI$se.exposure^2) # Fstat = 72.8377
} else if(BMI == 2015){
  # Locke et al. (2015) - European only
  exp_BMI <- TwoSampleMR::extract_instruments("ieu-a-835")
  nrow(exp_BMI) # 69 IVs
  mean(exp_BMI$beta.exposure^2/exp_BMI$se.exposure^2) # Fstat = 66.77919
} else {
  stop("BMI should be 2015 or 2018")
}

## condition names
conditions = readr::read_tsv("_data/index_gemini.txt") |> dplyr::select(file_prefix) |> dplyr::pull()

## perform MR for each of the outcomes
all_results = data.frame()

for(i in 1:length(conditions)){
  out = conditions[i]
  print(paste0(i, " : ", out))
  
  # load outcome data
  out_data = data.table::fread(paste0("GWASs_GEMINI/GEMINI.sumstats.", out, ".txt.gz"),
                               data.table = F)
  
  # format outcome data for TwoSampleMR
    out_data %>%
      mutate(P = 10^-neg_log_10_p_value) -> out_data
    out_formatted <- TwoSampleMR::format_data(out_data,
                                              type = "outcome",
                                              snps = exp_BMI$SNP,
                                              header = TRUE,
                                              snp_col = "rs_id",
                                              beta_col = "beta",
                                              se_col = "standard_error",
                                              eaf_col = "effect_allele_frequency",
                                              effect_allele_col = "effect_allele",
                                              other_allele_col = "other_allele",
                                              pval_col = "P")
  
  # harmonise exposure / outcome data
  data <- TwoSampleMR::harmonise_data(exposure_dat = exp_BMI, 
                                      outcome_dat = out_formatted)
  
  # perform MR
  mr_results <- mr(data)
  mr_results %>%
    transmute(exposure = "BMI",
              outcome_phenocode = out,
              method, nsnp, b, se, pval) -> mr_results
  
  # combine results for all outcomes
  all_results = rbind(all_results, mr_results)
  
  
}

## write all results
write.table(all_results, paste0("_results/MR_BMI", BMI, "_", my_date,".tsv"), sep = "\t", quote = F, row.names = F)

