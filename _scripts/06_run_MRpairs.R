#############################################
####                                     ####
####      GEMINI, run MR (pairs)         ####
####      to describe BMI role           ####
####                                     ####
####                                     ####
####      Ninon Mounier                  ####
####      20/11/2023                     ####
####                                     ####
#############################################

# assumed you are: 
#  1. running in the repo directory
#  2. have downloaded the "GEMINI GWAS" to directory GWASs_GEMINI/ 
#     get from Zenodo https://doi.org/10.5281/zenodo.14284046
#  3. ran the previous scripts

## to be launched on a linux command line
# Rscript 06_runMR_pairs.R 2015 201123

# load packages
library(tidyverse)
library(TwoSampleMR)
library(tidyverse)


BMI <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
my_date <- commandArgs(trailingOnly = TRUE)[2]

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

## pairs, list of outcomes
# need beta/se, so can't use munged files
pairs = readr::read_tsv("_data/ncases_15pairs.tsv") %>%
  mutate(name=paste0(x, "-", y))


## perform MR for each of the outcomes
all_results = data.frame()

for(i in 1:nrow(pairs)){
  out = pairs$name[i]
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
write.table(all_results, paste0("_results/MR_BMI", BMI, "_pairs_", my_date,".tsv"), sep = "\t", quote = F, row.names = F)


nice_results = all_results %>% filter(method=="Inverse variance weighted")

# add prevalence
# UKB
n_UKB = 502643
nice_results = full_join(nice_results,
                         pairs %>%
                           select(outcome_phenocode = name,
                                  ncases_UKB = ncases_xy) %>%
                           mutate(prev_UKB = ncases_UKB/n_UKB))
                         

# GWAS

get_sample_size <- function(my_file){
  D =  read_lines(my_file)
  return(D[str_detect(D, "cases")])
}

# transform this into a data.frame
# condition - ncases - ncontrols  - ntotal - prevalence

pairs %>%
  mutate(file = paste0("./", name, "/", name, "_REGENIE.log")) -> GWASs

GWASs$file %>%
  map(get_sample_size) %>%
  reduce(rbind) -> samplesize
nice_results = full_join(nice_results,
                         as.data.frame(samplesize) %>%
                           # get condition names
                           separate(V1, c("before", "condition", "after"), sep = "'", remove=T) %>%
                           # get sample sizes
                           separate(after, c("n1", "ncases", "n2", "n3", "ncontrols", "n4"), sep = " ", remove=T) %>%
                           # clean columns
                           transmute(outcome_phenocode = pairs$name,
                                     ncases_GWAS = as.numeric(ncases),
                                     ncontrols_GWAS = as.numeric(ncontrols),
                                     neff_GWAS = 2/(1/ncases_GWAS+1/ncontrols_GWAS),
                                     prev_GWAS = ncases_GWAS/(ncases_GWAS+ncontrols_GWAS)))
    

nice_results %>%
  mutate(R = 1000 * (prev_GWAS - (b * prev_GWAS)/(1 - prev_GWAS + b * prev_GWAS))) -> nice_results


write.table(nice_results, paste0("_results/BMIintervention", BMI, "_pairs_", my_date,".tsv"), sep = "\t", quote = F, row.names = F)


                          