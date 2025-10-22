#############################################
####                                     ####
####      GEMINI,                        ####
####      "direct" rg - BMI              ####
####                                     ####
####                                     ####
####      Ninon Mounier                  ####
####      04/10/2023                     ####
####                                     ####
#############################################

# assumed you are: 
#  1. running in the repo directory
#  2. have downloaded the "GEMINI GWAS" to directory GWASs_GEMINI/ 
#     get from Zenodo https://doi.org/10.5281/zenodo.14284046
#  3. downloaded the BMI summary statstics to a provided directory
#  4. have ran the MR script

# load packages

# remotes::install_github("n-mounier/bGWAS")
library(bGWAS)

# remotes::install_github("GEMINI-multimorbidity/partialLDSC")
library(partialLDSC) # version 0.1.1

library(tidyverse)
library(openxlsx)
library(dplyr)
library(tidyr)

#### create direct summary stats for all conditions ####
# could use an array job on isca but only takes a few minutes for each condition...

# use conditions for which we have evidence of BMI having a causal effect
MR_results = readr::read_tsv("_results/MR_BMI2015_070723.tsv") %>%
  filter(method == "Inverse variance weighted") 

partialrg_results = readr::read_tsv("_results/71LTCs_070823_200b_difference.tsv") 
partialrg_results = partialrg_results %>%
  mutate(diff.fdr = p.adjust(diff.P, method = "fdr")) 

partialrg_results %>%
  filter(diff.fdr < 0.05) ->partialrg_results_signif
conditions = unique(c(partialrg_results_signif$condition.1, partialrg_results_signif$condition.2))

MR_results  %>% 
  filter(outcome_phenocode %in% conditions) %>% filter(pval < 0.05) %>% nrow


# 23 conditions
MR_results  %>% 
       filter(outcome_phenocode %in% conditions) %>% filter(pval < 0.05/length(conditions)) -> MR_strong

# check if significant difference for these pairs
pairs = expand.grid(x = MR_strong$outcome_phenocode, y = MR_strong$outcome_phenocode, stringsAsFactors = F) %>%
  mutate(pair = case_when(
    x < y ~ paste0(x, "-", y),
    TRUE ~ paste0(y, "-", x)
  )) %>%
  filter(!duplicated(pair), x!=y)

partialrg_results_signif %>% mutate(pair = case_when(
  condition.1 < condition.2 ~ paste0(condition.1, "-", condition.2),
  TRUE ~ paste0(condition.2, "-", condition.1)
)) -> partialrg_results_signif
  
table(pairs$pair %in% partialrg_results_signif$pair)
# FALSE  TRUE 
#     7   246

conditions_bGWAS = MR_strong$outcome_phenocode

# user needs to produce ZMatrices 
my_Zmat = "/slade/home/nm572/Data/ZMatrices/"
my_studies = select_priorGWASs(include_traits=c("Body Mass Index"), Z_matrices = my_Zmat)  

# if folder does not exists, create _data/bGWAS
if(!dir.exists("bGWAS")) dir.create("bGWAS")

get_directSS <- function(condition){
  
  # load GEMINI GWAS file 
  data = data.table::fread(paste0("GWASs_GEMINI/GEMINI.sumstats.", condition, ".txt.gz"))

  # run bGWAS to remove causal effect of BMI
  bGWAS_res = bGWAS(name=condition,
                    GWAS = data,
                    prior_studies = my_studies, 
                    MR_threshold = 5e-8,
                    # MR_pruning_dist = 1000,
                    # MR_pruning_LD = 0.001,
                    Z_matrices = my_Zmat, 
                    stepwise_threshold = 0.05,
                    verbose = F)
  # each analysis takes a 2-4 minutes

  # save result file as RDS
  saveRDS(bGWAS_res, paste0("bGWAS/", condition, "_bGWAS.RDS"))

  if(!any(str_detect(bGWAS_res$log_info, "failed"))){

  # save observed/direct effects in LDSC format (make sure to align SNPs!)
  bGWAS_res_joined = full_join(data, bGWAS_res$all_BFs, by=c("SNP"="rsid"))
  aligned = bGWAS_res_joined$SNP[which(bGWAS_res_joined$A1 == bGWAS_res_joined$alt & bGWAS_res_joined$A2 == bGWAS_res_joined$ref)]
  swapped = bGWAS_res_joined$SNP[which(bGWAS_res_joined$A1 == bGWAS_res_joined$ref & bGWAS_res_joined$A2 == bGWAS_res_joined$alt)]
  
  bGWAS_res_joined %>%
    transmute(SNP,
              A1,
              A2,
              N,
              Z = case_when(
                SNP %in% aligned ~ z_direct,
                SNP %in% swapped ~ -z_direct,
                TRUE ~ NA_real_)) %>%
    filter(!is.na(Z)) -> bGWAS_direct
  data.table::fwrite(bGWAS_direct, paste0("bGWAS/", condition, "_direct.tsv.gz"),sep = "\t", quote = F, row.names = F )
  }
}


for(my_condition in conditions_bGWAS){
  print(my_condition)
  get_directSS(my_condition)
}


# once all files created, re-estimate rg 
# (use partial_ldsc script to make sure the SE are estimated correctly,
# even if we don't need the partial estimates)


# ld files
ld = "/path/to/LDSC/eur_w_ld_chr/"


conditions = paste0("bGWAS/", conditions_bGWAS, "_direct.tsv.gz")
condition.names = conditions_bGWAS


log_name <- "_results/23LTCs_direct_111023_200b"
n_blocks <- 200

# BMI file
confounder = "/path/to/GWAS/BMI_Yengo_2018.txt.sumstats.gz"
confounder.name = "BMI"

res = partial_ldsc(conditions = conditions,
                   confounder = confounder,
                   ld = ld, 
                   condition.names = condition.names, 
                   confounder.name = confounder.name,
                   n.blocks = n_blocks,
                   log.name = log_name)



### NOT NEEDED ####

# reformat results file to have partial / direct in a single file for the 253 pairs
res_table = readr::read_tsv("_results/23LTCs_direct_111023_200b_difference.tsv") %>% 
  mutate(pair = case_when(
    x < y ~ paste0(x, "-", y),
    TRUE ~ paste0(y, "-", x)
  )) %>% transmute(pair, x, y, direct_rg = rg, direct_rg.SE = rg.SE, 
                   direct_partial_rg = partial_rg, direct_partial_rg.SE = partial_rg.SE, 
                   direct_diff.P = diff.P, direct_diff.fdr = p.adjust(direct_diff.P, method = "fdr"))

partialrg_results %>% mutate(pair = case_when(
  x < y ~ paste0(x, "-", y),
  TRUE ~ paste0(y, "-", x)
)) -> partialrg_results

left_join(res_table, partialrg_results) -> all_res

readr::write_tsv(all_res, "_results/23LTCs_direct_111023_200b_difference_formatted.tsv")
