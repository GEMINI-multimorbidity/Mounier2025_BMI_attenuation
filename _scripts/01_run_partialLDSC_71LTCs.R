#############################################
####                                     ####
####      GEMINI,                        ####
####      partial rg - BMI               ####
####                                     ####
####                                     ####
####      Ninon Mounier                  ####
####      07/08/2023                     ####
####                                     ####
#############################################

# assumed you are: 
#  1. running in the repo directory
#  2. have downloaded the "GEMINI GWAS" to directory GWASs_GEMINI/ 
#     get from Zenodo https://doi.org/10.5281/zenodo.14284046
#  3. downloaded the BMI summary statstics to a provided directory

## to be launched on a linux command line
#     e.g., for BMI - 2018, different number of blocks
# Rscript 01_run_partialLDSC_71LTCs.R 71LTCs_070823_200b 200 /path/to/BMI_Yengo_2018.txt.sumstats.gz BMI

#remotes::install_github("GEMINI-multimorbidity/partialLDSC")
library(partialLDSC) # version 0.1.1

library(openxlsx)
library(dplyr)
library(tidyr)

# ld files
ld = "/path/to/LDSC/eur_w_ld_chr/"

# conditions files
conditions = readr::read_tsv("_data/index_gemini.txt") |> dplyr::select(file_prefix) |> dplyr::pull()


my_conditions = data.frame(condition = conditions) %>%
  filter(!conditions == "obesity") %>%
  mutate(folder = paste0("GWASs_GEMINI/",condition, "/"),
         file = paste0(folder, "/", condition, "_GEMINI.sumstats.gz"))

conditions = my_conditions$file
condition.names = my_conditions$condition


log_name <- commandArgs(trailingOnly = TRUE)[1]
n_blocks <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# BMI file
confounder = commandArgs(trailingOnly = TRUE)[3]
confounder.name = commandArgs(trailingOnly = TRUE)[4]


res = partial_ldsc(conditions = conditions,
                   confounder = confounder,
                   ld = ld,
                   condition.names = condition.names, 
                   confounder.name = confounder.name,
                   n.blocks = n_blocks,
                   log.name = log_name)
