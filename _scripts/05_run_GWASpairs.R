#############################################
####                                     ####
####      GEMINI,                        ####
####      GWAS - pairs                   ####
####                                     ####
####                                     ####
####      Ninon Mounier                  ####
####      14/11/2023                     ####
####                                     ####
#############################################




library(tidyverse)

setwd("/slade/projects/GEMINI/BMI_attenuation/")

library(tidyverse)

library(openxlsx)

library(rslurm)


# get the pairs, use already formatted results file
all_res = openxlsx::read.xlsx("_results/_tables/71LTCs_070823_200b_BMI.xlsx")

pairs = all_res %>%
  arrange(diff.FDR) %>%
  filter(diff.FDR < 0.05, rg.FDR < 0.05, partial_rg.FDR > 0.05) %>%
  slice(1:15)


# first, create pheno files for the 15 pairs for GWAS ####

if(!dir.exists("_data/GWASs_pairs")) dir.create("_data/GWASs_pairs")
setwd("_data/GWASs_pairs")

GWAS_info = readr::read_tsv("/slade/projects/GEMINI/UKBiobank/UKB_GWAS_version.txt",
                            col_names = c("condition", "version"))

pairs %>%
  transmute(x,y,
            ncases_x = NA_real_,
            ncases_y = NA_real_,
            ncases_xy = NA_real_) -> pairs


for(i in 1:nrow(pairs)){
  c1 = pairs$x[i]
  c2 = pairs$y[i]
  pair = paste0(c1, "-", c2)
  dir.create(paste0("_data/GWASs_pairs/", pair))
  
  # read phenofiles
  pheno_c1 = readr::read_delim(paste0("/slade/projects/GEMINI/UKBiobank/regenie_GWAS/",
                             c1, "/", c1, "_", 
                             GWAS_info$version[GWAS_info$condition==c1], ".pheno"))
  
  pheno_c2 = readr::read_delim(paste0("/slade/projects/GEMINI/UKBiobank/regenie_GWAS/",
                                    c2, "/", c2, "_", 
                                    GWAS_info$version[GWAS_info$condition==c2], ".pheno"))
  
  
  # merge to create new phenofile
  pheno_pair = full_join(pheno_c1, pheno_c2)
  
  n_cases = c(sum(pheno_c1[3]), sum(pheno_c2[3]))
  
  pheno_pair %>%
    select(FID, IID,
           cc1 = as.name(c1), cc2 = as.name(c2)) %>%
    # do not use {{pair}} := as PLINK does not handle well "-" in pheno names
    mutate(pheno := case_when(
      cc1 == 1 & cc2 == 1 ~ 1,
      TRUE ~ 0
    )) %>%
    select(FID, IID, pheno) -> pheno_pair

  n_cases = c(n_cases, sum(pheno_pair[3]))
  pairs[i,3:5] = n_cases
  names(n_cases) = c(c1, c2, pair)
  
  readr::write_delim(pheno_pair,
                     paste0("_data/GWASs_pairs/", pair, "/", pair, ".phenofile"))

  readr::write_rds(n_cases,
                   paste0("_data/GWASs_pairs/", pair, "/", pair, "_ncases.RDS"))
} 

readr::write_tsv(pairs,
                 "_data/GWASs_pairs/ncases_15pairs.tsv")


# then launch rslurm jobs for GWASs for each pair ####

# main function, for a pair
# use functions in REGENIE.R file
# would have been faster to analyse all pairs within a single pheno file...
run_GWAS <- function(pair){
  source("/slade/projects/GEMINI/BMI_attenuation/_scripts/REGENIE.R")
  setwd(paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  
  # use our own version of get_job_status
  get_job_status <- function (slr_job) {
    if (!(class(slr_job) == "slurm_job"))
      stop("input must be a slurm_job")
    stat <- suppressWarnings(system(paste("squeue -n", slr_job$jobname),
                                    intern = TRUE))
    if (length(stat) > 1) {
      res = "Job running or in queue."
    }
    else {
      res = "Job completed or stopped."
    }
    return(res)
  }
  # clean (_rslurm folder + .RDS files)
  clean <- function(slr_job, wd, keep_folder=T){
    # _rslurm folder
    if(!keep_folder) system(paste0("rm -rf  _rslurm_", slr_job$jobname))
    # .RDS files (in working directory)
    system(paste0("rm ", wd, "/", "results_*.RDS"))
  }
  
  ## step 0
  my_slurm_options_step0 = list(partition = "mrcq",
                                time = "0-00:15:00",
                                nodes = 1,
                                `ntasks-per-node` = 16,
                                mem = "250G",
                                account = "Research_Project-MRC158833",
                                `mail-type` = "END",
                                `mail-user` = "n.mounier@exeter.ac.uk")
  
  my_params=list(phenofile =  paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/", pair, ".phenofile"),
                 phenoname = "pheno",
                 output = pair,
                 wd = paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  
  step0 = slurm_call(REGENIE_exclusion,
                      params = my_params,
                      jobname = paste0("step0_", pair),
                      slurm_options = my_slurm_options_step0,
                      submit = TRUE)
  
  while(get_job_status(step0)!="Job completed or stopped."){
    Sys.sleep(10)
  } 
  clean(step0, paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  
  
  ## step 1 
  my_slurm_options_step1 = list(partition = "mrcq",
                                time = "0-10:00:00",
                                nodes = 1,
                                `ntasks-per-node` = 16,
                                mem = "250G",
                                account = "Research_Project-MRC158833",
                                `mail-type` = "END",
                                `mail-user` = "n.mounier@exeter.ac.uk")
  
  my_params =list(snplist =  paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/", pair, ".snplist"),
                  phenofile =  paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/", pair, ".phenofile"),
                  output = paste0(pair, "_REGENIE"),
                  binary = T,
                  wd = paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  step1 = slurm_call(REGENIE_step1,
                     params = my_params,
                     jobname =  paste0("step1_", pair),
                     slurm_options = my_slurm_options_step1,
                     submit = TRUE)
  
  
  while(get_job_status(step1)!="Job completed or stopped."){
    Sys.sleep(100)
  } 
  clean(step1, paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  
  
  ## step 2
  my_slurm_options_step2 = list(partition = "mrcq",
                                time = "0-120:00:00",
                                nodes = 1,
                                `ntasks-per-node` = 4,
                                mem = "64G",
                                account = "Research_Project-MRC158833",
                                `mail-type` = "END",
                                `mail-user` = "n.mounier@exeter.ac.uk")
  my_params = data.frame(chr=1:22,
                         predlist = paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/", pair, "_REGENIE_pred.list"),
                         phenofile =  paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/", pair, ".phenofile"),
                         output = paste0(pair, "_REGENIE_OUT"),
                         binary = T,
                         wd = paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  step2 = slurm_apply(REGENIE_step2,
                      nodes = 22,
                      cpus_per_node = 1,
                      params = my_params,
                      jobname = paste0("step2_", pair),
                      slurm_options = my_slurm_options_step2,
                      submit = TRUE)
  while(get_job_status(step2)!="Job completed or stopped."){
    Sys.sleep(100)
  } 
  clean(step2, paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  

  # QC
  clean_SummaryStats <- function(wd){
    init_wd = getwd()
    on.exit(setwd(init_wd))
    
    setwd(wd)
    name = stringr::str_split(wd, "/", simplify = T)[length(stringr::str_split(wd, "/", simplify=T))]
    
    # get name of trait
    trait =  stringr::str_remove(
        system("ls *.phenofile", intern = T),
        ".phenofile")
    
    ## rename CHR1 to a new file
    file.copy(paste0(trait, "_REGENIE_OUT_chr1_pheno.regenie"), paste0(trait, "_REGENIE.txt"))
    
    ## for the rest of the chromosomes append! use tail +2 to skip headers
    for(chr in 2:22){
      system(paste0("tail -n +2 ", trait, "_REGENIE_OUT_chr", chr, "_pheno.regenie >> ", trait, "_REGENIE.txt"))
    }
      
      
      
    ## replace spaces with a single tab to make it a tab-separated file + subset by A1freq and info
    GWAS_tmp <- data.table::fread(paste0(trait, "_REGENIE.txt")) # 71,520,913 
    data.table::fwrite(GWAS_tmp, paste0(trait, "_REGENIE.txt.gz"), sep="\t")
    
    file.remove(paste0(trait, "_REGENIE.txt"))
    GWAS_tmp %>%
      filter(INFO >= 0.3,
             A1FREQ > 0.001,
             A1FREQ < 0.999) -> GWAS_tmp # 16,468,335 (same as in Beth's QC'ed files)
    data.table::fwrite(GWAS_tmp, paste0(trait, "_REGENIE.txt.qc.gz"), sep="\t")
    
    
  }
  my_slurm_options_QC = list(partition = "mrcq",
                             time = "0-10:00:00",
                             nodes = 1,
                             `ntasks-per-node` = 1,
                             mem = "40G",
                             account = "Research_Project-MRC158833",
                             `mail-type` = "END",
                             `mail-user` = "n.mounier@exeter.ac.uk")
  my_params =list(wd = paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  QC = slurm_call(clean_SummaryStats,
                  params = my_params,
                  jobname = paste0("QC_", pair),
                  slurm_options = my_slurm_options_QC,
                  submit = TRUE)
  while(get_job_status(QC)!="Job completed or stopped."){
    Sys.sleep(100)
  } 
  clean(QC, paste0("/slade/projects/GEMINI/BMI_attenuation/_data/GWASs_pairs/", pair, "/"))
  
  
}

for(i in 1:nrow(pairs)){
  my_pair = paste0(pairs$x[i], "-", pairs$y[i])
  my_slurm_options = list(partition = "mrcq",
                          time = "0-128:00:00",
                          nodes = 1,
                          `ntasks-per-node` = 1,
                          mem = "2G",
                          account = "Research_Project-MRC158833",
                          `mail-type` = "END",
                          `mail-user` = "n.mounier@exeter.ac.uk")
  my_params = data.frame(pair = my_pair)
  GWAS = slurm_call(run_GWAS,
                    params = my_params,
                    jobname = paste0("GWAS_", my_pair),
                    slurm_options = my_slurm_options,
                    submit = TRUE)
}
