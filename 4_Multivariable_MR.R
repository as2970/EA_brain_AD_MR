# Setup
#------------------------------------------------------------------------------#

library(tidyverse)
library(data.table)
library(ivpack)
library(meta)
library(devtools)
library(pacman)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(phenoscanner)
library(LDlinkR)
library(mr.raps)
library(MRPRESSO)
library(extrafont)
library(anchors)

# Functions
#------------------------------------------------------------------------------#

MVMR_prep <- function(IDP) {
  
  # Extracting instruments from both exposure GWAS
  exposure_dat <- mv_extract_exposures_local(
    filenames_exposure = c(paste0(deparse(substitute(IDP)), "_exposure.txt", sep = ""), "EA_combined.txt"), 
    sep = " ",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "AF1",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    ncase_col = "N_case",
    ncontrol_col = "N_control",
    samplesize_col = "N",
    min_pval = 1e-200, 
    log_pval = FALSE, 
    pval_threshold = 5e-8, 
    clump_r2 = 0.001,
    clump_kb = 10000,
    harmonise_strictness = 2
  )

  
  # Printing number of IVs for exposure
  print(paste0("Number of IVs: ", as.character(length(exposure_dat$SNP)/2)))
  
  # Extracting instruments from outcome GWAS
  outcome_dat <- read_outcome_data(
    filename = "LOAD_outcome.txt", 
    snps = exposure_dat$SNP, 
    sep = " ",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "AF1",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    ncase_col = "N_case",
    ncontrol_col = "N_control",
    samplesize_col = "N",
    min_pval = 1e-200, 
    log_pval = FALSE, 
    chr_col = "CHR",
    pos_col = "POS"
  )
  
  # Harmonising exposure and outcome datasets
  dat <- mv_harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = outcome_dat, 
    harmonise_strictness = 2
  )
  
  # Creating a dataframe with data for multivariable MR
  mvmr_dat <- data.frame(
    "rsid" = rownames(dat[["exposure_beta"]]),
    "IDP_beta" = dat[["exposure_beta"]][ , 2],
    "IDP_se" = dat[["exposure_pval"]][ , 2], 
    "IDP_p" = dat[["exposure_se"]][ , 2], 
    "EA_beta" = dat[["exposure_beta"]][ , 1],
    "EA_se" = dat[["exposure_pval"]][ , 1],
    "EA_p" = dat[["exposure_se"]][ , 1],
    "LOAD_beta" = dat[["outcome_beta"]],
    "LOAD_se" = dat[["outcome_se"]],
    "LOAD_p" = dat[["outcome_pval"]]
  )
  rownames(mvmr_dat) <- c()

}


MVMR_analysis <- function(IDP) {
  
  load(paste0(deparse(substitute(IDP)), "_EA_LOAD.RData"))
  
  #identifier for the genetic variants
  rs <- mvmr_dat$rsid
  
  #summary data with LOAD
  beta_LOAD <- mvmr_dat$LOAD_beta
  se_LOAD <- mvmr_dat$LOAD_se
  
  #summary data with educational attainment and the IDP
  beta_EA <- mvmr_dat$EA_beta
  beta_IDP <- mvmr_dat$IDP_beta
  se_EA <- mvmr_dat$EA_se
  se_IDP <- mvmr_dat$IDP_se
  
  # Creating an MVMR object
  p_unload(TwoSampleMR)
  library(MendelianRandomization)
  
  mvmr_obj <- mr_mvinput(bx = cbind(beta_EA, beta_IDP), 
                         bxse = cbind(se_EA, se_IDP), 
                         by = beta_LOAD, 
                         byse = se_LOAD,
                         exposure = c("EA", deparse(substitute(IDP))),
                         outcome = "LOAD", 
                         snps = rs
  )
  
  # Performing a multivariable MR analysis considering the educational attainment and the IDP of interest
  mvmr_res = mr_mvivw(mvmr_obj)
  
  mvmr_res@Estimate = exp(mvmr_res@Estimate)
  mvmr_res@CILower = exp(mvmr_res@CILower)
  mvmr_res@CIUpper = exp(mvmr_res@CIUpper)
  
  df = data.frame(exposure = mvmr_res@Exposure, 
                  Estimate = paste0(sprintf("%.2f", round(mvmr_res@Estimate, 2)), " (", sprintf("%.2f", round(mvmr_res@CILower, 2)), ", ", sprintf("%.2f", round(mvmr_res@CIUpper, 2)), ")"),
                  P1 = format(mvmr_res@Pvalue, digits = 3, scientific = TRUE),
                  P2 = sprintf("%.3f", round(mvmr_res@Pvalue, 3))
  )
  
  rownames(df) = c()
  
  return(df)
}

#------------------------------------------------------------------------------#

# Data preparation
#------------------------------------------------------------------------------#

MVMR_prep(SA) #checked
MVMR_prep(vol) #checked
MVMR_prep(CT) #checked
MVMR_prep(LGI) #checked
MVMR_prep(MC) #checked
MVMR_prep(IC) #checked
MVMR_prep(FA_cortical) #no instruments for FA_cortical 
MVMR_prep(MD_cortical) #checked
MVMR_prep(ICVF) #checked
MVMR_prep(ODI) #checked
MVMR_prep(FA_WM) #checked
MVMR_prep(MD_WM) #checked
MVMR_prep(hippocampus_vol_left) #checked
MVMR_prep(hippocampus_vol_right) #checked
MVMR_prep(WMH_vol) #checked

# Run MVMR analysis 
#------------------------------------------------------------------------------#

mvmr_table = rbind(
  MVMR_analysis(SA),
  MVMR_analysis(vol),
  MVMR_analysis(CT),
  MVMR_analysis(LGI),
  MVMR_analysis(MC),
  MVMR_analysis(IC),
  MVMR_analysis(MD_cortical),
  MVMR_analysis(ICVF),
  MVMR_analysis(ODI),
  MVMR_analysis(FA_WM),
  MVMR_analysis(MD_WM),
  MVMR_analysis(hippocampus_vol_left),
  MVMR_analysis(hippocampus_vol_right),
  MVMR_analysis(WMH_vol)
)


