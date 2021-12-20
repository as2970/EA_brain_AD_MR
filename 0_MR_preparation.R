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

MR_prep <- function(exp, out) {

  # Extracting instruments
  exp_dat <- read_exposure_data(
    filename = paste0(deparse(substitute(exp)), "_exposure.txt", sep = ""),
    clump = FALSE,
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
  
  # Clumping instruments
  exp_dat <- exp_dat %>% 
    rename(
      rsid = SNP,
      pval = pval.exposure
    )
  
  exp_dat_clumped <- ld_clump(
    dat = exp_dat,
    clump_kb = 10000, 
    clump_r2 = 0.001, 
    clump_p = 5e-8,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "" #path to LD reference dataset
  )
  
  exp_dat_clumped <- exp_dat_clumped %>% 
    rename(
      SNP = rsid,
      pval.exposure = pval
    )

  
  # Printing number of IVs for exposure
  print(paste0("Number of IVs: ", as.character(length(exp_dat_clumped$SNP))))
  
  # Extracting instruments from outcome GWAS
  out_dat <- read_outcome_data(
    filename = paste0(deparse(substitute(out)), "_outcome.txt", sep = ""), 
    snps = exp_dat_clumped$SNP, 
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
  
  # Identifying & printing exposure instruments missing from outcome GWAS
  missing_IVs <- exp_dat_clumped$SNP[!(exp_dat_clumped$SNP %in% out_dat$SNP)]
  print(paste0("Number of IVs missing from outcome GWAS: ", as.character(length(missing_IVs))))
  print("List of IVs missing from outcome GWAS:")
  for (i in 1:length(missing_IVs)) {
    print(paste0(missing_IVs[i]))
  }
  
  # Replacing missing instruments from outcome GWAS with proxies
  if(length(missing_IVs) == 0) {
    
    print("All exposure IVs found in outcome GWAS.")
    
  } else {
    
    print("Some exposure IVs missing from outcome GWAS.")
    out_full <- fread(paste0(deparse(substitute(out)), "_outcome.txt", sep = ""))
    
    for (i in 1:length(missing_IVs)) {
      
      proxies <- LDproxy(snp = missing_IVs[i], pop = "EUR", r2d = "r2", token = "6fb632e022ef", file = FALSE)
      proxies <- proxies[proxies$R2 > 0.8, ]
      proxy_present = FALSE
      
      if(length(proxies$RS_Number) == 0){
        
        print(paste0("No proxy SNP available for ", missing_IVs[i]))
        
      } else {
        
        for (j in 1:length(proxies$RS_Number)) {
          
          proxy_present <- proxies$RS_Number[j] %in% out_full$SNP
          
          if (proxy_present) {
            proxy_SNP = proxies$RS_Number[j]
            proxy_SNP_allele_1 = str_sub(proxies$Alleles[j], 2, 2)
            proxy_SNP_allele_2 = str_sub(proxies$Alleles[j], 4, 4)
            original_SNP_allele_1 = str_sub(proxies$Alleles[1], 2, 2)
            original_SNP_allele_2 = str_sub(proxies$Alleles[1], 4, 4)
            break
          }
        }
      }
      
      if(proxy_present == TRUE) {
        print(paste0("Proxy SNP found. ", missing_IVs[i], " replaced with ", proxy_SNP))
        proxy_row <- out_dat[1, ]
        proxy_row$SNP = missing_IVs[i]
        proxy_row$beta.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "BETA"])
        proxy_row$se.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "SE"])
        if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome = original_SNP_allele_1
        if (out_full[out_full$SNP == proxy_SNP, "A1"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome = original_SNP_allele_2
        if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome = original_SNP_allele_1
        if (out_full[out_full$SNP == proxy_SNP, "A2"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome = original_SNP_allele_2
        proxy_row$pval.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "P"])
        proxy_row$samplesize.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N"])
        if("N_case" %in% colnames(out_full)) proxy_row$ncase.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N_case"])
        if("N_control" %in% colnames(out_full))proxy_row$ncontrol.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "N_control"])
        proxy_row$chr.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "chr.exposure"])
        proxy_row$pos.outcome = as.numeric(exp_dat_clumped[exp_dat_clumped$SNP == missing_IVs[i], "pos.exposure"])
        if("AF1" %in% colnames(out_full)) proxy_row$eaf.outcome = as.numeric(out_full[out_full$SNP == proxy_SNP, "AF1"])
        out_dat <- rbind(out_dat, proxy_row)
      }
      
      if(proxy_present == FALSE) {
        print(paste0("No proxy SNP available for ", missing_IVs[i], " in outcome GWAS."))
      }
    }
    
  }
  
  # Harmonising exposure and outcome datasets
  dat <- harmonise_data(
    exposure_dat = exp_dat_clumped, 
    outcome_dat = out_dat, 
    action = 2
  )

}

# EA -> LOAD
#------------------------------------------------------------------------------#
MR_prep(exp = EA, out = LOAD)

# LOAD -> EA
#------------------------------------------------------------------------------#
MR_prep(exp = LOAD, out = EA)

# EA -> IDPs
#------------------------------------------------------------------------------#
MR_prep(exp = EA, out = SA)
MR_prep(exp = EA, out = vol)
MR_prep(exp = EA, out = CT)
MR_prep(exp = EA, out = IC)
MR_prep(exp = EA, out = LGI)
MR_prep(exp = EA, out = MC)
MR_prep(exp = EA, out = ICVF)
MR_prep(exp = EA, out = ODI)
MR_prep(exp = EA, out = FA_cortical)
MR_prep(exp = EA, out = MD_cortical)
MR_prep(exp = EA, out = FA_WM)
MR_prep(exp = EA, out = MD_WM)
MR_prep(exp = EA, out = hippocampus_vol_left)
MR_prep(exp = EA, out = hippocampus_vol_right)
MR_prep(exp = EA, out = WMH_vol)

# LOAD -> IDPs
#------------------------------------------------------------------------------#
MR_prep(exp = LOAD, out = SA)
MR_prep(exp = LOAD, out = vol)
MR_prep(exp = LOAD, out = CT)
MR_prep(exp = LOAD, out = IC)
MR_prep(exp = LOAD, out = LGI)
MR_prep(exp = LOAD, out = MC)
MR_prep(exp = LOAD, out = ICVF)
MR_prep(exp = LOAD, out = ODI)
MR_prep(exp = LOAD, out = FA_cortical)
MR_prep(exp = LOAD, out = MD_cortical)
MR_prep(exp = LOAD, out = FA_WM)
MR_prep(exp = LOAD, out = MD_WM)
MR_prep(exp = LOAD, out = hippocampus_vol_left)
MR_prep(exp = LOAD, out = hippocampus_vol_right)
MR_prep(exp = LOAD, out = WMH_vol)

# IDPs -> EA
#------------------------------------------------------------------------------#
MR_prep(exp = SA, out = EA)
MR_prep(exp = vol, out = EA)
MR_prep(exp = CT, out = EA)
MR_prep(exp = IC, out = EA)
MR_prep(exp = LGI, out = EA)
MR_prep(exp = MC, out = EA)
MR_prep(exp = ICVF, out = EA)
MR_prep(exp = ODI, out = EA)
MR_prep(exp = FA_cortical, out = EA)
MR_prep(exp = MD_cortical, out = EA)
MR_prep(exp = FA_WM, out = EA)
MR_prep(exp = MD_WM, out = EA)
MR_prep(exp = hippocampus_vol_left, out = EA)
MR_prep(exp = hippocampus_vol_right, out = EA)
MR_prep(exp = WMH_vol, out = EA)

# IDPs -> LOAD
#------------------------------------------------------------------------------#
MR_prep(exp = SA, out = LOAD)
MR_prep(exp = vol, out = LOAD)
MR_prep(exp = CT, out = LOAD)
MR_prep(exp = IC, out = LOAD)
MR_prep(exp = LGI, out = LOAD)
MR_prep(exp = MC, out = LOAD)
MR_prep(exp = ICVF, out = LOAD)
MR_prep(exp = ODI, out = LOAD)
MR_prep(exp = FA_cortical, out = LOAD)
MR_prep(exp = MD_cortical, out = LOAD)
MR_prep(exp = FA_WM, out = LOAD)
MR_prep(exp = MD_WM, out = LOAD)
MR_prep(exp = hippocampus_vol_left, out = LOAD)
MR_prep(exp = hippocampus_vol_right, out = LOAD)
MR_prep(exp = WMH_vol, out = LOAD)

