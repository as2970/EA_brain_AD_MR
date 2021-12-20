# Setup
#------------------------------------------------------------------------------#
library(tidyverse)
library(data.table)
library(pacman)
library(TwoSampleMR)

het_table = function(exp, out, exp_label, out_label) {
  load(paste0(deparse(substitute(exp)), "_", deparse(substitute(out)), ".RData"))
  dat = dat[dat$mr_keep == TRUE, ]
  het = mr_heterogeneity(dat)[ ,c("exposure", "outcome", "method", "Q", "Q_df", "Q_pval")]
  het$exposure = exp_label
  het$outcome = out_label
  return(het)
}


Steiger_table = function(exposure, outcome, exposure_label, outcome_label) {
  p_unload(MendelianRandomization)
  library(TwoSampleMR)
  
  load(paste0(deparse(substitute(exposure)), "_", deparse(substitute(outcome)), ".RData"))
  
  dat = dat[dat$mr_keep == TRUE, ]
  
  if (deparse(substitute(exposure)) == "LOAD") {
    
    dat = dat[!is.na(dat$ncase.exposure), ]
    dat$units.outcome = "SD"
    
    dat$units.exposure = "log odds"
    for (i in 1:length(dat$SNP)) {
      if (is.na(dat$eaf.exposure[i])) dat$eaf.exposure[i] = dat$eaf.outcome[i]
    }
    dat$prevalence.exposure = 0.05
    
  } else {
    if (deparse(substitute(outcome)) == "LOAD") {
      
      dat$units.exposure = "SD"
      
      dat$units.outcome = "log odds"
      for (i in 1:length(dat$SNP)) {
        if (is.na(dat$eaf.outcome[i])) dat$eaf.outcome[i] = dat$eaf.exposure[i]
      }
      dat$prevalence.outcome = 0.05
    } else {
      dat$units.exposure = "SD"
      dat$units.outcome = "SD"
    }
  }

  
  steiger_results = steiger_filtering(dat)
  dat_valid = steiger_results[steiger_results$steiger_dir == TRUE, ]
  
  p_unload(TwoSampleMR)
  library(MendelianRandomization)
  mr_object_all = mr_input(bx = dat$beta.exposure, 
                           bxse = dat$se.exposure, 
                           by = dat$beta.outcome, 
                           byse = dat$se.outcome, 
                           exposure = deparse(substitute(exposure)), 
                           outcome = deparse(substitute(outcome)), 
                           snps = dat$SNP)
  
  ivw_res_all = mr_ivw(mr_object_all)
  
  if (deparse(substitute(outcome)) == "LOAD") {
    ivw_res_all@Estimate = exp(ivw_res_all@Estimate)
    ivw_res_all@CILower = exp(ivw_res_all@CILower)
    ivw_res_all@CIUpper = exp(ivw_res_all@CIUpper)
  }
  
  mr_object_valid = mr_input(bx = dat_valid$beta.exposure, 
                             bxse = dat_valid$se.exposure, 
                             by = dat_valid$beta.outcome, 
                             byse = dat_valid$se.outcome, 
                             exposure = deparse(substitute(exp)), 
                             outcome = deparse(substitute(out)), 
                             snps = dat_valid$SNP)
  
  ivw_res_valid = mr_ivw(mr_object_valid)
  
  table = data.frame(exposure = exposure_label, 
                     outcome = outcome_label, 
                     `Total SNPs` = length(dat$SNP),
                     `Valid SNPs` = sum(steiger_results$steiger_dir), 
                     `Invalid SNPs` = sum(!steiger_results$steiger_dir), 
                     `MR-IVW (all SNPs)` = paste0(sprintf("%.2f", round(ivw_res_all@Estimate, 2)), " (", sprintf("%.2f", round(ivw_res_all@CILower, 2)), ", ", sprintf("%.2f", round(ivw_res_all@CIUpper, 2)), ")"),
                     `P-value (all SNPs)` = ivw_res_all@Pvalue,
                     `MR-IVW (valid SNPs)` = paste0(sprintf("%.2f", round(ivw_res_valid@Estimate, 2)), " (", sprintf("%.2f", round(ivw_res_valid@CILower, 2)), ", ", sprintf("%.2f", round(ivw_res_valid@CIUpper, 2)), ")"),
                     `P-value (valid SNPs)` = ivw_res_valid@Pvalue
  )
  
  return(table)
}

Egger_intercept_table = function(exp, out) {
  
  p_unload(TwoSampleMR)
  library(MendelianRandomization)
  
  load(paste0(deparse(substitute(exp)), "_", deparse(substitute(out)), ".RData"))
  
  dat = dat[dat$mr_keep == TRUE, ]
  
  mr_object = mr_input(bx = dat$beta.exposure, 
                       bxse = dat$se.exposure, 
                       by = dat$beta.outcome, 
                       byse = dat$se.outcome, 
                       exposure = deparse(substitute(exp)), 
                       outcome = deparse(substitute(out)), 
                       snps = dat$SNP)
  
  est_egger = mr_egger(mr_object)
  
  Intercept = sprintf("%.2f",round((est_egger@Intercept),2))
  CILower.Int = sprintf("%.2f",round((est_egger@CILower.Int),2))
  CIUpper.Int = sprintf("%.2f",round((est_egger@CIUpper.Int),2))
  
  df = data.frame("Exposure" = deparse(substitute(exp)),
                  "Outcome" = deparse(substitute(out)), 
                  "Egger_intercept" = paste0(Intercept, " (", CILower.Int, ", ", CIUpper.Int, ")"), 
                  "Pval" = est_egger@Pvalue.Int
                  )
  
  return(df)
}


SNPs_table = function(exp, out, exp_label, out_label) {
  load(paste0(deparse(substitute(exp)), "_", deparse(substitute(out)), ".RData"))
  table = dat[dat$mr_keep == TRUE, c("SNP", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "beta.outcome", "se.outcome", "pval.outcome")]
  table$`F-statistic` = (table$beta.exposure/table$se.exposure)^2
  
  table$eaf.exposure = sprintf("%.2f", round(table$eaf.exposure, 2))
  table$`F-statistic` = sprintf("%.2f", round(table$`F-statistic`, 2))
  
  
  table$beta.exposure = sprintf("%.3f", round(as.numeric(table$beta.exposure), 3))
  table$se.exposure = sprintf("%.3f", round(as.numeric(table$se.exposure), 3))
  table$beta.outcome = sprintf("%.3f", round(as.numeric(table$beta.outcome), 3))
  table$se.outcome = sprintf("%.3f", round(as.numeric(table$se.outcome), 3))
  
  
  table$pval.exposure = format(table$pval.exposure, digits = 3, scientific = TRUE)
  table$pval.outcome = format(table$pval.outcome, digits = 3, scientific = TRUE)

  for (i in 1:length(table$SNP)) {
    if (as.numeric(table$pval.exposure[i]) >= 0.001) table$pval.exposure[i] = sprintf("%.3f", round(as.numeric(table$pval.exposure[i]), 3))
    if (as.numeric(table$pval.outcome[i]) >= 0.001) table$pval.outcome[i] = sprintf("%.3f", round(as.numeric(table$pval.outcome[i]), 3))
  }
  
  table = table %>% rename(
    `Effect allele` = effect_allele.exposure, 
    `Other allele` = other_allele.exposure, 
    `Effect allele frequency` = eaf.exposure
  )

}

# Sensitivity analyses: MR-Egger intercept test
#------------------------------------------------------------------------------#
Egger_intercept_table(exp = EA, out = LOAD)
Egger_intercept_table(exp = EA, out = SA)
Egger_intercept_table(exp = EA, out = vol)
Egger_intercept_table(exp = EA, out = LGI)
Egger_intercept_table(exp = EA, out = IC)
Egger_intercept_table(exp = EA, out = ICVF)
Egger_intercept_table(exp = EA, out = WMH_vol)
Egger_intercept_table(exp = SA, out = EA)
Egger_intercept_table(exp = vol, out = EA)
Egger_intercept_table(exp = IC, out = EA)
Egger_intercept_table(exp = LOAD, out = ODI)
Egger_intercept_table(exp = LOAD, out = MD_WM)
Egger_intercept_table(exp = LOAD, out = hippocampus_vol_left)
Egger_intercept_table(exp = LOAD, out = hippocampus_vol_right)



# Sensitivity analyses: heterogeneity statistics
#------------------------------------------------------------------------------#
heterogeneity_table = rbind(
het_table(exp = EA, out = LOAD, exp_label = "Educational attainment", out_label = "Late-onset Alzheimer's disease"), 
het_table(exp = EA, out = SA, exp_label = "Educational attainment", out_label = "Surface area"), 
het_table(exp = EA, out = vol, exp_label = "Educational attainment", out_label = "Volume"), 
het_table(exp = EA, out = LGI, exp_label = "Educational attainment", out_label = "Local gyrification index"), 
het_table(exp = EA, out = IC, exp_label = "Educational attainment", out_label = "Intrinsic curvature"), 
het_table(exp = EA, out = ICVF, exp_label = "Educational attainment", out_label = "Intracellular volume fraction"), 
het_table(exp = EA, out = WMH_vol, exp_label = "Educational attainment", out_label = "Total volume of white matter hyperintensities"), 
het_table(exp = SA, out = EA, exp_label = "Surface area", out_label = "Educational attainment"), 
het_table(exp = vol, out = EA, exp_label = "Volume", out_label = "Educational attainment"), 
het_table(exp = IC, out = EA, exp_label = "Intrinsic curvature", out_label = "Educational attainment"), 
het_table(exp = LOAD, out = ODI, exp_label = "Late-onset Alzheimer's disease", out_label = "Orientation dispersion index"), 
het_table(exp = LOAD, out = MD_WM, exp_label = "Late-onset Alzheimer's disease", out_label = "Mean diffusivity (white matter tracts)"), 
het_table(exp = LOAD, out = hippocampus_vol_left, exp_label = "Late-onset Alzheimer's disease", out_label = "Volume of left hippocampus"), 
het_table(exp = LOAD, out = hippocampus_vol_right, exp_label = "Late-onset Alzheimer's disease", out_label = "Volume of right hippocampus")
)

heterogeneity_table$Q = sprintf("%.2f", round(heterogeneity_table$Q, 2))

for (i in 1:28) {
  if (as.numeric(heterogeneity_table$Q_pval[i]) >= 0.001) heterogeneity_table$Q_pval[i] = sprintf("%.3f", round(as.numeric(heterogeneity_table$Q_pval[i]), 3))
  else heterogeneity_table$Q_pval[i] = format(as.numeric(heterogeneity_table$Q_pval[i]), digits = 3, scientific = TRUE)
}


# Sensitivity analyses: MR Steiger directionality test
#------------------------------------------------------------------------------#
Steiger_results = rbind(
  Steiger_table(exposure = EA, outcome = LOAD, exposure_label = "Educational attainment", outcome_label = "Late-onset Alzheimer's disease"),
  Steiger_table(exposure = EA, outcome = SA, exposure_label = "Educational attainment", outcome_label = "Surface area"),
  Steiger_table(exposure = EA, outcome = vol, exposure_label = "Educational attainment", outcome_label = "Volume"),
  Steiger_table(exposure = EA, outcome = LGI, exposure_label = "Educational attainment", outcome_label = "Local gyrification index"),
  Steiger_table(exposure = EA, outcome = IC, exposure_label = "Educational attainment", outcome_label = "Intrinsic curvature"),
  Steiger_table(exposure = EA, outcome = ICVF, exposure_label = "Educational attainment", outcome_label = "Intracellular volume fraction"),
  Steiger_table(exposure = EA, outcome = WMH_vol, exposure_label = "Educational attainment", outcome_label = "Total volume of white matter hyperintensities"),
  Steiger_table(exposure = SA, outcome = EA, exposure_label = "Surface area", outcome_label = "Educational attainment"),
  Steiger_table(exposure = vol, outcome = EA, exposure_label = "Volume", outcome_label = "Educational attainment"),
  Steiger_table(exposure = IC, outcome = EA, exposure_label = "Intrinsic curvature", outcome_label = "Educational attainment"),
  Steiger_table(exposure = LOAD, outcome = ODI, exposure_label = "Late-onset Alzheimer's disease", outcome_label = "Orientation dispersion index"),
  Steiger_table(exposure = LOAD, outcome = MD_WM, exposure_label = "Late-onset Alzheimer's disease", outcome_label = "Mean diffusivity (white matter tracts)"),
  Steiger_table(exposure = LOAD, outcome = hippocampus_vol_left, exposure_label = "Late-onset Alzheimer's disease", outcome_label = "Volume of left hippocampus"),
  Steiger_table(exposure = LOAD, outcome = hippocampus_vol_right, exposure_label = "Late-onset Alzheimer's disease", outcome_label = "Volume of right hippocampus")
)

Steiger_results$P.value..all.SNPs. = format(Steiger_results$P.value..all.SNPs., digits = 3, scientific = TRUE)
Steiger_results$P.value..valid.SNPs. = format(Steiger_results$P.value..valid.SNPs., digits = 3, scientific = TRUE)

for (i in 1:length(Steiger_results$exposure)) {
  if (as.numeric(Steiger_results$P.value..all.SNPs.[i]) >= 0.001) Steiger_results$P.value..all.SNPs.[i] = sprintf("%.3f", round(as.numeric(Steiger_results$P.value..all.SNPs.[i]), 3))
  if (as.numeric(Steiger_results$P.value..valid.SNPs.[i]) >= 0.001) Steiger_results$P.value..valid.SNPs.[i] = sprintf("%.3f", round(as.numeric(Steiger_results$P.value..valid.SNPs.[i]), 3))
}

# List of SNPs supplementary tables
#------------------------------------------------------------------------------#

# Primary analysis I: EA -> LOAD
SNPs_table(exp = EA, out = LOAD)

# Primary analysis II: EA -> IDPs
SNPs_table(exp = EA, out = SA)
SNPs_table(exp = EA, out = vol)
SNPs_table(exp = EA, out = CT)
SNPs_table(exp = EA, out = LGI)
SNPs_table(exp = EA, out = MC)
SNPs_table(exp = EA, out = IC)
SNPs_table(exp = EA, out = FA_cortical)
SNPs_table(exp = EA, out = MD_cortical)
SNPs_table(exp = EA, out = ICVF)
SNPs_table(exp = EA, out = ODI)
SNPs_table(exp = EA, out = FA_WM)
SNPs_table(exp = EA, out = MD_WM)
SNPs_table(exp = EA, out = hippocampus_vol_left)
SNPs_table(exp = EA, out = hippocampus_vol_right)
SNPs_table(exp = EA, out = WMH_vol)

# Primary analysis III: IDPs -> LOAD
SNPs_table(exp = SA, out = LOAD)
SNPs_table(exp = vol, out = LOAD)
SNPs_table(exp = CT, out = LOAD)
SNPs_table(exp = LGI, out = LOAD)
SNPs_table(exp = MC, out = LOAD)
SNPs_table(exp = IC, out = LOAD)
SNPs_table(exp = FA_cortical, out = LOAD)
SNPs_table(exp = MD_cortical, out = LOAD)
SNPs_table(exp = ICVF, out = LOAD)
SNPs_table(exp = ODI, out = LOAD)
SNPs_table(exp = FA_WM, out = LOAD)
SNPs_table(exp = MD_WM, out = LOAD)
SNPs_table(exp = hippocampus_vol_left, out = LOAD)
SNPs_table(exp = hippocampus_vol_right, out = LOAD)
SNPs_table(exp = WMH_vol, out = LOAD)

# Secondary analysis I: LOAD -> EA
SNPs_table(exp = LOAD, out = EA)

# Secondary analysis II: IDPs -> EA
SNPs_table(exp = SA, out = EA)
SNPs_table(exp = vol, out = EA)
SNPs_table(exp = CT, out = EA)
SNPs_table(exp = LGI, out = EA)
SNPs_table(exp = MC, out = EA)
SNPs_table(exp = IC, out = EA)
SNPs_table(exp = FA_cortical, out = EA)
SNPs_table(exp = MD_cortical, out = EA)
SNPs_table(exp = ICVF, out = EA)
SNPs_table(exp = ODI, out = EA)
SNPs_table(exp = FA_WM, out = EA)
SNPs_table(exp = MD_WM, out = EA)
SNPs_table(exp = hippocampus_vol_left, out = EA)
SNPs_table(exp = hippocampus_vol_right, out = EA)
SNPs_table(exp = WMH_vol, out = EA)

# Secondary analysis III: LOAD -> IDPs
SNPs_table(exp = LOAD, out = SA)
SNPs_table(exp = LOAD, out = vol)
SNPs_table(exp = LOAD, out = CT)
SNPs_table(exp = LOAD, out = LGI)
SNPs_table(exp = LOAD, out = MC)
SNPs_table(exp = LOAD, out = IC)
SNPs_table(exp = LOAD, out = FA_cortical)
SNPs_table(exp = LOAD, out = MD_cortical)
SNPs_table(exp = LOAD, out = ICVF)
SNPs_table(exp = LOAD, out = ODI)
SNPs_table(exp = LOAD, out = FA_WM)
SNPs_table(exp = LOAD, out = MD_WM)
SNPs_table(exp = LOAD, out = hippocampus_vol_left)
SNPs_table(exp = LOAD, out = hippocampus_vol_right)
SNPs_table(exp = LOAD, out = WMH_vol)
