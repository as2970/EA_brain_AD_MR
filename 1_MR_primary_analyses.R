# Setup
#------------------------------------------------------------------------------#

library(tidyverse)
library(data.table)
library(ivpack)
library(meta)
library(devtools)
library(pacman)
library(MendelianRandomization)
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

#'Creates table of IVW results for the exposure-outcome combination
ivw_table = function(exposure, outcome) {
  load(paste0(deparse(substitute(exposure)), "_", deparse(substitute(outcome)), ".RData"))
  dat <- dat[dat$mr_keep == TRUE, ]
  mr_object = mr_input(bx = dat$beta.exposure, 
                       bxse = dat$se.exposure, 
                       by = dat$beta.outcome, 
                       byse = dat$se.outcome, 
                       exposure = deparse(substitute(exp)), 
                       outcome = deparse(substitute(out)), 
                       snps = dat$SNP)
  ivw_res = mr_ivw(mr_object)
  if (deparse(substitute(outcome) == "LOAD")) {
    ivw_res@Estimate = exp(ivw_res@Estimate)
    ivw_res@CILower = exp(ivw_res@CILower)
    ivw_res@CIUpper = exp(ivw_res@CIUpper)
  }
  
  if (deparse(substitute(exposure) == "LOAD")) {
    ivw_res@Estimate = log(2)*(ivw_res@Estimate)
    ivw_res@CILower = log(2)*(ivw_res@CILower)
    ivw_res@CIUpper = log(2)*(ivw_res@CIUpper)
  }
  
  if (ivw_res@Pvalue < 0.001) {
    df = data.frame(exposure = deparse(substitute(exposure)), 
                    outcome = deparse(substitute(outcome)), 
                    N_SNPs = ivw_res@SNPs,
                    Estimate = paste0(sprintf("%.2f", round(ivw_res@Estimate, 2)), " (", sprintf("%.2f", round(ivw_res@CILower, 2)), ", ", sprintf("%.2f", round(ivw_res@CIUpper, 2)), ")"),
                    P = format(ivw_res@Pvalue, digits = 3, scientific = TRUE)
    )
  } else {
    df = data.frame(exposure = deparse(substitute(exposure)), 
                    outcome = deparse(substitute(outcome)), 
                    N_SNPs = ivw_res@SNPs,
                    Estimate = paste0(sprintf("%.2f", round(ivw_res@Estimate, 2)), " (", sprintf("%.2f", round(ivw_res@CILower, 2)), ", ", sprintf("%.2f", round(ivw_res@CIUpper, 2)), ")"),
                    P = sprintf("%.3f", round(ivw_res@Pvalue, 3))
    )
  }
  return(df)
}

# Primary analysis I: EA -> LOAD
#------------------------------------------------------------------------------#
# Univariable IVW MR
EA_LOAD_ivw = ivw_table(EA, LOAD)

# Primary analysis II: EA -> IDPs
#------------------------------------------------------------------------------#
# Univariable IVW MR
EA_IDPs_ivw = rbind(
  ivw_table(EA, SA),
  ivw_table(EA, vol),
  ivw_table(EA, CT),
  ivw_table(EA, LGI),
  ivw_table(EA, MC),
  ivw_table(EA, IC),
  ivw_table(EA, FA_cortical),
  ivw_table(EA, MD_cortical),
  ivw_table(EA, ICVF),
  ivw_table(EA, ODI),
  ivw_table(EA, FA_WM),
  ivw_table(EA, MD_WM),
  ivw_table(EA, hippocampus_vol_left),
  ivw_table(EA, hippocampus_vol_right),
  ivw_table(EA, WMH_vol)
)

EA_IDPs_ivw$outcome = factor(
  x = EA_IDPs_ivw$outcome,
  levels = c(
    "SA", 
    "vol", 
    "CT", 
    "LGI", 
    "MC", 
    "IC", 
    "FA_cortical", 
    "MD_cortical", 
    "ICVF", 
    "ODI", 
    "FA_WM", 
    "MD_WM", 
    "hippocampus_vol_left",
    "hippocampus_vol_right", 
    "WMH_vol"
    ), 
  labels = c(
    "Surface area", 
    "Volume",
    "Cortical thickness", 
    "Local gyrification index", 
    "Mean curvature", 
    "Intrinsic curvature", 
    "Fractional anisotropy (cortical)", 
    "Mean diffusivity (cortical)", 
    "Intracellular volume fraction", 
    "Orientation dispersion index", 
    "Fractional anisotropy (white matter tracts)", 
    "Mean diffusivity (white matter tracts)", 
    "Volume of left hippocampus", 
    "Volume of right hippocampus", 
    "White matter hyperintensities volume"
  )
)

# Primary analysis III: IDPs -> LOAD
#------------------------------------------------------------------------------#
# Univariable IVW MR
IDPs_LOAD_ivw = rbind(
  ivw_table(SA, LOAD),
  ivw_table(vol, LOAD),
  ivw_table(CT, LOAD),
  ivw_table(LGI, LOAD),
  ivw_table(MC, LOAD),
  ivw_table(IC, LOAD),
  ivw_table(MD_cortical, LOAD),
  ivw_table(ICVF, LOAD),
  ivw_table(ODI, LOAD),
  ivw_table(FA_WM, LOAD),
  ivw_table(MD_WM, LOAD),
  ivw_table(hippocampus_vol_left, LOAD),
  ivw_table(hippocampus_vol_right, LOAD),
  ivw_table(WMH_vol, LOAD)
)

IDPs_LOAD_ivw$exposure = factor(
  x = IDPs_LOAD_ivw$exposure,
  levels = c(
    "SA", 
    "vol", 
    "CT", 
    "LGI", 
    "MC", 
    "IC", 
    "FA_cortical", 
    "MD_cortical", 
    "ICVF", 
    "ODI", 
    "FA_WM", 
    "MD_WM", 
    "hippocampus_vol_left",
    "hippocampus_vol_right", 
    "WMH_vol"
  ), 
  labels = c(
    "Surface area", 
    "Volume",
    "Cortical thickness", 
    "Local gyrification index", 
    "Mean curvature", 
    "Intrinsic curvature", 
    "Fractional anisotropy (cortical)", 
    "Mean diffusivity (cortical)", 
    "Intracellular volume fraction", 
    "Orientation dispersion index", 
    "Fractional anisotropy (white matter tracts)", 
    "Mean diffusivity (white matter tracts)", 
    "Volume of left hippocampus", 
    "Volume of right hippocampus", 
    "White matter hyperintensities volume"
  )
)

# Secondary analysis I: LOAD -> EA
#------------------------------------------------------------------------------#
# Univariable IVW MR
LOAD_EA_ivw = ivw_table(LOAD, EA)

# Secondary analysis II: IDPs -> EA
#------------------------------------------------------------------------------#
# Univariable IVW MR
IDPs_EA_ivw = rbind(
  ivw_table(SA, EA),
  ivw_table(vol, EA),
  ivw_table(CT, EA),
  ivw_table(LGI, EA),
  ivw_table(MC, EA),
  ivw_table(IC, EA),
  ivw_table(MD_cortical, EA),
  ivw_table(ICVF, EA),
  ivw_table(ODI, EA),
  ivw_table(FA_WM, EA),
  ivw_table(MD_WM, EA),
  ivw_table(hippocampus_vol_left, EA),
  ivw_table(hippocampus_vol_right, EA),
  ivw_table(WMH_vol, EA)
)

IDPs_EA_ivw$exposure = factor(
  x = IDPs_EA_ivw$exposure,
  levels = c(
    "SA", 
    "vol", 
    "CT", 
    "LGI", 
    "MC", 
    "IC", 
    "FA_cortical", 
    "MD_cortical", 
    "ICVF", 
    "ODI", 
    "FA_WM", 
    "MD_WM", 
    "hippocampus_vol_left",
    "hippocampus_vol_right", 
    "WMH_vol"
  ), 
  labels = c(
    "Surface area", 
    "Volume",
    "Cortical thickness", 
    "Local gyrification index", 
    "Mean curvature", 
    "Intrinsic curvature", 
    "Fractional anisotropy (cortical)", 
    "Mean diffusivity (cortical)", 
    "Intracellular volume fraction", 
    "Orientation dispersion index", 
    "Fractional anisotropy (white matter tracts)", 
    "Mean diffusivity (white matter tracts)", 
    "Volume of left hippocampus", 
    "Volume of right hippocampus", 
    "White matter hyperintensities volume"
  )
)

# Secondary analysis III: LOAD -> IDPs
#------------------------------------------------------------------------------#
# Univariable IVW MR
LOAD_IDPs_ivw = rbind(
  ivw_table(LOAD, SA),
  ivw_table(LOAD, vol),
  ivw_table(LOAD, CT),
  ivw_table(LOAD, LGI),
  ivw_table(LOAD, MC),
  ivw_table(LOAD, IC),
  ivw_table(LOAD, FA_cortical),
  ivw_table(LOAD, MD_cortical),
  ivw_table(LOAD, ICVF),
  ivw_table(LOAD, ODI),
  ivw_table(LOAD, FA_WM),
  ivw_table(LOAD, MD_WM),
  ivw_table(LOAD, hippocampus_vol_left),
  ivw_table(LOAD, hippocampus_vol_right),
  ivw_table(LOAD, WMH_vol)
)

LOAD_IDPs_ivw$outcome = factor(
  x = LOAD_IDPs_ivw$outcome,
  levels = c(
    "SA", 
    "vol", 
    "CT", 
    "LGI", 
    "MC", 
    "IC", 
    "FA_cortical", 
    "MD_cortical", 
    "ICVF", 
    "ODI", 
    "FA_WM", 
    "MD_WM", 
    "hippocampus_vol_left",
    "hippocampus_vol_right", 
    "WMH_vol"
  ), 
  labels = c(
    "Surface area", 
    "Volume",
    "Cortical thickness", 
    "Local gyrification index", 
    "Mean curvature", 
    "Intrinsic curvature", 
    "Fractional anisotropy (cortical)", 
    "Mean diffusivity (cortical)", 
    "Intracellular volume fraction", 
    "Orientation dispersion index", 
    "Fractional anisotropy (white matter tracts)", 
    "Mean diffusivity (white matter tracts)", 
    "Volume of left hippocampus", 
    "Volume of right hippocampus", 
    "White matter hyperintensities volume"
  )
)
