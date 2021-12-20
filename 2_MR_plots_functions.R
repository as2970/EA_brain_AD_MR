# Setup
#------------------------------------------------------------------------------#

library(tidyverse)
library(pacman)
library(ggplot2)
library(plyr)
library(forestplot)
library(anchors)
library(grid)
library(magrittr)
library(checkmate)
library(extrafont)
loadfonts()

# Functions
#------------------------------------------------------------------------------#

#' MR forest plot
MR_ForestPlot <- function(exposure, outcome, exposure_label, outcome_label, ytextsize, exponentiate = FALSE) {
  
  p_unload(MendelianRandomization)
  library(TwoSampleMR)
  
  load(paste0(deparse(substitute(exposure)), "_", deparse(substitute(outcome)), ".RData"))
  dat = dat[dat$mr_keep == TRUE, ]
  
  
  res_sin = mr_singlesnp(dat, all_method = c("mr_ivw_mre"))
  res_sin$exposure = exposure_label
  res_sin$outcome = outcome_label
  
  res_sin$SNP[res_sin$SNP == "All - Inverse variance weighted (multiplicative random effects)"] = "MR-IVW estimate"
  res_sin$UCL = res_sin$b + qnorm(0.975) * res_sin$se
  res_sin$LCL = res_sin$b - qnorm(0.975) * res_sin$se
  res_sin = res_sin[ , c("exposure", "outcome", "SNP", "b", "UCL", "LCL")]
  
  SNPs = res_sin$SNP[!(res_sin$SNP == "MR-IVW estimate")]
  SNPs_ordered = SNPs[order(res_sin$b)]
  
  res_sin = rbind(res_sin, res_sin[nrow(res_sin), ])
  res_sin$SNP[nrow(res_sin)-1] = ""
  res_sin$b[nrow(res_sin)-1] = NA
  res_sin$UCL[nrow(res_sin)-1] = NA
  res_sin$LCL[nrow(res_sin)-1] = NA
  
  res_sin$SNP = ordered(res_sin$SNP, levels = c("MR-IVW estimate", "", SNPs_ordered))
  
  res_sin$dummy = 0
  res_sin$dummy[res_sin$SNP == "MR-IVW estimate"] = 1
  
  null_value = 0
  
  
  if (deparse(substitute(exposure) == "LOAD")) {
    res_sin$b = log(2)*res_sin$b
    res_sin$UCL = log(2)*res_sin$UCL
    res_sin$LCL = log(2)*res_sin$LCL
  }
  
    plot = ggplot(res_sin, aes(y = SNP, x = b)) +
      geom_vline(xintercept = null_value, linetype = "dotted") +
      geom_errorbarh(aes(xmin = LCL, xmax = UCL, size = as.factor(dummy), colour = as.factor(dummy)), height=0) +
      geom_point(aes(colour = as.factor(dummy))) +
      geom_hline(aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
      scale_colour_manual(values=c("black", "darksalmon")) +
      scale_size_manual(values=c(0.3, 1)) +
      labs(y = "", x = paste0("\nMR effect size for\n","'", res_sin$exposure[1], "' on\n'", res_sin$outcome[1], "'")) +
      theme_bw() +
      theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +    
      theme(legend.position = "none") + 
      theme(axis.text.x=element_text(size = 10, family = "Latin Modern Roman 10 Regular")) + 
      theme(axis.title.x = element_text(size = 12, family = "Latin Modern Roman 10 Regular")) +
      theme(axis.text.y=element_text(size = ytextsize, family = "Latin Modern Roman 10 Regular")) 
  
  return(plot)
}

#' MR leave-one-out-plot
MR_LeaveOneOutPlot <- function(exp, out, exp_label, out_label, ytextsize, exponentiate = FALSE) {
  
  p_unload(MendelianRandomization)
  library(TwoSampleMR)

  load(paste0(deparse(substitute(exp)), "_", deparse(substitute(out)), ".RData"))
  dat = dat[dat$mr_keep == TRUE, ]
  
  res_loo = mr_leaveoneout(dat, method = mr_ivw_mre)
  res_loo$exposure = exp_label
  res_loo$outcome = out_label
  
  res_loo$UCL = res_loo$b + qnorm(0.975) * res_loo$se
  res_loo$LCL = res_loo$b - qnorm(0.975) * res_loo$se
  res_loo = res_loo[ , c("exposure", "outcome", "SNP", "b", "UCL", "LCL")]
  
  SNPs = res_loo$SNP[!(res_loo$SNP == "All")]
  SNPs_ordered = SNPs[order(res_loo$b)]
  
  res_loo = rbind(res_loo, res_loo[nrow(res_loo), ])
  res_loo$SNP[nrow(res_loo)-1] = ""
  res_loo$b[nrow(res_loo)-1] = NA
  res_loo$UCL[nrow(res_loo)-1] = NA
  res_loo$LCL[nrow(res_loo)-1] = NA
  
  res_loo$SNP = ordered(res_loo$SNP, levels = c("All", "", SNPs_ordered))
  
  res_loo$dummy = 0
  res_loo$dummy[res_loo$SNP == "All"] = 1
  
  null_value = 0
  
  if(exponentiate == TRUE) {
    res_loo$b = base::exp(res_loo$b)
    res_loo$UCL = base::exp(res_loo$UCL)
    res_loo$LCL = base::exp(res_loo$LCL)
    null_value = 1
  }
  
  if (deparse(substitute(exp) == "LOAD")) {
    res_loo$b = log(2)*res_loo$b
    res_loo$UCL = log(2)*res_loo$UCL
    res_loo$LCL = log(2)*res_loo$LCL
  }
  
  plot = ggplot(res_loo, aes(y = SNP, x = b)) +
    geom_vline(xintercept = null_value, linetype = "dotted") +
    geom_errorbarh(aes(xmin = LCL, xmax = UCL, size = as.factor(dummy), colour = as.factor(dummy)), height=0) +
    geom_point(aes(colour = as.factor(dummy))) +
    geom_hline(aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
    scale_colour_manual(values=c("black", "darkseagreen3")) +
    scale_size_manual(values=c(0.3, 1)) +
    labs(y = "", x = paste0("\nMR effect size for\n","'", res_loo$exposure[1], "' on\n'", res_loo$outcome[1], "'")) +
    theme_bw() + 
    theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +  
    theme(legend.position = "none") + 
    theme(axis.text.x=element_text(size = 15, family = "Latin Modern Roman 10 Regular")) + 
    theme(axis.title.x = element_text(size = 20, family = "Latin Modern Roman 10 Regular")) +
    theme(axis.text.y=element_text(size = ytextsize, family = "Latin Modern Roman 10 Regular"))
  
  return(plot)
}


#' MR Funnel plot
MR_FunnelPlot <- function(exp, out, exp_label, out_label)
{
  
  p_unload(MendelianRandomization)
  library(TwoSampleMR)
  
  load(paste0(deparse(substitute(exp)), "_", deparse(substitute(out)), ".RData"))
  dat = dat[dat$mr_keep == TRUE, ]
  
  res_sin = mr_singlesnp(dat, all_method = c("mr_ivw_mre"))
  res_sin$exposure = exp_label
  res_sin$outcome = out_label
  
  res_sin$SNP[res_sin$SNP == "All - Inverse variance weighted (multiplicative random effects)"] = "MR-IVW estimate"
  
  plot = ggplot(res_sin[res_sin$SNP != "MR-IVW estimate", ], aes(y = 1/se, x = b)) +
    geom_point() +
    geom_vline(xintercept = res_sin[res_sin$SNP == "MR-IVW estimate", "b"], linetype = "solid", color = "slategray3", size = 1) + 
    labs(y = expression(1/SE[IV]), x = expression(beta[IV])) +
    theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +  
    theme_bw() +
    theme(axis.text = element_text(size = 10, family = "Latin Modern Roman 10 Regular"), 
          axis.title = element_text(size = 13, family = "Latin Modern Roman 10 Regular"))
  
  return(plot)
}

#' MR Scatter plot
MR_ScatterPlot <- function(exp, out, exp_label, out_label) {
  
  p_unload(TwoSampleMR)
  library(MendelianRandomization)
  
  load(paste0(deparse(substitute(exp)), "_", deparse(substitute(out)), ".RData"))
  dat = dat[dat$mr_keep == TRUE, ]
  
  dat$exposure = exp_label
  dat$outcome = out_label
  
  for (i in 1:length(dat$SNP)) {
    if(dat$beta.exposure[i] < 0) {
      dat$beta.exposure[i] = (-1) * dat$beta.exposure[i]
      dat$beta.outcome[i]= (-1) * dat$beta.outcome[i]
    }
  }
  
  mr_object = mr_input(bx = dat$beta.exposure, 
                       bxse = dat$se.exposure, 
                       by = dat$beta.outcome, 
                       byse = dat$se.outcome, 
                       exposure = deparse(substitute(exp)), 
                       outcome = deparse(substitute(out)), 
                       snps = dat$SNP)
  
  est_ivw = mr_ivw(mr_object)
  
  ggplot(data = dat, aes(x = beta.exposure, y = beta.outcome)) +
    geom_errorbar(aes(ymin = beta.outcome - qnorm(0.975) * se.outcome, ymax = beta.outcome + qnorm(0.975) * se.outcome), colour="grey", width=0) +
    geom_errorbarh(aes(xmin = beta.exposure - qnorm(0.975) * se.exposure, xmax = beta.exposure + qnorm(0.975) * se.exposure), colour="grey", height=0) +
    geom_point(aes(text = paste("SNP:", SNP))) +
    geom_abline(intercept = 0, slope = est_ivw$Estimate, linetype = "solid", color = "goldenrod3", size = 1) +
    labs(x = paste("\nSNP effect on", dat$exposure[1]), y=paste("SNP effect on", dat$outcome[1], "\n")) +
    guides(colour = guide_legend(ncol=2)) + 
    theme_bw() + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
    geom_vline(xintercept=0, linetype = "dashed", color = "black") +
    theme(panel.background = element_blank(), panel.grid.minor = element_blank()) +  
    theme(axis.title.y = element_text(size = 15, angle = 90, family = "Latin Modern Roman 10 Regular")) +
    theme(axis.title.x = element_text(size = 15, family = "Latin Modern Roman 10 Regular"))+ 
    guides(col = guide_legend(nrow = 5)) +
    theme(axis.text = element_text(size = 10, family = "Latin Modern Roman 10 Regular"))
}

#' MR Estimates plot

MR_EstimatesPlot <- function(exp, out, ticks, clips) {
  
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
  
  
  est_ivw = mr_ivw(mr_object)
  est_med = mr_median(mr_object)
  est_conmix = mr_conmix(mr_object)
  est_egger = mr_egger(mr_object)
  presso_df = data.frame("betaY" = mr_object@betaY, "betaX" = mr_object@betaX, "betaYse" = mr_object@betaYse, "betaXse" = mr_object@betaXse)
  est_presso = mr_presso(BetaOutcome = "betaY", BetaExposure = "betaX", SdOutcome = "betaYse", SdExposure = "betaXse",
                         data = presso_df, OUTLIERtest = TRUE, DISTORTIONtest = TRUE, NbDistribution = 12000,
                         SignifThreshold = 0.05)

  
  if (is.na(est_presso$'Main MR results'[2, 3])){
    MR_estimates = data.frame("Method" = factor(c("MR-IVW", "MR-Median", "MR-ConMix", "MR-PRESSO", "MR-Egger"), levels = c("MR-Egger", "MR-PRESSO", "MR-ConMix", "MR-Median", "MR-IVW")),
                              "SNPs" = length(mr_object$snps),
                              "Estimate" = c(est_ivw$Estimate, est_med$Estimate, est_conmix$Estimate, est_presso$'Main MR results'[1, 3], est_egger$Estimate),
                              "CILower" = c(est_ivw$CILower, est_med$CILower, est_conmix$CILower, est_presso$'Main MR results'[1, 3] - qt(0.975, dim(presso_df)[1]-2) * est_presso$'Main MR results'[1, 4], est_egger$CILower.Est),
                              "CIUpper" = c(est_ivw$CIUpper, est_med$CIUpper, est_conmix$CIUpper, est_presso$'Main MR results'[1, 3] + qt(0.975, dim(presso_df)[1]-2) * est_presso$'Main MR results'[1, 4], est_egger$CIUpper.Est), 
                              "Pvalue" = c(est_ivw$Pvalue, est_med$Pvalue, est_conmix$Pvalue, est_presso$'Main MR results'[1, 6], est_egger$Pvalue.Est)
    )
  } else {
    MR_estimates = data.frame("Method" = factor(c("MR-IVW", "MR-Median", "MR-ConMix", "MR-PRESSO", "MR-Egger"), levels = c("MR-Egger", "MR-PRESSO", "MR-ConMix", "MR-Median", "MR-IVW")),
                              "SNPs" = length(mr_object$snps),
                              "Estimate" = c(est_ivw$Estimate, est_med$Estimate, est_conmix$Estimate, est_presso$'Main MR results'[2, 3], est_egger$Estimate),
                              "CILower" = c(est_ivw$CILower, est_med$CILower, est_conmix$CILower, est_presso$'Main MR results'[2, 3] - qt(0.975, dim(presso_df)[1]-2) * est_presso$'Main MR results'[2, 4], est_egger$CILower.Est),
                              "CIUpper" = c(est_ivw$CIUpper, est_med$CIUpper, est_conmix$CIUpper, est_presso$'Main MR results'[2, 3] + qt(0.975, dim(presso_df)[1]-2) * est_presso$'Main MR results'[2, 4], est_egger$CIUpper.Est),
                              "Pvalue" = c(est_ivw$Pvalue, est_med$Pvalue, est_conmix$Pvalue, est_presso$'Main MR results'[2, 6], est_egger$Pvalue.Est)
    )
  }
  
  if (deparse(substitute(exp) == "LOAD")) {
    MR_estimates$Estimate = log(2) * MR_estimates$Estimate
    MR_estimates$CILower = log(2) * MR_estimates$CILower
    MR_estimates$CIUpper = log(2) * MR_estimates$CIUpper
  }
  
  MR_estimates$Estimate = sprintf("%.2f",round((MR_estimates$Estimate),2))
  MR_estimates$CILower = sprintf("%.2f",round((MR_estimates$CILower),2))
  MR_estimates$CIUpper = sprintf("%.2f",round((MR_estimates$CIUpper),2))
  MR_estimates$Overall = paste0(MR_estimates$Estimate, " (", MR_estimates$CILower, ", ", MR_estimates$CIUpper, ")")
  MR_estimates$Pvalue = format(MR_estimates$Pvalue, digits = 3, scientific = TRUE)
  
  for (i in 1:5) {
    if (as.numeric(MR_estimates$Pvalue[i]) >= 0.001) MR_estimates$Pvalue[i] = sprintf("%.3f", round(as.numeric(MR_estimates$Pvalue[i]), 3))
  }

  headers_row = c(NA, NA, NA, NA, NA, NA, NA, NA)
  MR_estimates <- rbind(headers_row, MR_estimates)
  
  tabletext = list(
    list("MR Method", "MR-IVW", "MR-Median", "MR-ConMix", "MR-PRESSO", "MR-Egger"),
    append(list("SNPs"), MR_estimates[!is.na(MR_estimates$SNPs), "SNPs"]),
    append(list("Estimate (95% CI)"), MR_estimates[!is.na(MR_estimates$Overall), "Overall"]), 
    append(list("p-value"), list(
      eval(bquote(expression(.(str_sub(MR_estimates$Pvalue[2], 1, 4))%*%10^.(str_sub(MR_estimates$Pvalue[2], -3, -1))))), 
      eval(bquote(expression(.(str_sub(MR_estimates$Pvalue[3], 1, 4))%*%10^.(str_sub(MR_estimates$Pvalue[3], -3, -1))))), 
      eval(bquote(expression(.(str_sub(MR_estimates$Pvalue[4], 1, 4))%*%10^.(str_sub(MR_estimates$Pvalue[4], -3, -1))))), 
      eval(bquote(expression(.(str_sub(MR_estimates$Pvalue[5], 1, 4))%*%10^.(str_sub(MR_estimates$Pvalue[5], -3, -1))))), 
      eval(bquote(expression(.(str_sub(MR_estimates$Pvalue[6], 1, 4))%*%10^.(str_sub(MR_estimates$Pvalue[6], -3, -1))))))
      )
    )
  
  for (i in 2:6) {
    if (as.numeric(MR_estimates$Pvalue[i]) >= 0.001) tabletext[[4]][[i]] = sprintf("%.3f", round(as.numeric(MR_estimates$Pvalue[i]), 3))
    else {
      if (str_sub(MR_estimates$Pvalue[i], -2, -2) == "0") tabletext[[4]][[i]] = eval(bquote(expression( .(str_sub(MR_estimates$Pvalue[i], 1, 4))%*%10^-.(str_sub(MR_estimates$Pvalue[i], -1, -1)) )))
    }
  }
    grid.newpage()
    forestplot(labeltext = tabletext,
               mean = as.numeric(MR_estimates$Estimate),
               lower = as.numeric(MR_estimates$CILower),
               upper = as.numeric(MR_estimates$CIUpper),
               align = c("l","c","c","c"),
               is.summary = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
               graph.pos = "right",
               clip = clips,
               xlab = expression(paste(beta, " (95% CI)")), 
               zero = 0,
               graphwidth = unit(10, "cm"),
               colgap = unit(5, "mm"),
               lineheight = unit(1, "cm"),
               line.margin = unit(0, "mm"),
               col = fpColors(box = "black", line= "black"),
               txt_gp = fpTxtGp(label = gpar(fontfamily = "Latin Modern Roman 10"), xlab = gpar(cex = 1), ticks = gpar(cex = 0.9), cex = 1),
               xlog = FALSE, 
               xticks = ticks,
               ci.vertices = TRUE,
               boxsize = 0.25,
    )
}
