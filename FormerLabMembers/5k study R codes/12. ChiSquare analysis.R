# this script generates contingency table for cell types

remove(list = ls())

library(readr, quietly = TRUE)
library(readxl, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(data.table, quietly = TRUE)
library(foreach, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(stringr)
library(lme4)
library(ggpubr)
library(openxlsx)
library(glmnet)
library(tidyr)
library(broom)
library(REdaS)
library(vioplot)
library(epiR)
library(gridExtra)
library(ggpmisc)

# load data table

results <- read_csv("./input/MASTERSPREADSHEET_5K_analysis.csv")


# filter out only markers with assigned cell-specificity

results_clean <- results %>% 
  filter(!is.na(`Single cell RNA expression top cell type`))

outcome <- c("cel","t2ll","combiwise","brain_damage","GMSdis","drb1_15_0103","bmi","bmi_cat2","smoking","ethnicity")

# generate vector of cell types

cell_types <- unique(results_clean$`Single cell RNA expression top cell type`)


cel_type_res <- as.data.frame(matrix(nrow = 0,ncol = 4))
names(cel_type_res) <- c("cell_type","count","directionality","outcome")



for(i in 1:length(outcome)){
  
  # i <- 2
  
  
  dat <- results_clean %>% 
    filter(.data[[paste0(outcome[i],"_coef_pval")]]<=0.05 & .data[[paste0(outcome[i],"_pval_FDR")]]<=0.05)
  
  dat_pos <- dat %>% 
    filter(.data[[paste0(outcome[i],"_coef")]]>0)
  
  dat_neg <- dat %>% 
    filter(.data[[paste0(outcome[i],"_coef")]]<0)
  
  
  cell_pos <- tibble::enframe(table(dat_pos$`Single cell RNA expression top cell type`),name = "cell_type",value = "count")
  cell_pos$directionality <- "positive"
  
  cell_neg <- tibble::enframe(table(dat_neg$`Single cell RNA expression top cell type`),name = "cell_type",value = "count")
  cell_neg$directionality <- "negative"
  
  # merge positive and negative
  
  cell_all <- rbind(cell_pos,cell_neg)
  
  cell_all$outcome <- outcome[i]
  
  # add to all results 
  
  cel_type_res <- rbind(cel_type_res,cell_all)
  
  }
  
###########################
  # add CEL_nfl residuals 


dat <- results_clean %>% 
  filter(cel_coef_pval<=0.05 & cel_pval_FDR<=0.05) %>% 
  filter(NFL_CEL_residuals_coef_pval<=0.05 & NFL_CEL_residuals_pval_FDR<=0.05)

dat_pos <- dat %>% 
  filter(NFL_CEL_residuals_coef>0)

dat_neg <- dat %>% 
  filter(NFL_CEL_residuals_coef<0)


cell_pos <- tibble::enframe(table(dat_pos$`Single cell RNA expression top cell type`),name = "cell_type",value = "count")
cell_pos$directionality <- "positive"

cell_neg <- tibble::enframe(table(dat_neg$`Single cell RNA expression top cell type`),name = "cell_type",value = "count")
cell_neg$directionality <- "negative"

# merge positive and negative

cell_all <- rbind(cell_pos,cell_neg)

cell_all$outcome <- "NFL_CEL_residuals"

# add to all results 

cel_type_res <- rbind(cel_type_res,cell_all)

###########################
# add sc_disability (must significanlty correlate with COmbiWISE AND must significantly correlate with SC disability)


dat <- results_clean %>% 
  filter(combiwise_coef_pval<=0.05 & combiwise_pval_FDR<=0.05) %>% 
  filter(sc_disability_coef_pval<=0.05 & sc_disability_pval_FDR<=0.05)

dat_pos <- dat %>% 
  filter(sc_disability_coef>0) %>% 
  filter(combiwise_coef>0)

dat_neg <- dat %>% 
  filter(sc_disability_coef<0) %>% 
  filter(combiwise_coef<0)


cell_pos <- tibble::enframe(table(dat_pos$`Single cell RNA expression top cell type`),name = "cell_type",value = "count")
cell_pos$directionality <- "positive"

cell_neg <- tibble::enframe(table(dat_neg$`Single cell RNA expression top cell type`),name = "cell_type",value = "count")
cell_neg$directionality <- "negative"

# merge positive and negative

cell_all <- rbind(cell_pos,cell_neg)

cell_all$outcome <- "sc_disability"

# add to all results 

cel_type_res <- rbind(cel_type_res,cell_all)



###########################
# add MS_vs_HD residuals 


dat <- results_clean %>% 
  filter(MS_vs_HD_diff_Wilcox_FDR <=0.05) 

dat_pos <- dat %>% 
  filter(MS_HD_median_diff>0)

dat_neg <- dat %>% 
  filter(MS_HD_median_diff<0)


cell_pos <- tibble::enframe(table(dat_pos$`Single cell RNA expression top cell type`),name = "cell_type",value = "count")
cell_pos$directionality <- "positive"

cell_neg <- tibble::enframe(table(dat_neg$`Single cell RNA expression top cell type`),name = "cell_type",value = "count")
cell_neg$directionality <- "negative"

# merge positive and negative

cell_all <- rbind(cell_pos,cell_neg)

cell_all$outcome <- "MS_vs_HD"

# add to all results 

cel_type_res <- rbind(cel_type_res,cell_all)

# save the file 
write_csv(cel_type_res,"./output/chiSquare_stat/all_outcomes_cell_specificity.csv")


# spread the count into positive and negative

cel_type_res_wide <- cel_type_res %>% 
  pivot_wider(names_from = directionality,values_from = count)

# replace NAs with 0

cel_type_res_wide$positive <- ifelse(is.na(cel_type_res_wide$positive),yes = 0,no = cel_type_res_wide$positive)
cel_type_res_wide$negative <- ifelse(is.na(cel_type_res_wide$negative),yes = 0,no = cel_type_res_wide$negative)

# run chisquare on individual lines

chi_sq_results <- cel_type_res_wide
chi_sq_results$chi_sqr_pval <- NA

for(i in 1:length(cel_type_res_wide$cell_type)){
  
  mat <- cel_type_res_wide[i,]
  mat <- mat %>% 
    dplyr::select(-c(outcome,cell_type))
  chisq <- chisq.test(mat)
  chi_sq_results$chi_sqr_pval[i] <- chisq$p.value
  
}

chi_sq_results$positive <- as.vector(chi_sq_results$positive )
chi_sq_results$negative <- as.vector(chi_sq_results$negative )

# FDR adjust
chi_sq_results$chi_sqr_pval_fdr <- p.adjust(chi_sq_results$chi_sqr_pval,method = "fdr")


# save the file 
write_csv(chi_sq_results,"./output/chiSquare_stat/all_outcomes_cell_specificity_chi_square_results.csv")


# convert it into wide table 

final_results <- tibble::enframe(cell_types,name = NULL,value = "cell_type")

outcomes <- unique(chi_sq_results$outcome)

for(i in 1:length(outcomes)){
  # i <- 1
  
  dat <- chi_sq_results %>% 
    filter(outcome==outcomes[i]) %>% 
    dplyr::select(-outcome)
  
  
  names(dat)[-1] <- paste(outcomes[i],names(dat)[-1],sep = "_")
  
  # merge with final results
  
  final_results <- left_join(final_results,dat)

  
}

# save the file 
write_csv(final_results,"./output/chiSquare_stat/cell_type_chi_square_results.csv")

