# this script calculates slopes for longitudinal analysis
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
library(tibble)

#####################################################################
#
# load functions
#
#####################################################################

remove_outlier <-
  function(x, mult, ...) {
    upper <-
      quantile(x, 3 / 4, na.rm = TRUE) + mult * IQR(x, na.rm = TRUE)   
    downer <-
      quantile(x, 1 / 4, na.rm = TRUE) - mult * IQR(x, na.rm = TRUE)   
    out <- ifelse(x > upper | x < downer, NA, x)
    return(out)
  }

remout_fun <- function(dat, outcomes, mult, ...) {
  dat_noout <- dat %>%
    mutate_at(outcomes, remove_outlier, mult = mult)
  return(dat_noout)
}

#####################################################################
#
# load files
#
#####################################################################

neuro <- read_csv("./input/clinic_scales.csv")
names(neuro) <- tolower(names(neuro))


mri <- read_csv("~/Desktop/MRI_export.csv")
names(mri) <- tolower(names(mri))

demo <- read_csv("./input/outcomes_5k.csv")

# isolate only MS and HD samples of interest from neuro file

samples_of_interest <- subset(demo,subset = diagnosis %in% c("HD","RR-MS","SP-MS","PP-MS"))
samples_of_interest <- unique(samples_of_interest$patientcode)


neuro <- neuro %>% 
  filter(patientcode %in% samples_of_interest) %>% 
  distinct()

# load somamers (adjusted and scaled and curbed)

soma <- fread("./input/soma_highSNR_adj_scaled_3_iqr_curbed.csv")

# load pathways 
path <- fread("./input/all_pathways.csv")

# merge soma and pathway

soma <- left_join(soma,path)


# load master spreadsheet

master <- read_csv("./input/MASTERSPREADSHEET_5K_analysis.csv")

# load translation file
trans <- read_excel("./input/translation file 5k new.xlsx")

trans <- trans[,-1]
names(trans)[1] <- "marker"


# add pathway names into translation file 
pathways <- tibble::enframe(names(path)[-1], name = NULL,value = "marker")
pathways$Gene <- pathways$marker
pathways$Target <- NA
pathways$Target_description <- NA

trans <- rbind(trans,pathways)


#######################################
# generate GMSdis disability outcome 
#######################################

# redo GMSdis outcomes in the whole cohort 

neuro_scaled <- neuro
neuro_scaled[,c("bd_disability","combiwise")] <- apply(neuro_scaled[,c("bd_disability","combiwise")] ,MARGIN = 2,FUN = scale)

# dis outcome
neuro_scaled <- neuro_scaled %>% 
  mutate(GMSdis_loc=combiwise+bd_disability) 


# build a model of GMSDis_loc from brain_damage

mod_gms_bd <- lm(GMSdis_loc~bd_disability,neuro_scaled)

summary(mod_gms_bd)

# build a model of GMSDis_loc from combiwise

mod_gms_combi <- lm(GMSdis_loc~combiwise,neuro_scaled)

summary(mod_gms_combi)

# predict GMSDis_loc from BD
neuro_scaled$GMSdis_loc[which(is.na(neuro_scaled$GMSdis_loc) & !is.na(neuro_scaled$bd_disability))] <- predict(mod_gms_bd,
                                                                                                                    newdata = neuro_scaled[which(is.na(neuro_scaled$GMSdis_loc) & !is.na(neuro_scaled$bd_disability)),]) 


# predict GMSDis_loc from CombiWISE
neuro_scaled$GMSdis_loc[which(is.na(neuro_scaled$GMSdis_loc) & !is.na(neuro_scaled$combiwise))] <- predict(mod_gms_combi,
                                                                                                                 newdata = neuro_scaled[which(is.na(neuro_scaled$GMSdis_loc) & !is.na(neuro_scaled$combiwise)),]) 
imputted_gmsdis <- neuro_scaled %>% 
  dplyr::select(id,GMSdis_loc) %>% 
  rename(GMSDIS_impute=GMSdis_loc)

# merge with neuro file 
neuro$id <- paste0(neuro$patientcode,neuro$date)

neuro <- left_join(neuro,imputted_gmsdis)

hist(neuro$GMSDIS_impute)

# check SC-disability for matching samples

sc_dis_orig <- soma %>% 
  dplyr::select(patientcode,dateneuroexam,sc_disability) %>% 
  mutate(id=paste0(patientcode,dateneuroexam)) %>% 
  dplyr::select(-c(patientcode,dateneuroexam)) %>% 
  rename(sc_dis_orig=sc_disability)

sc_dis_orig <- left_join(sc_dis_orig,neuro)

plot(sc_dis_orig$sc_disability~sc_dis_orig$sc_dis_orig)
plot(sc_disability~combiwise,sc_dis_orig)

#********************************************************#
##########################################################
#
# ANALYSIS ONE
#
##########################################################
#********************************************************#

# isolate first untreated LP in MS cohort

demo_ms <- demo %>% 
  filter(diagnosis %in% c("RR-MS","SP-MS","PP-MS")) %>% 
  filter(nice_therapy=="Untreated")


ms_first <- demo_ms %>% 
  group_by(patientcode) %>% 
  arrange(age) %>% 
  slice_head()

# get a vector of patients to analyze
ptx <- ms_first$patientcode

clin <- neuro[0,]
t2 <- mri[0,]

slopes <- tibble::enframe(ptx,name = NULL,value = "patientcode")

slopes$combi_slope <- NA
slopes$combi_r2 <- NA
slopes$combi_fu_length <- NA
slopes$combi_nr_visits <- NA

slopes$bd_slope <- NA
slopes$bd_r2 <- NA
slopes$bd_fu_length <- NA
slopes$bd_nr_visits <- NA

slopes$sc_slope <- NA
slopes$sc_r2 <- NA
slopes$sc_fu_length <- NA
slopes$sc_nr_visits <- NA

slopes$gmsdis_slope <- NA
slopes$gmsdis_r2 <- NA
slopes$gmsdis_fu_length <- NA
slopes$gmsdis_nr_visits <- NA

slopes$t2ll_slope <- NA
slopes$t2ll_r2 <- NA
slopes$t2ll_fu_length <- NA
slopes$t2ll_nr_visits <- NA

neuro_outcomes <- c("combiwise","bd_disability","GMSDIS_impute","sc_disability")
slp <- c("combi_slope","bd_slope","gmsdis_slope","sc_slope")
rsqr <- c("combi_r2","bd_r2","gmsdis_r2","sc_r2")
fu <- c("combi_fu_length","bd_fu_length","gmsdis_fu_length","sc_fu_length")
nr_visits <- c("combi_nr_visits","bd_nr_visits","gmsdis_nr_visits","sc_nr_visits")

for(i in 1:length(ptx)){
  # i <- 1
  
  for(j in 1:length(neuro_outcomes)){
    # j <- 1
    
    # isolate samples of interest - on or after first LP
    dat <- neuro %>% 
      filter(patientcode==ptx[i])%>% 
      filter(!is.na(.data[[neuro_outcomes[j]]])) %>% 
      filter(date >= ms_first$dateneuroexam[i] )
    
    
    # generate slopes
    if(dim(dat)[1]>1){   # analyze only patients with more than 1 visit
      if(max(dat$age)-min(dat$age)>=1.5){  #analyze only patients with at least 1.5y f/u
        mod <- summary(lm(dat[[neuro_outcomes[j]]]~dat$age))
        slopes[[slp[j]]][i] <- mod$coefficients[2,1]
        slopes[[rsqr[j]]][i] <- mod$r.squared
        slopes[[fu[j]]][i] <- max(dat$age)-min(dat$age)
        slopes[[nr_visits[j]]][i] <- dim(dat)[1]
      }
    }
    
    else
    {
      slopes[[slp[j]]][i] <- NA
      slopes[[rsqr[j]]][i] <- NA
      slopes[[fu[j]]][i] <-NA
      slopes[[nr_visits[j]]][i] <- NA
    }
  }
  
  # add T2LL
  dat <- mri %>% 
    filter(patientcode==ptx[i]) %>% 
    filter(!is.na(volume_lesions)) %>% 
    filter(date >= ms_first$mri_date[i]) 
    
    if(dim(dat)[1]>1){
      if(max(dat$age)-min(dat$age)>=1.5){
        mod <- summary(lm(volume_lesions~age,dat))
        slopes$t2ll_slope[i] <- mod$coefficients[2,1]
        slopes$t2ll_r2[i] <- mod$r.squared
        slopes$t2ll_fu_length[i] <- max(dat$age)-min(dat$age)
        slopes$t2ll_nr_visits[i] <- dim(dat)[1]
      }
    }
  
  else
  {
    slopes$t2ll_slope[i] <- NA
    slopes$t2ll_r2[i] <- NA
    slopes$t2ll_fu_length[i] <-NA
    slopes$t2ll_nr_visits[i] <- NA
  }
  
}


## ELIMINATE +/- 3*IQR outliers

outcomes <- c("t2ll_slope","bd_slope","gmsdis_slope","combi_slope","sc_slope")
outcomes_nice <- c("T2LL slope","BD slope","GMS disability slope","CombiWISE slope","SC disability slope")
outcomes_fu <- c("t2ll_fu_length","bd_fu_length","gmsdis_fu_length","combi_fu_length","sc_fu_length")
outcomes_nr_visits <- c("t2ll_nr_visits","bd_nr_visits","gmsdis_nr_visits","combi_nr_visits","sc_nr_visits")


slopes_clean <- remout_fun(dat = slopes,outcomes = outcomes,mult = 3)

# eliminate outlier values for fu_length
for(i in 1:length(outcomes_fu)){
  slopes_clean[[outcomes_fu[i]]] <- ifelse(is.na(slopes_clean[[outcomes[i]]]),yes = NA,no = slopes_clean[[outcomes_fu[i]]])
  slopes_clean[[outcomes_nr_visits[i]]] <- ifelse(is.na(slopes_clean[[outcomes[i]]]),yes = NA,no = slopes_clean[[outcomes_nr_visits[i]]])
}

write_csv(slopes_clean,"./output/longitudinal_analyses/Analysis_one_Outcome_slopes_data.csv")

# summarize SLOPES

slopes_summary <- slopes_clean %>%
  summarise(across(all_of(outcomes), summary)) %>% 
  mutate_at(.vars = outcomes,.funs = as.numeric) 
slopes_summary <- round(slopes_summary,digits = 2)


slopes_summary$Statistic <- names(summary(slopes_clean$combi_slope))

slopes_summary <- slopes_summary %>% 
  dplyr::select(Statistic,everything())

write_csv(slopes_summary,"./output/longitudinal_analyses/Analysis_one_Outcome_slopes_vs_Somamers_SLOPES_summary.csv")


# summarize fu_length

slopes_summary <- slopes_clean %>%
  summarise(across(all_of(outcomes_fu), summary)) %>% 
  mutate_at(.vars = outcomes_fu,.funs = as.numeric) 
slopes_summary <- round(slopes_summary,digits = 1)


slopes_summary$Statistic <- names(summary(slopes_clean$t2ll_fu_length))

slopes_summary <- slopes_summary %>% 
  dplyr::select(Statistic,everything())

write_csv(slopes_summary,"./output/longitudinal_analyses/Analysis_one_Outcome_slopes_vs_Somamers_fu_summary.csv")


# summarize nr_visits

slopes_summary <- slopes_clean %>%
  summarise(across(all_of(outcomes_nr_visits), summary)) %>% 
  mutate_at(.vars = outcomes_nr_visits,.funs = as.numeric) 
slopes_summary <- round(slopes_summary,digits = 1)


slopes_summary$Statistic <- names(summary(slopes_clean$t2ll_nr_visits))

slopes_summary <- slopes_summary %>% 
  dplyr::select(Statistic,everything())

write_csv(slopes_summary,"./output/longitudinal_analyses/Analysis_one_Outcome_slopes_vs_Somamers_nr_visits_summary.csv")


data_to_plot <- slopes_clean

list_of_plots <- list()

for(i in 1:length(outcomes)){
  # i <- 4
  
  data_to_plot <- slopes_clean
  
  tert <- quantile(x = slopes_clean[[outcomes[i]]],na.rm = TRUE,probs = seq(0,1,1/3))
  
  t1 <- as.numeric(tert[2])
  t2 <- as.numeric(tert[3])
  
  data_to_plot$tertile <- NA
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]<=t1)] <- "T1"
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]>t1 & data_to_plot[[outcomes[i]]]<t2)] <- "T2"
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]>=t2)] <- "T3"
  table(data_to_plot$tertile)
  
  
  # Perform the one-sample t.test
 ttest_result <- t.test(data_to_plot[[outcomes[i]]], mu = 0)
  p_value <- ttest_result$p.value
  
  p <- data_to_plot %>% 
    mutate(point="slope") %>%
    ggplot(aes(x=point,y=.data[[outcomes[i]]]))+
    geom_hline(yintercept = 0,linetype="dashed",color="gray",linewidth=0.45)+
    geom_jitter(aes(fill=tertile),shape=21,color="black",stroke=0.2,alpha=0.5, size=2, position = position_jitter(width = 0.1))+
    geom_boxplot(outliers = FALSE, fill="transparent")+
    geom_violin(linewidth=0.15, fill="transparent")+
    scale_fill_manual(values = c("blue","transparent","red"),breaks = c("T1","T2","T3"),"Slope tertile")+
    
    EnvStats::stat_n_text()+
    ggtitle(paste0("T-test p = ", format(p_value, scientific = TRUE,digits = 2)))+
    ylab(outcomes_nice[i])+
    theme_classic2()+
    theme(
      # legend.background = element_rect(fill="transparent"),
      plot.title = element_text(size = 10),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=12,colour="black",hjust = 0.5,vjust = .5),
      axis.title.y = element_text(color="black",size=14),
      axis.title.x = element_blank(),
      # legend.box.background = element_rect(colour = "black",fill="white"),
      legend.position = "none",
      axis.line = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      # panel.grid.major.y = element_line(color="gray",size=.25),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      axis.ticks.x = element_blank(),
      plot.background = element_blank())
p
  
  
list_of_plots[[paste0("p",i)]] <- p

}

# Extract the legend from one of the plots
legend <- ggpubr::get_legend(list_of_plots$p2 + theme(legend.position = "bottom"))

# Combine the plots without their legends
combined_plots <- grid.arrange(grobs = lapply(list_of_plots, function(p) p + theme(legend.position = "none")),
                               ncol = 5)

# Add the legend below the plots
final_plot <- arrangeGrob(combined_plots, legend, nrow = 2, heights = c(5, 1))

# Save the combined plot with the common legend
ggsave("./output/longitudinal_analyses/first_analysis_longitudinal_all_samples.png", plot = final_plot, width = 10, height = 3)



#########################################################################
#
# split cohort into tertiles and calculate differences 
# between first and third tertile for somamers
# measured at first LP that significantly correlate with outcome 
#
#########################################################################

sig_outcomes <- c("t2ll","brain_damage","GMSdis","combiwise","sc_disability")

analysis_one <- tibble::enframe(trans$marker,name = NULL,value = "marker")

for(i in 1:length(outcomes)){
  # i <- 4
  
  # isolate slopes for the outcome
  data_to_plot <- slopes_clean
  
  # calcluate tertiles
  tert <- quantile(x = slopes_clean[[outcomes[i]]],na.rm = TRUE,probs = seq(0,1,1/3))
  
  t1 <- as.numeric(tert[2])
  t2 <- as.numeric(tert[3])
  
  # identify samples in tertiles
  data_to_plot$tertile <- NA
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]<=t1)] <- "T1"
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]>t1 & data_to_plot[[outcomes[i]]]<t2)] <- "T2"
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]>=t2)] <- "T3"
  table(data_to_plot$tertile)
  
  # isolate T1 patients
  early <- data_to_plot %>% 
    filter(tertile=="T1")
  early_ptx <- early$patientcode
  
  # isolate T3 patients 
  late <- data_to_plot %>% 
    filter(tertile=="T3")
  late_ptx <- late$patientcode
  
  
  # generate vector of significant somamers 
    master_sel <- master %>% 
    filter(master[[paste0(sig_outcomes[i],"_coef_pval")]]<=0.05 & master[[paste0(sig_outcomes[i],"_pval_FDR")]]<=0.05)
  
  marker_sel <- master_sel$marker
  
  #isolate 1st untreated CSF sample in 1st and 3rd tertile cohort
  soma_sel <- soma %>% 
    filter(patientcode %in% c(early_ptx,late_ptx)) %>% 
    filter(sampleid %in% ms_first$sampleid)
  
  soma_sel$group <- ifelse(soma_sel$patientcode %in% early_ptx,yes = "t1",no = "t3")
  
  
# group differences

# t1 vs t3

differences <- tibble::enframe(marker_sel,name = NULL,value = "marker")
differences[[paste0(outcomes[i],"_t3_t1_median_diff")]] <- NA

for(l in 1:length(marker_sel)){
  # l <- 1
  differences[[paste0(outcomes[i],"_t3_t1_median_diff")]][l] <- median(subset(soma_sel,subset = group=="t3")[[marker_sel[l]]]) - median(subset(soma_sel,subset = group=="t1")[[marker_sel[l]]])
  
}

var <- "group"
dat <- soma_sel
wilc_p <- function(col_num){
  # col_num <- 1
  wilc <- wilcox.test(dat[[marker_sel[col_num]]]~dat[[var]], alternative = "two.sided")
  
  p <- tibble::enframe(wilc$p.value,name = NULL, value = marker_sel[col_num])
  out <- (p)
  
  return(out)
}


seq <- c(1:length(marker_sel))

fin <- bind_cols(lapply(1:length(seq), function(x) wilc_p(seq[x])))

final <- as.data.frame(t(fin))
final$marker <- marker_sel
final <- final %>% 
  select(marker,everything())
final$V3 <- p.adjust(final$V1,method = "fdr")


names(final) <- c("marker",paste0(outcomes[i],"_t3_t1_diff_Wilcox_raw"),paste0(outcomes[i],"_t3_t1_diff_Wilcox_FDR"))



results <- left_join(differences, final)

analysis_one <- left_join(analysis_one,results)

}

write_csv(analysis_one,"./output/longitudinal_analyses/Analysis_one_Outcome_slopes_vs_Somamers_atFirstLP_in_tertials_10152024.csv")




#********************************************************#
##########################################################
#
# ANALYSIS TWO
#
##########################################################
#********************************************************#

# isolate first untreated LP

demo_ms <- demo %>% 
  filter(diagnosis %in% c("RR-MS","SP-MS","PP-MS")) 

ms_first <- demo_ms %>% 
  filter(nice_therapy=="Untreated") %>% 
  group_by(patientcode) %>% 
  arrange(age) %>% 
  slice_head()

# identify last available LP
ms_last<- demo_ms %>%
  group_by(patientcode) %>%
  arrange(desc(age)) %>%
  slice_head() %>% 
  filter(patientcode %in% ms_first$patientcode)

ms_first$patientcode==ms_last$patientcode

# calculate age difference
ms_first$age_diff <- ms_last$age - ms_first$age

# filter out patients with no follow-up

ms_first_clean <- ms_first %>%
  filter(age_diff>1.5) %>%
  arrange(patientcode)

ms_last_clean <- ms_last %>%
  filter(patientcode %in% ms_first_clean$patientcode) %>%
  arrange(patientcode)


# check whether first and last have the same order of patients

ms_first_clean$patientcode==ms_last_clean$patientcode



# calculate vector of length of fu
length_fu <- as.numeric(( ms_last_clean$lpdate - ms_first_clean$lpdate ) /365.25)

hist(length_fu)

# calculate min, max, range

summary(length_fu)

# isolate samples of interest

ptx <- ms_first_clean$patientcode

clin <- neuro[0,]
t2 <- mri[0,]

slopes <- tibble::enframe(ptx,name = NULL,value = "patientcode")

slopes$combi_slope <- NA
slopes$combi_r2 <- NA
slopes$combi_fu_length <- NA
slopes$combi_nr_visits <- NA

slopes$bd_slope <- NA
slopes$bd_r2 <- NA
slopes$bd_fu_length <- NA
slopes$bd_nr_visits <- NA

slopes$sc_slope <- NA
slopes$sc_r2 <- NA
slopes$sc_fu_length <- NA
slopes$sc_nr_visits <- NA

slopes$gmsdis_slope <- NA
slopes$gmsdis_r2 <- NA
slopes$gmsdis_fu_length <- NA
slopes$gmsdis_nr_visits <- NA

slopes$t2ll_slope <- NA
slopes$t2ll_r2 <- NA
slopes$t2ll_fu_length <- NA
slopes$t2ll_nr_visits <- NA

neuro_outcomes <- c("combiwise","bd_disability","GMSDIS_impute","sc_disability")
slp <- c("combi_slope","bd_slope","gmsdis_slope","sc_slope")
rsqr <- c("combi_r2","bd_r2","gmsdis_r2","sc_r2")
fu <- c("combi_fu_length","bd_fu_length","gmsdis_fu_length","sc_fu_length")
nr_visits <- c("combi_nr_visits","bd_nr_visits","gmsdis_nr_visits","sc_nr_visits")

for(i in 1:length(ptx)){
  # i <- 1
  
  for(j in 1:length(neuro_outcomes)){
    # j <- 1
    
    # isolate samples of interest
    dat <- neuro %>% 
      filter(patientcode==ptx[i])%>% 
      filter(!is.na(.data[[neuro_outcomes[j]]])) %>% 
      filter(date >= ms_first_clean$dateneuroexam[i] & date <= ms_last_clean$dateneuroexam[i])
    
    
    # generate slopes
    if(dim(dat)[1]>1){
      if(max(dat$age)-min(dat$age)>=1.5){
        mod <- summary(lm(dat[[neuro_outcomes[j]]]~dat$age))
        slopes[[slp[j]]][i] <- mod$coefficients[2,1]
        slopes[[rsqr[j]]][i] <- mod$r.squared
        slopes[[fu[j]]][i] <- max(dat$age)-min(dat$age)
        slopes[[nr_visits[j]]][i] <- dim(dat)[1]
      }
    }
    
    else
    {
      slopes[[slp[j]]][i] <- NA
      slopes[[rsqr[j]]][i] <- NA
      slopes[[fu[j]]][i] <-NA
      slopes[[nr_visits[j]]][i] <- NA
    }
  }
  
  # add T2LL
  dat <- mri %>% 
    filter(patientcode==ptx[i]) %>% 
    filter(!is.na(volume_lesions)) %>% 
    filter(date >= ms_first_clean$mri_date[i] & date <= ms_last_clean$mri_date[i]) 
  
  if(dim(dat)[1]>1){
    if(max(dat$age)-min(dat$age)>=1.5){
      mod <- summary(lm(volume_lesions~age,dat))
      slopes$t2ll_slope[i] <- mod$coefficients[2,1]
      slopes$t2ll_r2[i] <- mod$r.squared
      slopes$t2ll_fu_length[i] <- max(dat$age)-min(dat$age)
      slopes$t2ll_nr_visits[i] <- dim(dat)[1]
    }
  }
  
  else
  {
    slopes$t2ll_slope[i] <- NA
    slopes$t2ll_r2[i] <- NA
    slopes$t2ll_fu_length[i] <-NA
    slopes$t2ll_nr_visits[i] <- NA
  }
  
}


## ELIMINATE +/- 3*IQR outliers

outcomes <- c("t2ll_slope","bd_slope","gmsdis_slope","combi_slope","sc_slope")
outcomes_nice <- c("T2LL slope","BD slope","GMS disability slope","CombiWISE slope","SC disability slope")
outcomes_fu <- c("t2ll_fu_length","bd_fu_length","gmsdis_fu_length","combi_fu_length","sc_fu_length")
outcomes_nr_visits <- c("t2ll_nr_visits","bd_nr_visits","gmsdis_nr_visits","combi_nr_visits","sc_nr_visits")


slopes_clean <- remout_fun(dat = slopes,outcomes = outcomes,mult = 3)

# eliminate outlier values for fu_length
for(i in 1:length(outcomes_fu)){
  slopes_clean[[outcomes_fu[i]]] <- ifelse(is.na(slopes_clean[[outcomes[i]]]),yes = NA,no = slopes_clean[[outcomes_fu[i]]])
  slopes_clean[[outcomes_nr_visits[i]]] <- ifelse(is.na(slopes_clean[[outcomes[i]]]),yes = NA,no = slopes_clean[[outcomes_nr_visits[i]]])
}


write_csv(slopes_clean,"./output/longitudinal_analyses/Analysis_two_Outcome_slopes_data.csv")

# summarize SLOPES

slopes_summary <- slopes_clean %>%
  summarise(across(all_of(outcomes), summary)) %>% 
  mutate_at(.vars = outcomes,.funs = as.numeric) 
slopes_summary <- round(slopes_summary,digits = 2)


slopes_summary$Statistic <- names(summary(slopes_clean$t2ll_slope))

slopes_summary <- slopes_summary %>% 
  dplyr::select(Statistic,everything())

write_csv(slopes_summary,"./output/longitudinal_analyses/Analysis_two_Outcome_longitudinal_firsVlast_samples_SLOPES_summary.csv")

# summarize fu_length

slopes_summary <- slopes_clean %>%
  summarise(across(all_of(outcomes_fu), summary)) %>% 
  mutate_at(.vars = outcomes_fu,.funs = as.numeric) 
slopes_summary <- round(slopes_summary,digits = 1)


slopes_summary$Statistic <- names(summary(slopes_clean$t2ll_fu_length))

slopes_summary <- slopes_summary %>% 
  dplyr::select(Statistic,everything())

write_csv(slopes_summary,"./output/longitudinal_analyses/Analysis_two_Outcome_longitudinal_firsVlast_samples_fu_summary.csv")


# summarize nr_visits

slopes_summary <- slopes_clean %>%
  summarise(across(all_of(outcomes_nr_visits), summary)) %>% 
  mutate_at(.vars = outcomes_nr_visits,.funs = as.numeric) 
slopes_summary <- round(slopes_summary,digits = 1)


slopes_summary$Statistic <- names(summary(slopes_clean$t2ll_nr_visits))

slopes_summary <- slopes_summary %>% 
  dplyr::select(Statistic,everything())

write_csv(slopes_summary,"./output/longitudinal_analyses/Analysis_two_Outcome_longitudinal_firsVlast_samples_nr_visits_summary.csv")




data_to_plot <- slopes_clean

list_of_plots <- list()

for(i in 1:length(outcomes)){
  # i <- 5
  
  data_to_plot <- slopes_clean
  
  tert <- quantile(x = slopes_clean[[outcomes[i]]],na.rm = TRUE,probs = seq(0,1,1/3))
  
  t1 <- as.numeric(tert[2])
  t2 <- as.numeric(tert[3])
  
  data_to_plot$tertile <- NA
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]<=t1)] <- "T1"
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]>t1 & data_to_plot[[outcomes[i]]]<t2)] <- "T2"
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]>=t2)] <- "T3"
  table(data_to_plot$tertile)
  
  
  # Perform the one-sample t.test
  ttest_result <- t.test(data_to_plot[[outcomes[i]]], mu = 0)
  p_value <- ttest_result$p.value
  
  p <- data_to_plot %>% 
    mutate(point="slope") %>%
    ggplot(aes(x=point,y=.data[[outcomes[i]]]))+
    geom_hline(yintercept = 0,linetype="dashed",color="gray",linewidth=0.45)+
    geom_jitter(aes(fill=tertile),shape=21,color="black",stroke=0.2,alpha=0.5, size=2, position = position_jitter(width = 0.1))+
    geom_boxplot(outliers = FALSE, fill="transparent")+
    geom_violin(linewidth=0.15, fill="transparent")+
    scale_fill_manual(values = c("blue","transparent","red"),breaks = c("T1","T2","T3"),"Slope tertile")+
    
    EnvStats::stat_n_text()+
    ggtitle(paste0("T-test p = ", format(p_value, scientific = TRUE,digits = 2)))+
    ylab(outcomes_nice[i])+
    theme_classic2()+
    theme(
      # legend.background = element_rect(fill="transparent"),
      plot.title = element_text(size = 10),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=12,colour="black",hjust = 0.5,vjust = .5),
      axis.title.y = element_text(color="black",size=14),
      axis.title.x = element_blank(),
      # legend.box.background = element_rect(colour = "black",fill="white"),
      legend.position = "none",
      axis.line = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      # panel.grid.major.y = element_line(color="gray",size=.25),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      axis.ticks.x = element_blank(),
      plot.background = element_blank())
  p
  
  
  list_of_plots[[paste0("p",i)]] <- p
  
}

# Extract the legend from one of the plots
legend <- ggpubr::get_legend(list_of_plots$p2 + theme(legend.position = "bottom"))

# Combine the plots without their legends
combined_plots <- grid.arrange(grobs = lapply(list_of_plots, function(p) p + theme(legend.position = "none")),
                               ncol = 5)

# Add the legend below the plots
final_plot <- arrangeGrob(combined_plots, legend, nrow = 2, heights = c(5, 1))

# Save the combined plot with the common legend
ggsave("./output/longitudinal_analyses/second_analysis_longitudinal_firsVlast_samples.png", plot = final_plot, width = 10, height = 3)




################################################
# calculate somamer and pathway yearly slopes 
################################################

high_markers <- readLines("./input/high_SNR_markers.txt")

org <- read_csv("./input/organism.csv")

org <- org %>% 
  filter(organism=="Human")

high_markers <- high_markers[which(high_markers%in%org$marker)]  


high_pathway <- names(path)[-1]

markers <- c(high_markers,high_pathway)



##################!!!!!!!!!!!!!!!!!!!!!!##################################
#
#.  THIS STEP TAKES TIME -> RUN ONCE AND THEN JUST UPLOAD RESULTS
#
##################!!!!!!!!!!!!!!!!!!!!!!!#################################

# generate patient-specific slopes for relevant markers (those that correlate 
# cross-sectionally with at least one outcome) and for patients that have 
# relevant slope (two LPs with at least 1.5y apart)

##############################################################
# generate a list of relevant patients and their sample id
# (all LP between first and last LP)

rel_ptx <- vector()

for(i in 1:length(outcomes)){
  # i <- 1
  dat <- slopes_clean %>%
    filter(!is.na(slopes_clean[[outcomes[i]]]))
  rel_ptx <- c(rel_ptx,dat$patientcode)
  }

# select only unique patients
rel_ptx <- sort(unique(rel_ptx))

# clean ms_first and ms_last for relevant patients only

ms_first_rel <- ms_first_clean %>%
  filter(patientcode %in% rel_ptx) %>%
  arrange(patientcode)

ms_last_rel <- ms_last_clean %>%
  filter(patientcode %in% rel_ptx) %>%
  arrange(patientcode)

# check whether first and last match
ms_first_rel$patientcode==rel_ptx
ms_last_rel$patientcode==rel_ptx

# isolate samples between first and last LP

rel_samples <- soma[0,]

for(i in 1:length(rel_ptx)){
  # i <- 91
dat <- soma %>%
  filter(patientcode==rel_ptx[i]) %>%
  filter(lpdate>=ms_first_rel$lpdate[i] & lpdate<=ms_last_rel$lpdate[i])

rel_samples <- rbind(rel_samples,dat)
}

# isolate relevant markers

# generate vector of significant somamers

# generate vector of significant somamers

rel_markers <- vector()
for(i in 1:length(sig_outcomes)){

    master_sel <- master %>%
    filter(master[[paste0(sig_outcomes[i],"_coef_pval")]]<=0.05 & master[[paste0(sig_outcomes[i],"_pval_FDR")]]<=0.05)
rel_markers <- c(rel_markers,master_sel$marker)
}
# isolate unique markers
rel_markers <- unique(rel_markers)

##########################################
# generate patient-specific slope - this was ChatGPT-generated

get_patient_specific_coefficients <- function(df, subject_col, age_col, variables) {
  # Create an empty list to store results
  results_list <- list()

  # Get unique subject IDs
  subject_ids <- unique(df[[subject_col]])

  total_subjects <- length(subject_ids)

  # Loop over each subject
  for (i in seq_along(subject_ids)) {
    subject <- subject_ids[i]

    # Print progress update
    cat(sprintf("Processing subject %d of %d\n", i, total_subjects))

    # Subset data for the current subject
    subject_data <- df %>% filter(.[[subject_col]] == subject)

    # Initialize a list to store coefficients for the current subject
    subject_coeffs <- setNames(numeric(length(variables)), variables)

    # Loop over each variable
    for (var in variables) {
      # Fit the linear model
      model <- lm(as.formula(paste0("`",var,"`","~", age_col)), data = subject_data)

      # Extract the age coefficient
      if ("age" %in% names(coef(model))) {
        subject_coeffs[var] <- coef(model)["age"]
      } else {
        subject_coeffs[var] <- NA
      }
    }

    # Store the coefficients in the results list
    results_list[[as.character(subject)]] <- subject_coeffs
  }

  # Combine results into a data frame
  results_df <- do.call(rbind, lapply(results_list, as.data.frame))

  # Add subject ID column
  results_df <- cbind(Subject_ID = rownames(results_df), results_df)

  return(results_df)
}


# Example variables list, replace with your actual variable names
variables <- rel_markers  # or specify your variable names explicitly

# Your DataFrame
df <- rel_samples

# Get patient-specific coefficients
patient_coeffs <- get_patient_specific_coefficients(df, subject_col = "patientcode", age_col = "age", variables = variables)

# View the results
head(patient_coeffs)

patient_coeffs_clean <- patient_coeffs %>%
  separate(Subject_ID, into = c("patientcode","markers"),sep="\\.")

# rename outcome
names(patient_coeffs_clean)[3] <- "slope"

patient_coeffs_clean_wide <- patient_coeffs_clean %>%
  pivot_wider(names_from = markers, values_from = slope)

write_csv(patient_coeffs_clean_wide,"./output/longitudinal_analyses/patient_specific_slopes.csv")

# patient_coeffs_clean_wide <- read_csv("./output/longitudinal_analyses/patient_specific_slopes.csv")

#########################################################################
#
# split cohort into tertiles and calculate differences 
# between first and third tertile for somamer and pathway yearly slopes 
# that significantly correlate with outcome 
#
#########################################################################

sig_outcomes <- c("t2ll","brain_damage","GMSdis","combiwise","sc_disability")

analysis_two <- tibble::enframe(trans$marker,name = NULL,value = "marker")

for(i in 1:length(outcomes)){
  # i <- 1
  
  data_to_plot <- slopes_clean
  
  tert <- quantile(x = slopes_clean[[outcomes[i]]],na.rm = TRUE,probs = seq(0,1,1/3))
  
  t1 <- as.numeric(tert[2])
  t2 <- as.numeric(tert[3])
  
  data_to_plot$tertile <- NA
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]<=t1)] <- "T1"
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]>t1 & data_to_plot[[outcomes[i]]]<t2)] <- "T2"
  data_to_plot$tertile[which(data_to_plot[[outcomes[i]]]>=t2)] <- "T3"
  table(data_to_plot$tertile)
  
  early <- data_to_plot %>% 
    filter(tertile=="T1")
  early_ptx <- early$patientcode
  
  late <- data_to_plot %>% 
    filter(tertile=="T3")
  late_ptx <- late$patientcode
  
  
  # generate vector of significant somamers 
  
  master_sel <- master %>% 
    filter(master[[paste0(sig_outcomes[i],"_coef_pval")]]<=0.05 & master[[paste0(sig_outcomes[i],"_pval_FDR")]]<=0.05) 
  
  marker_sel <- master_sel$marker
  
  soma_sel <- patient_coeffs_clean_wide %>% 
    filter(patientcode %in% c(early_ptx,late_ptx))
  
  soma_sel$group <- ifelse(soma_sel$patientcode %in% early_ptx,yes = "t1",no = "t3")
  
  
  # group differences
  
  # t1 vs t3
  
  differences <- tibble::enframe(marker_sel,name = NULL,value = "marker")
  differences[[paste0(outcomes[i],"_change_t3_t1_median_diff")]] <- NA
  
  for(l in 1:length(marker_sel)){
    # l <- 1
    differences[[paste0(outcomes[i],"_change_t3_t1_median_diff")]][l] <- median(subset(soma_sel,subset = group=="t3")[[marker_sel[l]]]) - median(subset(soma_sel,subset = group=="t1")[[marker_sel[l]]])
    
  }
  
  var <- "group"
  dat <- soma_sel
  wilc_p <- function(col_num){
    # col_num <- 1
    wilc <- wilcox.test(dat[[marker_sel[col_num]]]~dat[[var]], alternative = "two.sided")
    
    p <- tibble::enframe(wilc$p.value,name = NULL, value = marker_sel[col_num])
    out <- (p)
    
    return(out)
  }
  
  
  seq <- c(1:length(marker_sel))
  
  fin <- bind_cols(lapply(1:length(seq), function(x) wilc_p(seq[x])))
  
  final <- as.data.frame(t(fin))
  final$marker <- marker_sel
  final <- final %>% 
    select(marker,everything())
  final$V3 <- p.adjust(final$V1,method = "fdr")
  
  
  names(final) <- c("marker",paste0(outcomes[i],"_change_t3_t1_diff_Wilcox_raw"),paste0(outcomes[i],"_change_t3_t1_diff_Wilcox_FDR"))
  
  
  
  results <- left_join(differences, final)
  
  analysis_two <- left_join(analysis_two,results)
  
}

write_csv(analysis_two,"./output/longitudinal_analyses/Analysis_TWO_Outcome_slopes_vs_SomamerSLopes_in_tertials.csv")

