## 0) Script description:#######################################################
#
# Date: 24/03/2022
# Author: Pol Jimenez-Arenas
# Topic: PlaDNAm-neuro analysis
# Objective: To select variables, load and preprocess DNAm data, generate
# descriptives and run analysis of the PlaDNAm-neuro analysis in the yamabuki server
##

##
# version: v.1.0
# latest version date: 24/03/2023
##

##
# changelist: NA
##


## 1) Description of the script:################################################

# The following R code will allow you to complete all the EWAS requested in the 
# PACE ext/int behavior and placenta DNA methylation analysis plan.
# The code also produces files summarising the variables included in the EWAS.
# You shouldn't have to rewrite or add to the following code, unless otherwise stated.
#
# If the following variables are named differently in your dataset, please rename the 
# variables accordingly details on how to code these variables are provided in the analysis plan.


## 2) Variables for me:#########################################################

# There are just two inputs required for this analysis:
# 1. phenodataframe: a dataframe containing all the "phenotype" data needed for this project. 
#    Each row is a sample (individual) and each column is a different variable: 
#    - variables of interest: 
#             - "pt_igd"
#             - "pt_motr"
#             - "pt_cadapt"
#             - "pt_socioem"
#             - "pt_cogn"
#             - "pt_com"
#    - Main covariates:
#             - "Sex" (from DNAm data)
#             - "ethnicity_c_4cat.y"
#             - "gestage_0y_c_weeks_fur"
#             - "edad"
#             - "parity_m_2cat"
#             - "smoke_sust_m"
#             - "educ_level_m_3cat"
#             - "CELL TYPE PROP"
#             - "hospital_del_pa"
#             - "covid_confinement_m"
#             - "age_dp3"
#    - Optional covariates:
#             - "MAT_MENTALH"
#             - BREASTFEED_6M
#    - Variables for sensitivity analysis: pregnancy_complications (presence of preeclampsia, diabetes 2, or other conditions considered by cohorts that
#                                          can imply pregnancy complications), preterm (gestational age<37 or >37), low_bw (low birthweight, <2500g)
#
# 2. Betasnooutliers: a matrix of methylation illumina beta values. Each column is a 
#    sample and each row is a probe on the array (450k or EPIC). Column names must correspond 
#    to the SampleID ("ID") in the phenodataframe.


## 3) Set wd and load packages:#################################################

rm(list=ls())

setwd("/PROJECTES/BISC_OMICS/analyses/BiSC_22/031_EWASpla_Neuro_PJ/results/") # Set working directory: go to the directory where to save results

library(PACEanalysis)


## 4) Data loading:#############################################################

# If your data have been preprocessed according to PACEanalysis package instructions, 
# then proceed to the next step (Data preparation). Otherwise, please follow the 
# preprocessing steps found in https://www.epicenteredresearch.com/pace/birthsize/

# Check your data
# If you encounter any issues, please check out our troubleshooting guide to see if there is any guidance that may help: https://www.epicenteredresearch.com/pace/troubleshooting/
# If you closed prior R session (step 1), you can load list of pre-processed objects that is automatically saved by the preprocessingofData function, e.g.

## METHYLATION DATA ##########

load("/PROJECTES/BISC_OMICS/data_final/methyl/child/0y/placenta_EPIC_QCPACE_20220823/BISC_withoutSiblings_20230201_PreprocessedBetas_nooutliers.RData")

## PHENODATAFRAME ############

# Metadata de les dades de metilacio: ID, sexe segons metilacio, basename (ID de metilacio) i contaminacio

pheno <- read.table("/PROJECTES/BISC_OMICS/data_final/methyl/child/0y/placenta_EPIC_QCPACE_20220823/BISC_withoutSiblings_20230201_Metadata.txt", 
                             sep="\t", header=TRUE) 
pheno <- as.data.frame(pheno[, c("id_participant.c.num", "Sex", "Basename", "Meanlog2oddsContamination")]) # subset de l'anterior nomes amb vars d'interes

# DP-3 i covars

load("/PROJECTES/BISC_OMICS/analyses/BiSC_22/031_EWASpla_Neuro_PJ/db/pheno_pt.RData") 

phenodataframe <- merge(pheno, phenodataframe, by.x = "id_participant.c.num", by.y = "id_bisc")

rm(pheno)

# Cell types

celltypes <- read.table("/PROJECTES/BISC_OMICS/data_final/methyl/child/0y/placenta_EPIC_QCPACE_20220823/BISC_withoutSiblings_20230201_Celltypes.txt", header = TRUE)

celltypes <- celltypes[rownames(celltypes)%in%phenodataframe$Basename,] # nomes dels IDs que tinc


## 5) Data preprocessing:#######################################################

# ensure var type

str(phenodataframe)
phenodataframe$id_participant.c.num <- as.character(phenodataframe$id_participant.c.num)
phenodataframe$Sex <- as.factor(phenodataframe$Sex)
phenodataframe$Basename <- as.character(phenodataframe$Basename)
phenodataframe$Meanlog2oddsContamination <- as.numeric(phenodataframe$Meanlog2oddsContamination)
phenodataframe$pt_igd <- as.numeric(phenodataframe$pt_igd)
phenodataframe$pt_motr <- as.numeric(phenodataframe$pt_motr)
phenodataframe$pt_cadapt <- as.numeric(phenodataframe$pt_cadapt)
phenodataframe$pt_socioem <- as.numeric(phenodataframe$pt_socioem)
phenodataframe$pt_cogn <- as.numeric(phenodataframe$pt_cogn)
phenodataframe$pt_com <- as.numeric(phenodataframe$pt_com)
phenodataframe$gestage_0y_c_weeks_fur <- as.numeric(phenodataframe$gestage_0y_c_weeks_fur)
phenodataframe$edad <- as.numeric(phenodataframe$edad)
phenodataframe$parity_m_2cat <- as.factor(phenodataframe$parity_m_2cat)
phenodataframe$smoke_sust_m_2cat <- as.factor(phenodataframe$smoke_sust_m)
phenodataframe$educ_level_m_3cat <- as.factor(phenodataframe$educ_level_m_3cat)
phenodataframe$covid_confinement_m_3cat <- as.factor(phenodataframe$covid_confinement_m)
phenodataframe$age_dp3 <- as.numeric(phenodataframe$age_dp3)
phenodataframe$dummy_igs <- as.factor(phenodataframe$dummy_igs)
phenodataframe$ethnicity_c_3cat <- as.factor(phenodataframe$ethnicity_c_3cat)
phenodataframe$hosp_part <- as.factor(phenodataframe$hosp_part)


## 6) Descriptives:#############################################################

# Calculate descriptives for your variables of interest: 

mean_pt_igd <- mean(phenodataframe$pt_igd) 
sd_pt_igd <-  sd(phenodataframe$pt_igd)
min_pt_igd<- min(phenodataframe$pt_igd)
max_pt_igd<- max(phenodataframe$pt_igd)

mean_pt_motr <- mean(phenodataframe$pt_motr) 
sd_pt_motr <-  sd(phenodataframe$pt_motr)
min_pt_motr <- min(phenodataframe$pt_motr)
max_pt_motr <- max(phenodataframe$pt_motr)

mean_pt_cadapt <- mean(phenodataframe$pt_cadapt)
sd_pt_cadapt <-  sd(phenodataframe$pt_cadapt)
min_pt_cadapt <- min(phenodataframe$pt_cadapt)
max_pt_cadapt <- max(phenodataframe$pt_cadapt)

mean_pt_socioem <- mean(phenodataframe$pt_socioem) 
sd_pt_socioem <-  sd(phenodataframe$pt_socioem)
min_pt_socioem<- min(phenodataframe$pt_socioem)
max_pt_socioem<- max(phenodataframe$pt_socioem)

mean_pt_cogn <- mean(phenodataframe$pt_cogn) 
sd_pt_cogn <-  sd(phenodataframe$pt_cogn)
min_pt_cogn <- min(phenodataframe$pt_cogn)
max_pt_cogn <- max(phenodataframe$pt_cogn)

mean_pt_com <- mean(phenodataframe$pt_com)
sd_pt_com <-  sd(phenodataframe$pt_com)
min_pt_com <- min(phenodataframe$pt_com)
max_pt_com <- max(phenodataframe$pt_com)

Tabledescriptives <- rbind(mean_pt_igd, sd_pt_igd, min_pt_igd, max_pt_igd, 
                           mean_pt_motr, sd_pt_motr, min_pt_motr, max_pt_motr, 
                           mean_pt_cadapt, sd_pt_cadapt, min_pt_cadapt, max_pt_cadapt, 
                           mean_pt_socioem, sd_pt_socioem, min_pt_socioem, max_pt_socioem,
                           mean_pt_cogn, sd_pt_cogn, min_pt_cogn, max_pt_cogn, 
                           mean_pt_com, sd_pt_com, min_pt_com, max_pt_com)

Tabledescriptives <- round(Tabledescriptives,2)
Tabledescriptives 

# Save results

cohort <- "BiSC"
write.csv(Tabledescriptives, file=paste0("PACE_Pla_DP3", cohort, 
                                         "_Descriptives.var_", 
                                         format(Sys.Date(), "%d%m%Y"), ".csv"))

# Histograms

pdf(file=paste0(cohort,"_pt_igd_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_igd)
dev.off()

pdf(file=paste0(cohort,"_pt_motr_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_motr)
dev.off()

pdf(file=paste0(cohort,"_pt_cadapt_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_cadapt)
dev.off()

pdf(file=paste0(cohort,"_pt_socioem_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_socioem)
dev.off()

pdf(file=paste0(cohort,"_pt_cogn_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_cogn)
dev.off()

pdf(file=paste0(cohort,"_pt_com_Histogram_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_com)
dev.off()

# Complete cases: as we cannot have NAs in the variables of the principal model

phenodataframe <- phenodataframe[complete.cases(phenodataframe), ]

# Calculate descriptives for your variables of interest ONLY FOR COMPLETE CASES 

mean_pt_igd <- mean(phenodataframe$pt_igd) 
sd_pt_igd <-  sd(phenodataframe$pt_igd)
min_pt_igd<- min(phenodataframe$pt_igd)
max_pt_igd<- max(phenodataframe$pt_igd)

mean_pt_motr <- mean(phenodataframe$pt_motr) 
sd_pt_motr <-  sd(phenodataframe$pt_motr)
min_pt_motr <- min(phenodataframe$pt_motr)
max_pt_motr <- max(phenodataframe$pt_motr)

mean_pt_cadapt <- mean(phenodataframe$pt_cadapt)
sd_pt_cadapt <-  sd(phenodataframe$pt_cadapt)
min_pt_cadapt <- min(phenodataframe$pt_cadapt)
max_pt_cadapt <- max(phenodataframe$pt_cadapt)

mean_pt_socioem <- mean(phenodataframe$pt_socioem) 
sd_pt_socioem <-  sd(phenodataframe$pt_socioem)
min_pt_socioem<- min(phenodataframe$pt_socioem)
max_pt_socioem<- max(phenodataframe$pt_socioem)

mean_pt_cogn <- mean(phenodataframe$pt_cogn) 
sd_pt_cogn <-  sd(phenodataframe$pt_cogn)
min_pt_cogn <- min(phenodataframe$pt_cogn)
max_pt_cogn <- max(phenodataframe$pt_cogn)

mean_pt_com <- mean(phenodataframe$pt_com)
sd_pt_com <-  sd(phenodataframe$pt_com)
min_pt_com <- min(phenodataframe$pt_com)
max_pt_com <- max(phenodataframe$pt_com)

Tabledescriptives_cc <- rbind(mean_pt_igd, sd_pt_igd, min_pt_igd, max_pt_igd, 
                           mean_pt_motr, sd_pt_motr, min_pt_motr, max_pt_motr, 
                           mean_pt_cadapt, sd_pt_cadapt, min_pt_cadapt, max_pt_cadapt, 
                           mean_pt_socioem, sd_pt_socioem, min_pt_socioem, max_pt_socioem,
                           mean_pt_cogn, sd_pt_cogn, min_pt_cogn, max_pt_cogn, 
                           mean_pt_com, sd_pt_com, min_pt_com, max_pt_com)

Tabledescriptives_cc <- round(Tabledescriptives_cc,2)
Tabledescriptives_cc 

# Save results

cohort <- "BiSC"
write.csv(Tabledescriptives_cc, file=paste0("PACE_Pla_DP3", cohort, 
                                         "_Descriptives.var_cc_", 
                                         format(Sys.Date(), "%d%m%Y"), ".csv"))

# Histograms

pdf(file=paste0(cohort,"_pt_igd_Histogram_cc_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_igd)
dev.off()

pdf(file=paste0(cohort,"_pt_motr_Histogram_cc_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_motr)
dev.off()

pdf(file=paste0(cohort,"_pt_cadapt_Histogram_cc_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_cadapt)
dev.off()

pdf(file=paste0(cohort,"_pt_socioem_Histogram_cc_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_socioem)
dev.off()

pdf(file=paste0(cohort,"_pt_cogn_Histogram_cc_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_cogn)
dev.off()

pdf(file=paste0(cohort,"_pt_com_Histogram_cc_", format(Sys.Date(), "%d%m%Y"),".pdf"))
hist(phenodataframe$pt_com)
dev.off()

save(phenodataframe,file="BiSC_phenofinal_DP320230330.Rdata") ## edit with your cohort name

###*Correlations between main variables and save them:################

## For data sets with continuous, polytomous and dichotomous variables.

# Select relevant variables

# relevant_variables <- c("Sex","Meanlog2oddsContamination","gestage_0y_c_weeks_fur","edad",
#                         "parity_m_2cat", "smoke_sust_m_2cat", "educ_level_m_3cat", "covid_confinement_m_3cat",
#                         "age_dp3", "ethnicity_c_3cat", "hosp_part", "pt_igd", "pt_motr", 
#                         "pt_cadapt", "pt_socioem", "pt_cogn", "pt_cogn", "pt_com") #DP3's scores must be the last columns!

# column_positions <- match(relevant_variables, colnames(phenodataframe))
# corrdata <- phenodataframe[column_positions]
# colnames(corrdata)

# Check data type and coerce to numeric (only for variables that not numeric and continuous as indicated)
#NOTE: do not call as.numeric() directly since as.numeric() gives the internal codes 
# str(corrdata)

# corrdata$mat_edu <- as.numeric(as.character(corrdata$mat_edu))
# corrdata$mat_smk <- as.numeric(as.character(corrdata$mat_smk))
# corrdata$cesarean <- as.numeric(as.character(corrdata$cesarean))
# corrdata$mat_health <- as.numeric(as.character(corrdata$mat_health))


#Calculate correlations  

# source("mycor.ci.r") #Before calculating correlations, run the attached separate script "mycor.ci" (a modified version of cor.ci from psych package, supplied with this example code)

# continuous_variables <- c("mat_age","ges_age","child_age", "unint","unext", "untotal") ## use the untransformed variables
# polychoric_variables <- c("mat_edu","mat_smk") # edit according to your binary/dichotomous or polychoric covariates (for example some cohort might have mat_smk as binary)
# dichotomous_variables <- c("cesarean", "mat_health") # as before, edit accordingly

# c_vars_positions <- match(continuous_variables, colnames(corrdata))
# p_vars_positions <- match(polychoric_variables, colnames(corrdata))
# d_vars_positions <- match(dichotomous_variables, colnames(corrdata))

# corrdata<-as.matrix(corrdata)

# correlations <- mycor.ci(x=corrdata, keys = NULL, n.iter = 1000, p = 0.05, overlap = FALSE, 
# poly = FALSE, method = "pearson", plot = FALSE, minlength = 5, 
# cvars= c_vars_positions, 
# pvars= p_vars_positions, 
# dvars= d_vars_positions)
# correlations

#Get table of correlation coefficients
# corr_table <- correlations["rho"]
# corr_table

#Get table of confidence intervals and only keep those with p<.0.05
# sig_table <- data.frame(correlations["ci"]) %>% mutate(Correlation=row.names(.)) %>% filter(ci.p<.05)
# sig_table

#Write results files and save as .csv file

# write.csv(corr_table, file=paste0(cohort,"_intexttotal_correlations_r_", 
#   format(Sys.Date(), "%d%m%Y"),".csv"), col.names=FALSE, row.names=TRUE, quote=FALSE)

# write.csv(sig_table, file=paste0(cohort,"_intexttotal_correlations_ci_", 
#       format(Sys.Date(), "%d%m%Y"),".csv"), col.names=FALSE, row.names=TRUE, quote=FALSE)


## 7) Data analysis:############################################################

# This stage of the analysis is specific to the chosen exposure/outcome and the specified adjustment variables. 
# Below is the code for all of the analyses to run for the mental health project. 
# Please be sure to update the cohort and date information in the below code for your analysis, as well as the destination path. 
#
# Finally, be sure to update the column names of the exposure/outcome(s) of interest, the adjustment variables, 
# and the table 1 variables. These should correspond to column names in the dataframe specified in the phenofinal 
# argument of the dataAnalysis function. Quick check to make sure the function runs in your cohort

# Given the modeling approaches used, the dataAnalysis function requires a good deal of time to run. 
# We recommend first checking whether the function runs on a relatively small subset of sites (i.e. 100 CpG loci), as follows:
# betafinal=Betasnooutliers[1:100,], 
#
# If you encounter any issues, please let us know. If not, proceed to the next step.


## 6 models, 5 dimensions and GDS #####

# Main models NOT stratified by sex

modelstorun <- data.frame(varofinterest=c("pt_igd", "pt_motr", "pt_cadapt", "pt_socioem", "pt_cogn", "pt_com"))
modelstorun$varofinterest <- as.character(modelstorun$varofinterest)
modelstorun$vartype <- "ExposureCont"
Betasnooutliers <- betafinal.nooutlier ## Change betafinal.nooutlier to the name that appears when you load the methylation data object
Betas <- Betasnooutliers[,colnames(Betasnooutliers)%in%phenodataframe$Basename] # in case your Betas dataframe has more individuals than your phenodataframe, filter them so that they can be matched
dim(Betas) # 273 individuals: OK
Betasnooutliers <- Betas

# If not matched, then match them as follows: 

phenodataframe <- phenodataframe[match(colnames(Betasnooutliers),phenodataframe$Basename),]
all(phenodataframe$Basename == colnames(Betasnooutliers)) # TRUE

# Same with cell types. To match them we have to use a dataframe but we need to run the analysis witht them as a matrix:

celltypes_new2 <- as.data.frame(celltypes)

# Match celltypes with phenodataframe: 

phenodataframe <- phenodataframe[match(rownames(celltypes_new2),phenodataframe$Basename),]
all(phenodataframe$Basename == rownames(celltypes_new2)) # TRUE

rm(celltypes_new2)
rm(celltypes)

celltypes <- read.table("/PROJECTES/BISC_OMICS/data_final/methyl/child/0y/placenta_EPIC_QCPACE_20220823/BISC_withoutSiblings_20230201_Celltypes.txt", header = TRUE)

celltypes <- celltypes[rownames(celltypes)%in%phenodataframe$Basename,] # nomes dels IDs que tinc

celltypes_new2 <- as.data.frame(celltypes)

phenodataframe <- phenodataframe[match(rownames(celltypes_new2),phenodataframe$Basename),]
all(phenodataframe$Basename == rownames(celltypes_new2)) # TRUE

celltypes_new2 <- as.matrix(celltypes_new2)

# Then we check that when matching celltypes and phenodataframe we have not disorganised the data that was previously matched:

all(phenodataframe$Basename == colnames(Betasnooutliers)) # TRUE

# And finally we transform celltypes to matrix for the analysis:

celltypes2 <- as.matrix(celltypes_new2)

### *Analysis:#############################################

for (i in 1:nrow(modelstorun)){
  
  cat("Exposure:",modelstorun$varofinterest[i],"\n")
  
  tempresults<-dataAnalysis(phenofinal=phenodataframe,
                            betafinal=Betasnooutliers, 
                            array="EPIC", ## edit if you have 450K data
                            maxit=100,
                            robust=TRUE,
                            Omega=celltypes2,
                            vartype=modelstorun$vartype[i],
                            varofinterest=modelstorun$varofinterest[i],
                            Table1vars=c("Sex", "ethnicity_c_3cat", "gestage_0y_c_weeks_fur", "edad", "parity_m_2cat", 
                                         "smoke_sust_m_2cat", "educ_level_m_3cat", "hosp_part", "covid_confinement_m_3cat", 
                                         "age_dp3", "Meanlog2oddsContamination"), 
                            # remember to add optional variables such as "sel", "batch", "anc" or ancestry "PCs" if appropiate for your cohort 
                            StratifyTable1=FALSE,
                            StratifyTable1var=NULL,             
                            adjustmentvariables=c("Sex", "ethnicity_c_3cat", "gestage_0y_c_weeks_fur", "edad", "parity_m_2cat", 
                                                  "smoke_sust_m_2cat", "educ_level_m_3cat", "hosp_part", "covid_confinement_m_3cat", 
                                                  "age_dp3", "Meanlog2oddsContamination"), 
                            # remember to add optional variables such as "sel", "batch", "anc" or ancestry "PCs" if appropiate for your cohort 
                            RunUnadjusted=FALSE,
                            RunAdjusted=TRUE,
                            RunCellTypeAdjusted=TRUE,
                            RunSexSpecific=FALSE,
                            RestrictToSubset = FALSE,
                            RestrictionVar = NULL,
                            RestrictToIndicator = NULL,
                            RunCellTypeInteract = FALSE,
                            destinationfolder=setwd("/PROJECTES/BISC_OMICS/analyses/BiSC_22/031_EWASpla_Neuro_PJ/results/") ,  # EDIT with your directory
                            savelog=TRUE,
                            cohort="BISC", # EDIT with your cohort name
                            analysisdate = "20230417",# EDIT
                            analysisname = "testMAIN_pt")
}
