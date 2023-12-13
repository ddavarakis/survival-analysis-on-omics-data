
# Assignment: Machine Learning Mini-project

# Dimitrios – Tilemachos Davarakis

# UCD student ID: 22200956

suppressPackageStartupMessages({
  library(survminer)
  library(dplyr)
  library(ggplot2)
  library(ggsurvfit)
  library(gtsummary)
  library("glmnet")
  library("glmpath")
  library("glmnet")
  library(survival)
  library(survminer)
  library(lubridate)
  library(tidycmprsk)
})
#set directory
setwd("C:\\Users\\mtheo\\Documents\\R Working Dir\\P Medicine\\Survival")

# Dataset Used: Lung Adenocarcinoma dataset (TCGA, PanCancer Atlas, cell 2108)
# Refer to https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018

folder <- "luad_tcga_pan_can_atlas_2018"

### Read clinical patient data
filename = paste0(folder,"\\data_clinical_patient.txt")
clinical_patient_data <- read.delim(filename, skip = 4)
rownames(clinical_patient_data) <- clinical_patient_data$PATIENT_ID

### Read clinical sample data
filename = paste0(folder,"\\data_clinical_sample.txt")
clinical_sample_data <- read.delim(filename, skip = 4)
rownames(clinical_sample_data) <- clinical_sample_data$PATIENT_ID

# get from Sample data Sample_ID and Tumor_Type columns
clinical_data = merge(clinical_patient_data, clinical_sample_data[,c("PATIENT_ID", "SAMPLE_ID","TUMOR_TYPE")], by="PATIENT_ID")
# replace '-' to '.' in SAMPLE_ID column
clinical_data$SAMPLE_ID<-gsub("-",".",as.character(clinical_data$SAMPLE_ID))
# Preprocess values related to Progression Free Survival  
# Remove na values!
clinical_data <- clinical_data[!(is.na(clinical_data$PFS_MONTHS)), ]
clinical_data <- clinical_data[!(is.na(clinical_data$PFS_STATUS)), ]
# remove <=0 values from PFS_MONTHS
clinical_data <- clinical_data[clinical_data$PFS_MONTHS > 0, ]
clinical_data$PFS_STATUS <- ifelse (clinical_data$PFS_STATUS == "1:PROGRESSION", 1, 0 )

# Preprocess values related to Overall Survival 
# Remove na values!
clinical_data <- clinical_data[!(is.na(clinical_data$OS_MONTHS)), ]
clinical_data <- clinical_data[!(is.na(clinical_data$OS_STATUS)), ]
# remove <=0 values from OS_MONTHS
clinical_data <- clinical_data[clinical_data$OS_MONTHS > 0, ]
clinical_data$OS_STATUS <- ifelse (clinical_data$OS_STATUS == "1:DECEASED", 1, 0 )

rownames(clinical_data) <- clinical_data$SAMPLE_ID
# Confirm than no NA values exists
count<-which(is.na(clinical_data$PFS_MONTHS))
count
count<-which(is.na(clinical_data$PFS_STATUS))
count
count<-which(is.na(clinical_data$OS_MONTHS))
count
count<-which(is.na(clinical_data$OS_STATUS))
count

### Load and pre-process DNA data
load_rna_data <- function() {
  ### Read RNA expression data of Case samples
  filename = paste0(folder,"\\data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
  # 
  RNA_data <- read.delim(filename)

  ## Pre-processing
  # Check how many empty Hugo_Symbol exists 
  print(paste0("Number of genes with empty value: ", sum(RNA_data$Hugo_Symbol == "")))
  # Replace empty string with NA 
  RNA_data["Hugo_Symbol"][RNA_data["Hugo_Symbol"] == ''] <- NA
  # Remove NA 
  RNA_data <- RNA_data[!(is.na(RNA_data$Hugo_Symbol)), ]
  # Remove Entrez_Gene_id column
  RNA_data <- RNA_data[,-2]
  # check for duplicates
  geneNames <- RNA_data$Hugo_Symbol
  # find duplicates in Hugo_Symbol column, if any
  table(duplicated(geneNames))
  # There are duplicates
  # in order to resolve the duplication issue
  # list genes that are duplicate
  duplicated_genes <- geneNames[duplicated(geneNames) == TRUE]
  print("Duplicates genes are: ")
  print(duplicated_genes)
  # Rename the 2nd occurrence of each duplicated gene as : gene_2
  for (dup_gene in duplicated_genes) {
    res <- grep(dup_gene, geneNames)
    geneNames[res[2]] <- paste0(dup_gene,"_2")
    
  }
  # Recheck if the duplication issue resolved:
  duplicated_genes <- geneNames[duplicated(geneNames) == TRUE]
  duplicated_genes
  row.names(RNA_data) <- geneNames
  # transpose dataframe
  RNA_data <- as.data.frame(t(RNA_data))
  RNA_data <- RNA_data[-1,]
  # add SAMPLE_ID column to perform the alignment with clinical data
  RNA_data$SAMPLE_ID <- rownames(RNA_data)  
  return (RNA_data)
}
### Read & preprocess RNA expression data of Case samples
RNA_cases_data <- load_rna_data()
# Check for NA values
any(is.na(RNA_cases_data))
count <- which(colSums(is.na(RNA_cases_data))>0)
#count
# Remove columns that contain NA values
RNA_cases_data <- RNA_cases_data[ , colSums(is.na(RNA_cases_data))==0]
any(is.na(RNA_cases_data))

# During transposing numeric type values turned to character type
# Convert character type to numeric!
t1<-sapply(RNA_cases_data, as.numeric)
any(is.na(t1))
t11 <- t1[ , colSums(is.na(t1))==0]
any(is.na(t11))
t2 <- as.data.frame(t11)
t2$SAMPLE_ID <- rownames(RNA_cases_data)  # for alignment
# check whether NA exist
any(is.na(t2))

### Align clinical data with RNA_cases data
RNA_Sample_Ids <- t2$SAMPLE_ID
Clinical_Sample_Ids <- clinical_data$SAMPLE_ID
# get the list of smaple_ids common in RNA and clinical data
common_sample_ids<-intersect(Clinical_Sample_Ids, RNA_Sample_Ids)

# get clinical data only for common sample ids 
aligned_clinical_data <- clinical_data[clinical_data$SAMPLE_ID %in% common_sample_ids, ]
# get RNA data only for common sample ids
aligned_RNA_cases_data <- t2[RNA_cases_data$SAMPLE_ID %in% common_sample_ids, ]


# Check whether the alignment performed well 
RNA_Sample_Ids <- aligned_clinical_data$SAMPLE_ID
Clinical_Sample_Ids <- aligned_RNA_cases_data$SAMPLE_ID
# get those samples ids exist in RNA data and do not exist in clinical data
setdiff(RNA_Sample_Ids, Clinical_Sample_Ids) 
# get those samples ids that exist in clinical data and do not exists in RNA data
setdiff(Clinical_Sample_Ids, RNA_Sample_Ids) 

# Re-Confirm that there are no NA values
any(is.na(aligned_clinical_data$PFS_MONTHS))
any(is.na(aligned_clinical_data$PFS_STATUS))
any(is.na(aligned_RNA_cases_data))

### Overall Survival Analysis
# Create a survival object consisting of times & censoring
os_surv_obj <- Surv(time = aligned_clinical_data$OS_MONTHS, 
                     event = aligned_clinical_data$OS_STATUS)

# fit univariate survival curve (Kaplan Meier analysis)
fit <- survfit(os_surv_obj ~ 1, data = aligned_clinical_data)
fit
# the median survival time is : 49.3 months
survfit(os_surv_obj ~ 1, data = aligned_clinical_data) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Average survival (95% CI)**"
  )

# visualize the Kaplan Meier curve
ggsurvplot(fit, data = aligned_clinical_data, xlab = "Month", ylab = "Overall survival")

ggsurvplot(fit,
           data = aligned_clinical_data,
           risk.table = T,
           xlab = "Month", ylab = "Overall survival")

#to see the numbers at risk in a table below the x-axis.
survfit2(os_surv_obj ~ 1, data = aligned_clinical_data) %>% 
  ggsurvfit() +
  labs(
    x = "Month",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  add_risktable()

# Estimating 3-years survival
summary(survfit(os_surv_obj ~ 1, data = aligned_clinical_data), times = 72)
# the 3-year probability of survival in this study is 38.3%.

survfit(os_surv_obj ~ 1, data = aligned_clinical_data) %>% 
  tbl_survfit(
    times = 72,
    label_header = "**3-year survival (95% CI)**"
  )

# RADIATION THERAPY
# fit univariate survival curve (Kaplan Meier analysis) based on RADIATION THERAPY
fit1 <- survfit(os_surv_obj ~ RADIATION_THERAPY, data = aligned_clinical_data)
fit1
#summary(fit1)
ggsurvplot(fit1, data = aligned_clinical_data,
           pval = T,
           risk.table = T,
           xlab = "Month", ylab = "Overall survival")


# compute the probability of surviving 3 years grouped by RADIATION THERAPY
survfit(os_surv_obj ~ RADIATION_THERAPY, data = aligned_clinical_data) %>% 
  tbl_survfit(
    times = 72,
    label_header = "**3-year survival (95% CI)**"
  )

# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(os_surv_obj ~ RADIATION_THERAPY, 
                   data = aligned_clinical_data)
summary(fit.coxph)
ggforest(fit.coxph, data = aligned_clinical_data)
fit.coxph %>% 
  tbl_regression(exp = TRUE)

# compute the probability of surviving 3 years grouped by RADIATION THERAPY
survfit(os_surv_obj ~ RADIATION_THERAPY, data = aligned_clinical_data) %>% 
  tbl_survfit(
    times = 72,
    label_header = "**3-year survival (95% CI)**"
  )


## Progression Free Survival analysis

# create a survival object consisting of times & censoring
pfs_surv_obj <- Surv(time = aligned_clinical_data$PFS_MONTHS, 
                     event = aligned_clinical_data$PFS_STATUS)

# fit univariate survival curve (Kaplan Meier analysis)
fit <- survfit(pfs_surv_obj ~ 1, data = aligned_clinical_data)
# visualize the Kaplan Meier curve
ggsurvplot(fit, data = aligned_clinical_data, 
           risk.table = T,
           xlab = "Month", ylab = "Progressoin free survival")

#to see the numbers at risk in a table below the x-axis.
survfit2(pfs_surv_obj ~ 1, data = aligned_clinical_data) %>% 
  ggsurvfit() +
  labs(
    x = "Month",
    y = "Progression free survival"
  ) + 
  add_confidence_interval() +
  add_risktable()

# compute the median progression free time
survfit(pfs_surv_obj ~ 1, data = aligned_clinical_data) %>% 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Average PFS (95% CI)**"
  )
# compute the 1-year probability of PFS
survfit(pfs_surv_obj ~ 1, data = aligned_clinical_data) %>% 
  tbl_survfit(
    times = 12,
    label_header = "**1-year PFS (95% CI)**"
  )

# fit univariate survival curve (Kaplan Meier analysis) based on RADIATION THERAPY
fit1 <- survfit(pfs_surv_obj ~ RADIATION_THERAPY, data = aligned_clinical_data)
#summary(fit1)
ggsurvplot(fit1, data = aligned_clinical_data, 
           pval = T,
           risk.table = T,
           xlab = "Month", ylab = "Progression Free survival")


# compute the probability of surviving 3 years grouped by RADIATION THERAPY
survfit(pfs_surv_obj ~ RADIATION_THERAPY, data = aligned_clinical_data) %>% 
  tbl_survfit(
    times = 12,
    label_header = "**1-year PFS (95% CI)**"
  )

# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(pfs_surv_obj ~ RADIATION_THERAPY, 
                   data = aligned_clinical_data)
summary(fit.coxph)
ggforest(fit.coxph, data = aligned_clinical_data)
fit.coxph %>% 
  tbl_regression(exp = TRUE)


# fit univariate survival curve (Kaplan Meier analysis) based on SEX
fit1 <- survfit(pfs_surv_obj ~ SEX, data = aligned_clinical_data)
ggsurvplot(fit1, data = aligned_clinical_data, 
           pval = T,
           risk.table = T,
           xlab = "Month", ylab = "Progression Free survival")
# No signification difference in PFS when stratified by Sex

# univariate Cox regression based on Sex
fit.coxph <- coxph(pfs_surv_obj ~ AGE, 
                   data = aligned_clinical_data)

fit.coxph %>% 
  tbl_regression(exp = TRUE)
# HR = 1.0 : no difference between groups

# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(pfs_surv_obj ~ SEX + RADIATION_THERAPY + AGE, 
                   data = aligned_clinical_data)
ggforest(fit.coxph, data = aligned_clinical_data)
fit.coxph %>% 
  tbl_regression(exp = TRUE)


# Overall Survival and Genes
### Does the genes, reported by the literature as the most influencing genes 
# in lung adenocarcinoma, also influence the OS ?  
# list of genes reported by literature
citated_genes <- c("EGFR", "KRAS", "BRAF", "NF1", "ALK", "ROS1", 
                   "PTEN", "PIK3CA", "CDKN2A", "MYC", "`NKX2-1`")
# dynamically built of formula
# Refer to https://stackoverflow.com/questions/64251165/r-how-to-put-all-the-column-names-of-a-dataframe-into-a-formula
fS <- os_surv_obj ~ . 
fs = reformulate(citated_genes, fS[[2]])
fs
# fit the multivariate Cox regression
fit.coxph <- coxph(fs, data = aligned_RNA_cases_data)
summary(fit.coxph)
fit.coxph %>% 
  tbl_regression(exp = TRUE)
ggforest(fit.coxph, data = aligned_RNA_cases_data)

# Is the dysregulation of genes associated with overall survival?
# For EGFR and KRAS genes check whether their dysregulation is associated with
# the overall survival
# Divide patients into low, medium and high expressing groups

# Define the cut-off for critical z-score values to be -1.96 and +1.96 
# (when using 95% CI)
cut_off = 1.96
# Divide patients into low (z-score<-1.96), medium (-1.96<=z-score<=1.96) and 
# high (z-score<=1.96) expressing groups

# EGFR
# denote which cases have high,low or medium expression
aligned_RNA_cases_data$EGFR_strata <- ifelse(aligned_RNA_cases_data$EGFR < cut_off, "LOW", 
                                             ifelse(aligned_RNA_cases_data$EGFR > cut_off, "HIGH", "MEDIUM"))
# fitting survival curve 
fit <- survfit(os_surv_obj ~ EGFR_strata, data = aligned_RNA_cases_data)
fit
ggsurvplot(fit,
           data = aligned_RNA_cases_data,
           pval = T,
           risk.table = T,
           xlab = "Month", ylab = "Overall survival")
# KRAS
# denote which cases have high,low or medium expression
aligned_RNA_cases_data$KRAS_strata <- ifelse(aligned_RNA_cases_data$KRAS < cut_off, "LOW", 
                                             ifelse(aligned_RNA_cases_data$KRAS > cut_off, "HIGH", "MEDIUM"))
# fitting survival curve 
fit <- survfit(os_surv_obj ~ KRAS_strata, data = aligned_RNA_cases_data)
fit
ggsurvplot(fit,
           data = aligned_RNA_cases_data,
           pval = T,
           risk.table = T,
           xlab = "Month", ylab = "Overall survival")
# drop
drops <- c("EGFR_strata","KRAS_strata")
aligned_RNA_cases_data <- aligned_RNA_cases_data[ , !(names(aligned_RNA_cases_data) %in% drops)]

# Do the mutated oncogenes reported in cBioportal for the TCGA-LUAD dataset 
# affect the overall survival rate of lung adenocarcinoma? 
# Download from cbioportal site the Mutated Genes table into Mutated_Genes.txt file
# Load the mutated genes
filename = paste0(folder,"\\Mutated_Genes.txt")
mutated_genes <- read.delim(filename)
colnames(mutated_genes) <- c("Gene", "Mul.Sig", "No.Mutations", 
                             "No.Patients.with.Mutation","No.Profiled.Samples",
                             "Freq", "Is.Cancer.Gene")
# remove % from Freq column
mutated_genes$Freq<-gsub("%","",as.character(mutated_genes$Freq))
# change Freq type from character to numeric
mutated_genes <- transform(mutated_genes, Freq = as.numeric(Freq))
## Feature selection
# get those genes that are tagged as 'Cancer' genes and their 
# frequency is more than or equal of 95% of all mutated genes
# and sort them by their frequency value (descending)
cancer_genes<-mutated_genes %>% 
  arrange(desc(Freq))  %>%
  filter(Is.Cancer.Gene == 'Yes') %>%
  mutate(c = cume_dist(Freq)) %>%
  filter(c >= 0.95) %>% 
  select(Gene, Freq)
#cancer_genes
head(cancer_genes, 5)
# create a subset of the RNA data to contain only oncogenes
a<-cancer_genes$Gene
subsetRNA_cases_data<-aligned_RNA_cases_data %>%
  select(any_of(a))

# fit multivariate model (COX proportional hazard) based on oncogenes
fit.coxph <- coxph(os_surv_obj ~ ., data = subsetRNA_cases_data)
ggforest(fit.coxph, data = subsetRNA_cases_data)
fit.coxph %>% 
  tbl_regression(exp = TRUE)

# Fit Regularized Cox regression (by using glmnet R package) to find out the 
# genes that statistically influence the survival rate. 
# Prepare the data for glmnet function
# Remove SAMPLE_ID column (lst column)
my_x<-aligned_RNA_cases_data[1:(length(aligned_RNA_cases_data)-1)]
# convert RNA data from dataframe to matrix
my_x <- as.matrix(my_x)
# get the  time
my_time <- aligned_clinical_data$OS_MONTHS
# get the status
my_status <- aligned_clinical_data$OS_STATUS
# create a matrix with time and status
my_y <- as.matrix(data.frame(time=my_time,status=my_status))

# check for NA values
any(is.na(aligned_RNA_cases_data))
any(is.na(my_time))
any(is.na(my_status))
# create a dataframe that contain time, status and RNA data
my_data <- data.frame(time=my_time,status=my_status,my_x)

# Create a survival object (for overall survival)
my_z <- Surv(my_time, my_status)

# Fit the glmnet - Regularized Cox regression
# It takes some time !
my_fit <- glmnet(my_x, my_z, family = "cox")
print(my_fit)
plot(my_fit)
title("Regularized Cox Regression on Overall Survival", line = 2.5, cex.main = 1)

# Analysis
my_cfs = coef(my_fit)
summary(my_cfs)
# extract the coefficients at a certain value of λ 
# λ = 0.102400 returns the 16 most significant genes
my_cfs = coef(my_fit, s = 0.102400) # gets 16 genes
head(my_cfs)
# get the relevant genes
my_meaning_coefs = rownames(my_cfs)[my_cfs[,1]!= 0]
my_meaning_coefs
# and their values
my_meaning_vals = my_cfs[my_cfs[,1]!=0,]
my_meaning_vals

# dynamically built of formula
fS <- my_z ~ . 
fs = reformulate(my_meaning_coefs, fS[[2]])
fs

# use the built formula to fit the Cox regression
fit.coxph <- coxph(fs , data = my_data)
ggforest(fit.coxph, data = my_data)
fit.coxph %>% 
  tbl_regression(exp = TRUE)

# Compute K-fold cross-validation (CV) for the Cox model
set.seed(1)
# It takes some time (~2-3 min) !
#cv.fit <- cv.glmnet(my_x, my_z, family = "cox", alpha = .5, lambda = NULL)
cv.fit <- cv.glmnet(my_x, my_z, family = "cox", type.measure = "C")

# view the optimal λ value and a cross validated error plot
plot(cv.fit)
title("Cross Validation of Regularized Cox Regression on Overall Survival", line = 2.5, cex.main = 1)

# optimal value of λ which minimizes the cross-validation error
cv.fit$lambda.min
cv.fit$lambda.1se
# check which covariates the model chose as active, 
# and see the coefficients of those covariates
est.coef = coef(cv.fit, s = cv.fit$lambda.min)
active.k = which(est.coef != 0)
active.k
meaning_coefs = rownames(est.coef)[active.k]
meaning_coefs
# and their values
active.k.vals = est.coef[active.k]
active.k.vals

# dynamically build of formula
fS <- my_z ~ . 
fs = reformulate(meaning_coefs, fS[[2]])
fs

# use the built formula to git the Cox regression
fit.coxph <- coxph(fs, data = my_data)
ggforest(fit.coxph, data = my_data)
fit.coxph %>% 
  tbl_regression(exp = TRUE)

