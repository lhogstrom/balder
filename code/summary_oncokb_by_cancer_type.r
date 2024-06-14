library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)
library(ggplot2)
source("R/snv_indel_annotation.R")
source("R/utils.R")

### connect to result DB and get variants
outDir <- "../../output/oncokb_cancer_type_summary_20240528"
bDir <- "../../data/processed/balderResultsDb"
dbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v20240412.sqlite")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), dbName)

svCompiled <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM patientObservedVariantTableWithSampleInfoOncokb')
aa.genome.representative <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM representativeClinicalAnnotatedPatientVariants')
aa.genome.representative.oncokb <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM representativeClinicalAnnotatedPatientVariantsWithOncoKB')

#################################
### summary of oncokb results ###
#################################

### find oncokb matched variants
onckbAnnotated <- svCompiled %>%
  dplyr::filter(ONCOKB_IS_SENS_DX_PX_OR_RESISTANCE == T)

table(onckbAnnotated$balderVariantID %in% aa.genome.representative$balderVariantID)
table(aa.genome.representative$balderVariantID %in% onckbAnnotated$balderVariantID)

## which variants are onckb annotated but not annotated via MIA/CiVIC? 
iOnckbOnlyGrp <- onckbAnnotated$balderVariantID %in% aa.genome.representative$balderVariantID
okBID <- onckbAnnotated[!iOnckbOnlyGrp,"balderVariantID"]
oncokbOnly <-svCompiled[svCompiled$balderVariantID %in% okBID,]  
  
###########################
### oncokb level counts ###
###########################

inF <- "../../data/AACR_Project_GENIE/Release_15p0_public/assay_information.txt"
assayInfo <- read.csv(inF,sep="\t")

# variant counts per subject
svAssaySummary <- svCompiled %>%
  filter(!is.na(SEQ_ASSAY_ID)) %>%
  dplyr::group_by(SEQ_ASSAY_ID,Tumor_Sample_Barcode) %>%
  dplyr::summarize(n.variants=dplyr::n_distinct(balderVariantID),
                   n.patients=dplyr::n_distinct(Tumor_Sample_Barcode),
                   n.oncokb.level1=sum(ONCOKB_HIGHEST_LEVEL_SUMMARY == "LEVEL_1"),
                   has.oncokb.level1=any(ONCOKB_HIGHEST_LEVEL_SUMMARY == "LEVEL_1")) %>%
  dplyr::left_join(assayInfo,by="SEQ_ASSAY_ID") %>%
  dplyr::arrange(desc(n.variants))
svAssaySummary[is.na(svAssaySummary$has.oncokb.level1),"has.oncokb.level1"] <- F



# nullSeq <- is.na(svAssaySummary$SEQ_ASSAY_ID)
# outF <-  paste0(outDir,"/variant_count_by_assay_type.pdf")
# ggplot(svAssaySummary[!nullSeq,])+
#   #geom_point(alpha=.2)+
#   geom_density_2d(aes(x=number_of_genes,y=n.variants,color=library_selection))+
#   theme_bw()+
#   #facet_grid(SAMPLE_TYPE~.,scale="free_y",)+
#   xlab("Number of genes measured")+
#   ylab("number of variants")+
#   ggtitle(paste0("Per patient variant count by assay type"))+
#   scale_fill_brewer(palette="Set1",drop=FALSE)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   theme(plot.title = element_text(hjust = 0.5))
# ggsave(outF,height = 8,width = 7)

###########################
### oncokb level counts ###
###########################

### create an incidence table by cancer type
cType.sample.counts <- svCompiled %>%
  dplyr::group_by(ONCOTREE_CODE) %>%
  dplyr::summarise(n.patients.total=dplyr::n_distinct(Tumor_Sample_Barcode),
                   CANCER_TYPE=paste0(unique(CANCER_TYPE),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.total))

otCodeOncokBSummary <- oncokbOnly %>%
  dplyr::group_by(ONCOTREE_CODE,ONCOKB_HIGHEST_LEVEL_SUMMARY) %>%
  dplyr::summarize(n.variants.matched=dplyr::n_distinct(balderVariantID),
                   n.patients.matched=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients.matched))

otCTypeSummary <- cType.sample.counts %>%
  dplyr::left_join(otCodeOncokBSummary,by="ONCOTREE_CODE") %>%
  dplyr::mutate(perc=100*(n.patients.matched/n.patients.total))
otCTypeSummary$ONCOTREE_CODE <- factor(otCTypeSummary$ONCOTREE_CODE, levels=unique(otCTypeSummary$ONCOTREE_CODE))
outF <- paste0(outDir,"/oncokb_match_level_summary_by_ot_code.csv")
write.table(otCTypeSummary,outF,row.names=F,quote=F,sep=",")

highest_c_counts <- unique(otCTypeSummary$ONCOTREE_CODE)[1:10]
iHighest <- otCTypeSummary$ONCOTREE_CODE %in% as.character(highest_c_counts)

outF <-  paste0(outDir,"/cancer_type_oncokb_evidence_level.pdf")
ggplot(otCTypeSummary[iHighest,],aes(x=ONCOTREE_CODE,y=perc,fill=ONCOKB_HIGHEST_LEVEL_SUMMARY))+
  geom_bar(stat = "identity", position = "stack",alpha=.6) +
  theme_bw()+
  #facet_grid(SAMPLE_TYPE~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("perc of patients")+
  ggtitle(paste0("Oncokb evidence level by cancer type"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)

otTissueOncokBSummary <- oncokbOnly %>%
  dplyr::group_by(ot_code_level_1,ONCOKB_HIGHEST_LEVEL_SUMMARY) %>%
  dplyr::summarize(n.variants.restricted.matched=dplyr::n_distinct(balderVariantID),
                   n.patients.restricted.matched=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients.restricted.matched))

#############################################
### oncokb level counts  - primary vs met ###
#############################################

iPrimOrMet <- (oncokbOnly$SAMPLE_TYPE=="Primary") | (oncokbOnly$SAMPLE_TYPE=="Metastasis")
primMetOncokBSummary <- oncokbOnly[iPrimOrMet,] %>%
  dplyr::group_by(ONCOTREE_CODE,ONCOKB_HIGHEST_LEVEL_SUMMARY,SAMPLE_TYPE) %>%
  dplyr::summarize(n.variants.matched=dplyr::n_distinct(balderVariantID),
                   n.patients.matched=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients.matched))
primMetOncokBSummary$ONCOTREE_CODE <- factor(primMetOncokBSummary$ONCOTREE_CODE, levels=unique(otCTypeSummary$ONCOTREE_CODE))

highest_c_counts <- unique(otCTypeSummary$ONCOTREE_CODE)[1:20]
iHighest <- primMetOncokBSummary$ONCOTREE_CODE %in% as.character(highest_c_counts)

iLevel <- grepl("1",primMetOncokBSummary$ONCOKB_HIGHEST_LEVEL_SUMMARY)

outF <-  paste0(outDir,"/cancer_type_prim_met_oncokb_evidence_level.pdf")
ggplot(primMetOncokBSummary[iHighest & iLevel,],aes(x=ONCOTREE_CODE,y=n.patients.matched,fill=ONCOKB_HIGHEST_LEVEL_SUMMARY))+
  geom_bar(stat = "identity", position = "stack",alpha=.6) +
  theme_bw()+
  facet_grid(SAMPLE_TYPE~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("perc of patients")+
  ggtitle(paste0("Oncokb evidence level by cancer type"))+
  #scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)

#######################################
### oncokb level counts  - assay ID ###
#######################################

# Which data sources see the largest differences in actionability rates? Is this because of cohort differences or technology differences? Which specific mutations explain the key differences?

oncokbOnly$SEQ_ASSAY_ID_mod <- oncokbOnly$SEQ_ASSAY_ID
oncokbOnly[is.na(oncokbOnly$SEQ_ASSAY_ID),"SEQ_ASSAY_ID_mod"] <- oncokbOnly[is.na(oncokbOnly$SEQ_ASSAY_ID),"SourceStudy"]
svCompiled$SEQ_ASSAY_ID_mod <- svCompiled$SEQ_ASSAY_ID
svCompiled[is.na(svCompiled$SEQ_ASSAY_ID),"SEQ_ASSAY_ID_mod"] <- svCompiled[is.na(svCompiled$SEQ_ASSAY_ID),"SourceStudy"]

assayCntRank <- svCompiled %>%
  dplyr::group_by(SEQ_ASSAY_ID_mod) %>%
  dplyr::summarise(n=dplyr::n_distinct(Tumor_Sample_Barcode)) %>% 
  dplyr::arrange(desc(n)) %>%
  dplyr::left_join(assayInfo,by=c("SEQ_ASSAY_ID_mod"="SEQ_ASSAY_ID")) 
outF <- paste0(outDir,"/institution_seq_platoform.csv")
write.table(assayCntRank,outF,row.names=F,quote=F,sep=",")

# select assays with >N patients
assaySelect <- assayCntRank[assayCntRank$n>5000,"SEQ_ASSAY_ID_mod"][[1]]

seqTypeCntAssay <- svCompiled[svCompiled$SEQ_ASSAY_ID_mod %in% assaySelect,] %>%
  dplyr::group_by(ONCOTREE_CODE,SEQ_ASSAY_ID_mod) %>% # SAMPLE_TYPE
  dplyr::summarise(n.patients.total=dplyr::n_distinct(Tumor_Sample_Barcode),
                   CANCER_TYPE=paste0(unique(CANCER_TYPE),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.total))

# compile clinical info by sequencing type
seqTypeClinicalInfo <- svCompiled[svCompiled$SEQ_ASSAY_ID_mod %in% assaySelect,] %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>% 
  dplyr::filter(dplyr::row_number() == 1) %>% # select a single representative entry per subject
  dplyr::ungroup() %>%
  dplyr::group_by(ONCOTREE_CODE,SEQ_ASSAY_ID_mod) %>% # SAMPLE_TYPE
  dplyr::summarise(n.patients.total=dplyr::n_distinct(Tumor_Sample_Barcode),
                   CANCER_TYPE=paste0(unique(CANCER_TYPE),collapse=";"),
                   n.Metastatic=sum(SAMPLE_TYPE == "Metastasis"),
                   n.Primary=sum(SAMPLE_TYPE == "Primary"),
                   n.LocalRecurrence=sum(SAMPLE_TYPE == "LocalRecurance"),
                   meanTumorPurity=mean(TUMOR_PURITY,na.rm=T),
                   n.NaHeme = sum(SAMPLE_TYPE == "NA or Heme"),
                   median.age=median(AGE_AT_SEQ_REPORT,na.rm=T),
                   n.dead=sum(DEAD == "TRUE" | DEAD == "True",na.rm=T),
                   median.IntSeqToDeath=median(INT_Seq_or_Biopsy_To_Death, na.rm=T),
                   n.StageI=sum(STAGE=="I" | STAGE=="1" | STAGE=="1a" | STAGE=="Ia" | STAGE=="Ib" | STAGE=="1b",na.rm=T),
                   n.StageII=sum(STAGE=="II" | STAGE=="2" | STAGE=="2a" | STAGE=="IIa" | STAGE=="IIb" | STAGE=="2b",na.rm=T)) %>%
  dplyr::arrange(desc(n.patients.total))
outF <- paste0(outDir,"/cancer_type_assay_clinical_info.csv")
write.table(seqTypeClinicalInfo,outF,row.names=F,quote=F,sep=",")


iPrimOrMet <- (oncokbOnly$SAMPLE_TYPE=="Primary") | (oncokbOnly$SAMPLE_TYPE=="Metastasis")
seqOncokBSummary <- oncokbOnly[iPrimOrMet,] %>%
  dplyr::group_by(ONCOTREE_CODE,SEQ_ASSAY_ID_mod,ONCOKB_HIGHEST_LEVEL_SUMMARY) %>%
  dplyr::summarize(n.variants.matched=dplyr::n_distinct(balderVariantID),
                   n.patients.matched=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients.matched))
primMetOncokBSummary$ONCOTREE_CODE <- factor(primMetOncokBSummary$ONCOTREE_CODE, levels=unique(otCTypeSummary$ONCOTREE_CODE))

seqOtCTypeSummary <- seqTypeCntAssay %>%
  dplyr::left_join(seqOncokBSummary,by=c("ONCOTREE_CODE","SEQ_ASSAY_ID_mod")) %>%
  dplyr::filter(n.patients.total>50) %>%
  dplyr::mutate(perc=100*(n.patients.matched/n.patients.total)) %>%
  dplyr::left_join(seqTypeClinicalInfo[,!colnames(seqTypeClinicalInfo) %in% c("n.patients.total","CANCER_TYPE")],
                   by=c("ONCOTREE_CODE","SEQ_ASSAY_ID_mod")) %>%
  dplyr::left_join(assayInfo,by=c("SEQ_ASSAY_ID_mod"="SEQ_ASSAY_ID" ))
seqOtCTypeSummary$ONCOTREE_CODE <- factor(seqOtCTypeSummary$ONCOTREE_CODE, levels=unique(otCTypeSummary$ONCOTREE_CODE))
outF <- paste0(outDir,"/oncokb_match_level_summary_by_seq_type.csv")
write.table(seqOtCTypeSummary,outF,row.names=F,quote=F,sep=",")

highest_c_counts <- unique(otCTypeSummary$ONCOTREE_CODE)[1:5]
iHighest <- seqOtCTypeSummary$ONCOTREE_CODE %in% as.character(highest_c_counts)
iLevel1 <- grepl("1",seqOtCTypeSummary$ONCOKB_HIGHEST_LEVEL_SUMMARY)
iLevel2 <- grepl("2",seqOtCTypeSummary$ONCOKB_HIGHEST_LEVEL_SUMMARY)
seqOtCTypeSummary$level_1 <- iLevel1

outF <-  paste0(outDir,"/cancer_type_Seq_Type_oncokb_evidence_level.pdf")
ggplot(seqOtCTypeSummary[iHighest & (iLevel1),],aes(x=ONCOTREE_CODE,y=perc,fill=ONCOKB_HIGHEST_LEVEL_SUMMARY))+
  geom_bar(stat = "identity", position = "stack",alpha=.6) +
  theme_bw()+
  facet_grid(SEQ_ASSAY_ID_mod~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("perc of patients")+
  ylim(0,60)+
  ggtitle(paste0("Oncokb evidence level by cancer type"))+
  #scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 30, width = 10)

highest_c_counts <- unique(otCTypeSummary$ONCOTREE_CODE)[1:40]
iHighest <- seqOtCTypeSummary$ONCOTREE_CODE %in% as.character(highest_c_counts)
outF <-  paste0(outDir,"/cancer_type_Seq_Type_oncokb_level1_jitter.pdf")
ggplot(seqOtCTypeSummary[iHighest & (iLevel1),],aes(x=ONCOTREE_CODE,y=perc,color=ONCOKB_HIGHEST_LEVEL_SUMMARY))+
  #geom_bar(stat = "identity", position = "stack",alpha=.6) +
  geom_jitter(width=.1,alpha=.8)+
  theme_bw()+
  facet_grid(level_1~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("perc of patients")+
  ylim(0,85)+
  ggtitle(paste0("Oncokb evidence level by cancer type"))+
  scale_color_brewer(palette="Set2",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8, width = 18)

outF <-  paste0(outDir,"/cancer_type_Seq_Type_oncokb_jitter_by_level.pdf")
ggplot(seqOtCTypeSummary[iHighest & (iLevel1 | iLevel2),],aes(x=ONCOTREE_CODE,y=perc,color=ONCOKB_HIGHEST_LEVEL_SUMMARY))+
  #geom_bar(stat = "identity", position = "stack",alpha=.6) +
  geom_jitter(width=.1,alpha=.8)+
  theme_bw()+
  facet_grid(level_1~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("perc of patients")+
  ylim(0,85)+
  ggtitle(paste0("Oncokb evidence level by cancer type"))+
  scale_color_brewer(palette="Set2",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8, width = 18)

#############################################
### level1 and level 2 hits by assay type ###
#############################################

iLevel1 <- grepl("1",seqOtCTypeSummary$ONCOKB_HIGHEST_LEVEL_SUMMARY)
iLevel2 <- grepl("2",seqOtCTypeSummary$ONCOKB_HIGHEST_LEVEL_SUMMARY)

level12_hits_per_subj <- oncokbOnly[iLevel1 | iLevel2,] %>% 
  group_by(Tumor_Sample_Barcode,SEQ_ASSAY_ID_mod) %>% 
  dplyr::summarise(n=n()) %>%
  dplyr::group_by(n,SEQ_ASSAY_ID_mod) %>%
  dplyr::summarise(nPerSub=n()) %>%
  dplyr::left_join(assayInfo,by=c("SEQ_ASSAY_ID_mod"="SEQ_ASSAY_ID" )) 

###########################
### feature engineering ###
###########################

# # Handle missing values
# preProcess_missingdata <- preProcess(data, method = 'knnImpute')
# data_imputed <- predict(preProcess_missingdata, newdata = data)
# 
# # Scale and center numeric variables
# preProcess_scale <- preProcess(data_imputed[, sapply(data_imputed, is.numeric)], method = c("center", "scale"))
# data_imputed[, sapply(data_imputed, is.numeric)] <- predict(preProcess_scale, newdata = data_imputed[, sapply(data_imputed, is.numeric)])
# 
# # Convert categorical variables to dummy variables
# dummies_model <- dummyVars(ResponseVar ~ ., data = data_imputed)
# data_dummies <- predict(dummies_model, newdata = data_imputed)
# 
# # Combine the response variable back with the transformed predictors
# data_final <- cbind(data_dummies, ResponseVar = data_imputed$ResponseVar)
# 
# # Convert to a data frame for further processing
# data_final <- as.data.frame(data_final)
# 
# # Display the final processed dataset
# print(data_final)

###############################################
### modeling level1 hits by sequencing site ###
###############################################

# model which sequence facility factors have the biggest impact on proportion of level1 variants per cancer type

library(rpart)
library(rpart.plot)
#library(mlflow)

highest_c_counts <- unique(otCTypeSummary$ONCOTREE_CODE)[1:15]
iHighest <- seqOtCTypeSummary$ONCOTREE_CODE %in% as.character(highest_c_counts)
iLevel1 <- grepl("LEVEL_1",seqOtCTypeSummary$ONCOKB_HIGHEST_LEVEL_SUMMARY)

modelCols <- c("perc",
               "ONCOTREE_CODE",
               "is_paired_end",
               "library_selection",
               "platform",
               "number_of_genes",
               "calling_strategy",
               "preservation_technique")

data <- seqOtCTypeSummary[iHighest & iLevel1,] %>%
  #dplyr::left_join(assayInfo,by=c("SEQ_ASSAY_ID_mod"="SEQ_ASSAY_ID" )) %>%
  dplyr::select(all_of(modelCols)) %>%
  dplyr::group_by(ONCOTREE_CODE) %>%
  dplyr::arrange(desc(perc)) %>% 
  dplyr::mutate(actionable_perc_rank=dplyr::row_number())
  
# Convert categorical variables to factors (if not already factors)
data$ONCOTREE_CODE <- factor(data$ONCOTREE_CODE, levels=unique(data$ONCOTREE_CODE))
data$is_paired_end <- as.factor(data$is_paired_end)
data$library_selection <- as.factor(data$library_selection)
data$platform <- as.factor(data$platform)
data$calling_strategy <- as.factor(data$calling_strategy)
data$preservation_technique <- as.factor(data$preservation_technique)

# Split the data into training and testing sets
set.seed(42)
train_index <- caret::createDataPartition(data$perc, p = 0.8, list = FALSE)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

# Define the formula for the model
#formula <- perc ~ ONCOTREE_CODE + is_paired_end + library_selection + preservation_technique
formula <- actionable_perc_rank ~ is_paired_end + library_selection + preservation_technique

# Train the Decision Tree Regressor
tree_model <- rpart(formula, data = train_data, method = "anova")

# Print the model summary
print(summary(tree_model))

# Visualize the tree
rpart.plot(tree_model)

# Extract the rules (cut points) from the tree
tree_rules <- as.data.frame(tree_model$frame)
tree_rules_no_leaves <- tree_rules[tree_rules$var != "<leaf>", ]  # Exclude leaf nodes

# Print the rules
print(tree_rules_no_leaves)

# plot model
plot(tree_model)
text(tree_model, pretty = 0)

# Predict on the test data
predictions <- predict(tree_model, newdata = test_data)

# Evaluate the model performance
mse <- mean((predictions - test_data$perc)^2)
cat("Mean Squared Error (MSE):", mse, "\n")

#### 

# feature analysis 
library(corrplot)
library(ggplot2)
library(reshape2)

seqOtFull <- seqOtCTypeSummary[iHighest & iLevel1,] %>%
  dplyr::left_join(assayInfo,by=c("SEQ_ASSAY_ID_mod"="SEQ_ASSAY_ID" )) 


# 1. Correlation matrix and heatmap for numeric variables
# numeric_vars <- df[, sapply(df, is.numeric)]
# cor_matrix <- cor(numeric_vars, use = "complete.obs")
# corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)

# 2. Chi-squared test of independence for categorical variables
aInfo <- data.frame(lapply(assayInfo, as.factor))
cSelect <- colnames(aInfo)
cSelect <- cSelect[!cSelect %in% c("SEQ_ASSAY_ID")]
categorical_vars <- aInfo[, cSelect]
chi_sq_matrix <- matrix(NA, ncol = ncol(categorical_vars), nrow = ncol(categorical_vars))
colnames(chi_sq_matrix) <- colnames(categorical_vars)
rownames(chi_sq_matrix) <- colnames(categorical_vars)

for (i in 1:ncol(categorical_vars)) {
  for (j in 1:ncol(categorical_vars)) {
    if (i != j) {
      chi_sq_matrix[i, j] <- chisq.test(categorical_vars[, i], categorical_vars[, j])$p.value
    } else {
      chi_sq_matrix[i, j] <- NA
    }
  }
}

chi_sq_melt <- melt(chi_sq_matrix, na.rm = TRUE)
colnames(chi_sq_melt) <- c("Var1", "Var2", "value")
chi_sq_melt$log_value <- -log10(chi_sq_melt$value)

# Plot chi-squared test results as heatmap
outF <-  paste0(outDir,"/sequencing_platoform_categorical_variables_chi_sq.pdf")
ggplot(data = chi_sq_melt, aes(x = Var1, y = Var2, fill = log_value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey50", name = "-log10 p-value") +
  #scale_fill_gradient(low = "red", high = "white", na.value = "grey50", name = "p-value") + # , trans = "reverse"
  theme_minimal() +
  labs(title = "Chi-squared Test of Independence for Categorical Variables", x = "Variables", y = "Variables") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(outF,height = 10, width = 10)

#### perform per-patient modeling of level1 status ### 

# create representative variant per row per subject - choose level 1 variant if avaible 
perPatientLevel1 <- svCompiled %>%
  dplyr::left_join(assayInfo,by=c("SEQ_ASSAY_ID_mod"="SEQ_ASSAY_ID")) %>%
  dplyr::mutate(isOncoKbLevel1=ONCOKB_HIGHEST_LEVEL_SUMMARY=="LEVEL_1" & !is.na(ONCOKB_HIGHEST_LEVEL_SUMMARY)) %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::arrange(desc(isOncoKbLevel1)) %>%
  dplyr::filter(dplyr::row_number()==1) %>%
  dplyr::ungroup()
  
# make data subset selections
highest_c_counts <- unique(otCTypeSummary$ONCOTREE_CODE)[1:5]
iHighest <- perPatientLevel1$ONCOTREE_CODE %in% as.character(highest_c_counts)

modelCols <- c("isOncoKbLevel1",
               "ONCOTREE_CODE",
               "is_paired_end",
               "library_selection",
               "platform",
               "number_of_genes",
               "calling_strategy",
               "preservation_technique")


data <- perPatientLevel1[iHighest,] %>%
  dplyr::select(all_of(modelCols))
iNull <- rowSums(!is.na(data[,modelCols])) != length(modelCols)
data <- data[!iNull,]

# Convert categorical variables to factors (if not already factors)
data$isOncoKbLevel1 <- factor(data$isOncoKbLevel1)
data$ONCOTREE_CODE <- factor(data$ONCOTREE_CODE, levels=unique(data$ONCOTREE_CODE))
data$is_paired_end <- as.factor(data$is_paired_end)
data$library_selection <- as.factor(data$library_selection)
data$platform <- as.factor(data$platform)
data$calling_strategy <- as.factor(data$calling_strategy)
data$preservation_technique <- as.factor(data$preservation_technique)

# Split the data into training and testing sets
set.seed(42)
train_index <- caret::createDataPartition(data$isOncoKbLevel1, p = 0.3, list = FALSE)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

# Define the formula for the model
formula <- isOncoKbLevel1 ~ ONCOTREE_CODE + is_paired_end + library_selection + preservation_technique

# Train the Decision Tree Regressor
tree_model <- rpart(formula, data = train_data, method = "class")

# Print the model summary
print(summary(tree_model))

# Visualize the tree
rpart.plot(tree_model)

# Extract the rules (cut points) from the tree
tree_rules <- as.data.frame(tree_model$frame)
tree_rules_no_leaves <- tree_rules[tree_rules$var != "<leaf>", ]  # Exclude leaf nodes

# Print the rules
print(tree_rules_no_leaves)

# plot model
plot(tree_model)
text(tree_model, pretty = 0)

# Predict on the test data
predictions <- predict(tree_model, newdata = test_data)

# Evaluate the model performance
mse <- mean((predictions - test_data$perc)^2)
cat("Mean Squared Error (MSE):", mse, "\n")


### Apply GLM model

# Fit a GLM model
glm_model <- glm(isOncoKbLevel1 ~ ONCOTREE_CODE + is_paired_end + library_selection + preservation_technique, data = train_data, family = gaussian)

# Display the summary of the model
summary(glm_model)

# Extract coefficients and convert them to a data frame
coefficients_df <- broom::tidy(glm_model)

# Plot the coefficients
outF <-  paste0(outDir,"/glm_coefficients.pdf")
ggplot(coefficients_df, aes(x = term, y = estimate)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  labs(title = "GLM Coefficients",
       x = "Predictors",
       y = "Coefficient Estimate") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 7, width = 7)

### try z-scoring proportion data by cancer type (or rank-based)


