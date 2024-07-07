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

# SEQ_ASSAY_ID and SEQ_PIPELINE_ID
seq_grp <- assayInfo %>%
  dplyr::group_by(SEQ_ASSAY_ID,SEQ_PIPELINE_ID) %>%
  dplyr::summarise(number_of_panels=n())
outF <- paste0(outDir,"/institution_seq_platoform_IDs.csv")
write.table(seq_grp,outF,row.names=F,quote=F,sep=",")

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
  #dplyr::group_by(ONCOTREE_CODE) %>%
  #group_by(!!sym("ONCOTREE_CODE")) %>%
  group_by(as.symbol("ONCOTREE_CODE")) %>%
  dplyr::summarise(n_patients_total=dplyr::n_distinct(Tumor_Sample_Barcode),
                   CANCER_TYPE=paste0(unique(CANCER_TYPE),collapse=";")) %>%
  dplyr::arrange(desc(n_patients_total))

otCodeOncokBSummary <- oncokbOnly %>%
  dplyr::group_by(ONCOTREE_CODE,ONCOKB_HIGHEST_LEVEL_SUMMARY) %>%
  dplyr::summarize(n.variants.matched=dplyr::n_distinct(balderVariantID),
                   n.patients.matched=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients.matched))

otCTypeSummary <- cType.sample.counts %>%
  dplyr::left_join(otCodeOncokBSummary,by="ONCOTREE_CODE") %>%
  dplyr::mutate(perc=100*(n.patients.matched/n_patients_total))
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
  dplyr::summarise(n_patients_total=dplyr::n_distinct(Tumor_Sample_Barcode),
                   CANCER_TYPE=paste0(unique(CANCER_TYPE),collapse=";")) %>%
  dplyr::arrange(desc(n_patients_total))

# compile clinical info by sequencing type
seqTypeClinicalInfo <- svCompiled[svCompiled$SEQ_ASSAY_ID_mod %in% assaySelect,] %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>% 
  dplyr::filter(dplyr::row_number() == 1) %>% # select a single representative entry per subject
  dplyr::ungroup() %>%
  dplyr::group_by(ONCOTREE_CODE,SEQ_ASSAY_ID_mod) %>% # SAMPLE_TYPE
  dplyr::summarise(n_patients_total=dplyr::n_distinct(Tumor_Sample_Barcode),
                   CANCER_TYPE=paste0(unique(CANCER_TYPE),collapse=";"),
                   n_Metastatic=sum(SAMPLE_TYPE == "Metastasis"),
                   n_Primary=sum(SAMPLE_TYPE == "Primary"),
                   met_to_primary_ratio=n_Metastatic/n_Primary,
                   n_LocalRecurrence=sum(SAMPLE_TYPE == "LocalRecurance"),
                   perc_LocalRecurrence=n_LocalRecurrence/n_patients_total,
                   meanTumorPurity=mean(TUMOR_PURITY,na.rm=T),
                   anyTumorPurityInfo=any(!is.na(TUMOR_PURITY)),
                   n_NaHeme = sum(SAMPLE_TYPE == "NA or Heme"),
                   perc_NaHeme=n_NaHeme/n_patients_total,
                   median_age=median(AGE_AT_SEQ_REPORT,na.rm=T),
                   n_dead=sum(DEAD == "TRUE" | DEAD == "True",na.rm=T),
                   perc_dead=n_dead/n_patients_total,
                   median_IntSeqToDeath=median(INT_Seq_or_Biopsy_To_Death, na.rm=T),
                   n_StageI=sum(STAGE=="I" | STAGE=="1" | STAGE=="1a" | STAGE=="Ia" | STAGE=="IA" |  STAGE=="Ib" | STAGE=="IB" | STAGE=="1b", na.rm=T),
                   perc_StageI=n_StageI/n_patients_total,
                   n_StageII=sum(STAGE=="II" | STAGE=="2" | STAGE=="2a" | STAGE=="IIa" | STAGE=="IIA" | STAGE=="IIb" | STAGE=="IIB" | STAGE=="2b", na.rm=T),
                   perc_StageII=n_StageII/n_patients_total,
                   n_StageIII=sum(STAGE=="III" | STAGE=="3" | STAGE=="3a" | STAGE=="IIIA" | STAGE=="IIIB" | STAGE=="3B", na.rm=T),
                   perc_StageIII=n_StageIII/n_patients_total,
                   n_StageIV=sum(STAGE=="IV" | STAGE=="M1" | STAGE=="M2" | STAGE=="M3", na.rm=T),
                   perc_StageIV=n_StageIV/n_patients_total) %>%
  dplyr::arrange(desc(n_patients_total))
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
  dplyr::filter(n_patients_total>50) %>%
  dplyr::mutate(perc=100*(n.patients.matched/n_patients_total)) %>%
  dplyr::left_join(seqTypeClinicalInfo[,!colnames(seqTypeClinicalInfo) %in% c("n_patients_total","CANCER_TYPE")],
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

##################################
### proportion of null entries ###
##################################

# what proportion of each column is null
perSubjVar <- svCompiled[svCompiled$SEQ_ASSAY_ID_mod %in% assaySelect,] %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>% 
  dplyr::filter(dplyr::row_number() == 1) # select a single representative entry per subject
#  dplyr::ungroup() %>%

nullCnts <- data.frame()  
i<-1
colList <- colnames(perSubjVar) 
for (x in colList[-length(colList)]) {
  for (y in unique(perSubjVar$SEQ_ASSAY_ID_mod)) {
    iSeq <- perSubjVar$SEQ_ASSAY_ID_mod == y
    selectCol <- perSubjVar[iSeq,x]
    iNull <- is.na(selectCol) | (selectCol=="")
    
    nullCnts[i,"column"] <- x
    nullCnts[i,"SEQ_ASSAY_ID_mod"] <- y
    nullCnts[i,"proportion_null"] <- mean(iNull)
    i<-i+1
  }
}

nullCntsWide <- nullCnts %>%
  pivot_wider(names_from = column,values_from=proportion_null)
outF <- paste0(outDir,"/null_entry_counts_v1.csv")
write.table(na_proportions_raw,outF,row.names=F,quote=F,sep=",")

na_proportions_raw <- perSubjVar %>%
  dplyr::group_by(SEQ_ASSAY_ID_mod) %>% 
  summarise(across(tidyselect::everything(), ~ mean(is.na(.)))) #%>%
  #summarise(dplyr::across(clinical_cols, ~ mean(is.na(.)))) %>%
  #pivot_longer(cols = everything(), names_to = "column", values_to = "proportion")
outF <- paste0(outDir,"/null_entry_counts_v2.csv")
write.table(na_proportions_raw,outF,row.names=F,quote=F,sep=",")



# na_proportions_summary <- seqOtCTypeSummary %>%
#   dplyr::group_by(SEQ_ASSAY_ID_mod) %>% 
#   #summarise(across(tidyselect::everything(), ~ mean(is.na(.)))) %>%
#   summarise(dplyr::across(clinical_cols, ~ mean(is.na(.)))) %>%
#   pivot_longer(cols = everything(), names_to = "column", values_to = "proportion")
# 


nullClinInfo <-  perSubjVar %>%
  dplyr::group_by(SEQ_ASSAY_ID_mod) %>% # SAMPLE_TYPE
  dplyr::summarise(n_patients_total=dplyr::n_distinct(Tumor_Sample_Barcode),
                   CANCER_TYPE=paste0(unique(CANCER_TYPE),collapse=";"),
                   n_Metastatic=sum(SAMPLE_TYPE == "Metastasis"),
                   n_Primary=sum(SAMPLE_TYPE == "Primary"),
                   met_to_primary_ratio=n_Metastatic/n_Primary,
                   n_LocalRecurrence=sum(SAMPLE_TYPE == "LocalRecurance"),
                   perc_LocalRecurrence=n_LocalRecurrence/n_patients_total,
                   meanTumorPurity=mean(TUMOR_PURITY,na.rm=T),
                   anyTumorPurityInfo=any(!is.na(TUMOR_PURITY)),
                   n_NaHeme = sum(SAMPLE_TYPE == "NA or Heme"),
                   perc_NaHeme=n_NaHeme/n_patients_total,
                   median_age=median(AGE_AT_SEQ_REPORT,na.rm=T),
                   n_dead=sum(DEAD == "TRUE" | DEAD == "True",na.rm=T),
                   perc_dead=n_dead/n_patients_total,
                   median_IntSeqToDeath=median(INT_Seq_or_Biopsy_To_Death, na.rm=T),
                   n_StageI=sum(STAGE=="I" | STAGE=="1" | STAGE=="1a" | STAGE=="Ia" | STAGE=="IA" |  STAGE=="Ib" | STAGE=="IB" | STAGE=="1b", na.rm=T),
                   perc_StageI=n_StageI/n_patients_total,
                   n_StageII=sum(STAGE=="II" | STAGE=="2" | STAGE=="2a" | STAGE=="IIa" | STAGE=="IIA" | STAGE=="IIb" | STAGE=="IIB" | STAGE=="2b", na.rm=T),
                   perc_StageII=n_StageII/n_patients_total,
                   n_StageIII=sum(STAGE=="III" | STAGE=="3" | STAGE=="3a" | STAGE=="IIIA" | STAGE=="IIIB" | STAGE=="3B", na.rm=T),
                   perc_StageIII=n_StageIII/n_patients_total,
                   n_StageIV=sum(STAGE=="IV" | STAGE=="M1" | STAGE=="M2" | STAGE=="M3", na.rm=T),
                   perc_StageIV=n_StageIV/n_patients_total) %>%
  dplyr::arrange(desc(n_patients_total))
outF <- paste0(outDir,"/cancer_type_assay_clinical_info.csv")
write.table(seqTypeClinicalInfo,outF,row.names=F,quote=F,sep=",")



############################
### exploratory plotting ###
############################

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

clinical_cols <- c("median_age",
                   "met_to_primary_ratio",
                   "meanTumorPurity",
                   "anyTumorPurityInfo",
                   "perc_LocalRecurrence",
                   "perc_NaHeme",
                   "perc_dead",
                   "perc_StageI",
                   "perc_StageII",
                   "perc_StageIII",
                   "perc_StageIV",
                   "median_IntSeqToDeath")
technical_cols <- c("is_paired_end",
                    "library_selection",
                    "platform",
                    "number_of_genes",
                    "calling_strategy",
                    "preservation_technique")

modelCols <- c("perc",
               "ONCOTREE_CODE",
               clinical_cols,
               technical_cols)

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

###########################################
### chi-sq test for categorical columns ###
###########################################

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

############################################
### corr structure of clinical variables ###
############################################
library(corrplot)

clin_train <- train_data[,clinical_cols]

# Compute Pearson correlation matrix
pearson_corr <- cor(clin_train, 
                    use = "pairwise.complete.obs",
                    method = "pearson")

# Visualize the Pearson correlation matrix
outF <-  paste0(outDir,"/clinical_factor_pairwise_pearson_corr.pdf")
pdf(outF, width = 10, height = 10)
corrplot(pearson_corr, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45, 
         title = "Pearson Correlation Matrix", mar = c(0,0,1,0))
dev.off()

# Compute spearman correlation matrix
spearman_corr <- cor(clin_train, 
                    #use = "pairwise.complete.obs",
                    method = "spearman")

# Visualize the Pearson correlation matrix
outF <-  paste0(outDir,"/clinical_factor_pairwise_spearman_corr.pdf")
pdf(outF, width = 10, height = 10)
corrplot(spearman_corr, method = "color", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45, 
         title = "Spearman Correlation Matrix", mar = c(0,0,1,0))
dev.off()

### proportion of null values
# na_proportions <- train_data %>%
#   summarise_all(~ mean(is.na(.))) %>%
#   gather(key = "column", value = "proportion")

na_proportions <- train_data[,clinical_cols] %>%
  #summarise(across(tidyselect::everything(), ~ mean(is.na(.)))) %>%
  summarise(dplyr::across(clinical_cols, ~ mean(is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "column", values_to = "proportion")

#####################################################
### perform per-patient modeling of level1 status ### 
#####################################################

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

##################
### assay info ###
##################

inF <- "../../data/AACR_Project_GENIE/Release_15p0_public/genomic_information.txt"
genomicInfo <- read.csv(inF,sep="\t")

# SEQ_ASSAY_ID and SEQ_PIPELINE_ID
seq_grp <- assayInfo %>%
  dplyr::group_by(SEQ_ASSAY_ID,SEQ_PIPELINE_ID) %>%
  dplyr::summarise(number_of_panels=n())
outF <- paste0(outDir,"/institution_seq_platoform_IDs.csv")
write.table(seq_grp,outF,row.names=F,quote=F,sep=",")

table(svCompiled$SEQ_ASSAY_ID %in% genomicInfo$SEQ_ASSAY_ID,svCompiled$SourceStudy)

#
length(unique(genomicInfo$Hugo_Symbol))

# rank genes according to which genees were measured most frequently
geneRank <- genomicInfo %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(n_assays=dplyr::n_distinct(SEQ_ASSAY_ID)) %>%
  dplyr::arrange(desc(n_assays))

geneRank$Hugo_Symbol <- factor(geneRank$Hugo_Symbol,levels=geneRank$Hugo_Symbol)
outF <-  paste0(outDir,"/gene_counts_by_seq_assay_ID_top.pdf")
ggplot(geneRank[1:60,], aes(x = Hugo_Symbol, y = n_assays)) +
  geom_bar(stat = "identity") +
  labs(title = "Count of assays measuring gene of interest (n=113 total)",
       x = "gene",
       y = "number of assays") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 7, width = 10)

geneRank$Hugo_Symbol <- factor(geneRank$Hugo_Symbol,levels=geneRank$Hugo_Symbol)
outF <-  paste0(outDir,"/gene_counts_by_seq_assay_ID_all.pdf")
ggplot(geneRank, aes(x = Hugo_Symbol, y = n_assays)) +
  geom_bar(stat = "identity") +
  labs(title = "Count of assays measuring gene of interest (n=113 total)",
       x = "gene",
       y = "number of assays") +
  theme_minimal()+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.x = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 7, width = 10)

## what proportion actionable mutation appear in top ranking genes?
iLevel1 <- grepl("1",oncokbOnly$ONCOKB_HIGHEST_LEVEL_SUMMARY)
lvl1_compiled <- oncokbOnly[iLevel1,]

lvl1_gene_cnts <- lvl1_compiled %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::summarise(n_patients=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n_patients))

rank_lvl1 <- geneRank %>%
  dplyr::left_join(lvl1_gene_cnts,by="Hugo_Symbol") 
rank_lvl1$n_patients[is.na(rank_lvl1$n_patients)] <- 0  
rank_lvl1$n_patients_cumsum <- cumsum(rank_lvl1$n_patients)
rank_lvl1$cumsum_perc_patients <- 100 * (rank_lvl1$n_patients_cumsum / length(unique(svCompiled$Tumor_Sample_Barcode)))
rank_lvl1$gene_rank <- rownames(rank_lvl1)

rank_lvl1$Hugo_Symbol <- factor(rank_lvl1$Hugo_Symbol,levels=geneRank$Hugo_Symbol)
outF <-  paste0(outDir,"/cumulative_sum_patients_with_level1_gene_mutation.pdf")
ggplot(rank_lvl1, aes(x = Hugo_Symbol, y = cumsum_perc_patients)) +
  geom_point()+
  #geom_bar(stat = "identity") +
  labs(title = "Cumulative sum of patients with a level 1 mutaiton \n genes ranked by most commonly measured in panels",
       x = "gene rank",
       y = "percent of patients") +
  theme_minimal()+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.x = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 7, width = 10)

rLvl1_long <- rank_lvl1[,c("Hugo_Symbol","n_assays","cumsum_perc_patients")] %>%
  tidyr::pivot_longer(!Hugo_Symbol, names_to = "type", values_to = "value") 
  
rLvl1_long$Hugo_Symbol <- factor(rLvl1_long$Hugo_Symbol,levels=geneRank$Hugo_Symbol)
outF <-  paste0(outDir,"/cumulative_sum_patients_with_level1_gene_mutation_and_prop.pdf")
ggplot(rLvl1_long, aes(x = Hugo_Symbol, y = value)) +
  geom_point()+
  facet_grid(type ~ ., scales="free_y")+
  labs(title = "Cumulative sum of patients with a level 1 mutaiton \n genes ranked by most commonly measured in panels",
       x = "gene rank",
       y = "percent of patients") +
  theme_minimal()+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.x = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 7, width = 10)

outF <-  paste0(outDir,"/cumulative_sum_patients_with_level1_gene_mutation_and_prop_top.pdf")
ggplot(rLvl1_long[1:180,], aes(x = Hugo_Symbol, y = value)) +
  geom_point()+
  facet_grid(type ~ ., scales="free_y")+
  labs(title = "Cumulative sum of patients with a level 1 mutaiton \n genes ranked by most commonly measured in panels",
       x = "gene rank",
       y = "percent of patients") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 7, width = 12)

# what are examples of genes with high frequency, but low panel rank? 
iHighRank <- as.numeric(rank_lvl1$gene_rank) >= 30
iHighPatient <- rank_lvl1$n_patients >= 300
unexpected_hits <- rank_lvl1[iHighRank & iHighPatient,]

# print info on unexpected hits
for (x in unexpected_hits$Hugo_Symbol) {
  print(x)
  iGeneMatch = lvl1_compiled$Hugo_Symbol==x
  print("actionability type:")
  print(table(lvl1_compiled[iGeneMatch,"ONCOKB_HIGHEST_LEVEL_SUMMARY"]))
  print("disease level counts:")
  print(table(lvl1_compiled[iGeneMatch,"level_1_disease"]))
  print(table(lvl1_compiled[iGeneMatch,"level_3_disease"]))
  print("study counts:")
  print(table(lvl1_compiled[iGeneMatch,"SourceStudy"]))
  print("panel types:")
  print(table(lvl1_compiled[iGeneMatch,"SEQ_ASSAY_ID_mod"]))
  print("--------")
}


### pairwise intersection of genes among panels
runPairwiseCalc <- 0
if (runPairwiseCalc==1){
  pairwiseIntersection <- data.frame()
  panels <- unique(genomicInfo$SEQ_ASSAY_ID)
  i <- 1
  for (ix in seq(length(panels))) {
    x <- panels[ix]
    xg <- genomicInfo[genomicInfo$SEQ_ASSAY_ID==x,]
    xgenes <- unique(xg$Hugo_Symbol)
    for (iy in seq(length(panels))) {
      y <- panels[iy]
      yg <- genomicInfo[genomicInfo$SEQ_ASSAY_ID==y,]
      ygenes <- unique(yg$Hugo_Symbol)
      if (ix > iy) {
        print("")
      }
      else {
        print(ix) 
        print(iy) 
        print("----")
        pairwiseIntersection[i,"panel_1"] <- x
        pairwiseIntersection[i,"panel_2"] <- y
        pairwiseIntersection[i,"panel_1_n_genes"] <- length(xgenes)
        pairwiseIntersection[i,"panel_2_n_genes"] <- length(ygenes)
        pairwiseIntersection[i,"panel_1_2_intersection"] <- length(intersect(xgenes,ygenes))
      }
      i <- i+1
    }
  }
  pairwiseIntersection$intersection_as_prop_of_panel1 <- pairwiseIntersection$panel_1_2_intersection / pairwiseIntersection$panel_1_n_genes
}

### genomic intersection plots
#geneList <- c("KRAS","TP53","ABL1","EGFR")
geneList <- c("FBXW7","PPP2R1A","CCNE1")
#geneList <- geneRank$Hugo_Symbol[1:16]

for (gene in geneList) {
  # hit info
  lvl1_comp_gene <- lvl1_compiled[lvl1_compiled$Hugo_Symbol==gene,]
  nPosTotal <- dim(lvl1_comp_gene)[1]
  lvl1PosCnts <- lvl1_comp_gene  %>%
    dplyr::group_by(Start_Position) %>%
    dplyr::summarise(n=n(),
                     geneLevelFreq=n/nPosTotal) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::rename(coordinate=Start_Position,
                  size=n) %>%
    dplyr::mutate(pos="actionable_count",
                  SEQ_ASSAY_ID="1_observed_events",
                  type="observed")
    
  # genomic info
  geneSubset <- genomicInfo[genomicInfo$Hugo_Symbol==gene,]
  geneSubset$region_ID <- seq(1,dim(geneSubset)[1])
  geneSubset$assay_ID <- factor(geneSubset$SEQ_ASSAY_ID,labels=seq(1,length(unique(geneSubset$SEQ_ASSAY_ID))))
  geneSubset$assay_ID <- as.numeric(geneSubset$assay_ID)
  geneSubset$region_length <- geneSubset$End_Position - geneSubset$Start_Position
  geneSubset$type <- "panel"
  geneSubset$size <- 1
  
  intervalPlotDf <- geneSubset[,c("Start_Position","End_Position","SEQ_ASSAY_ID","size","type")] %>%
    tidyr::pivot_longer(c("Start_Position","End_Position"), names_to = "pos", values_to = "coordinate")
  
  densityPlotDf <- rbind(intervalPlotDf,lvl1PosCnts[,colnames(intervalPlotDf)])
  
  outF <-  paste0(outDir,"/genomic_positions_of_panels_for_gene_",gene,".pdf")
  ggplot(intervalPlotDf, aes(x = coordinate, y = SEQ_ASSAY_ID, shape = pos)) +
    geom_point()+
    #facet_grid(type ~ ., scales="free_y")+
    labs(title = paste0(gene, "\nGenomic coordiantes tested by panel"),
         x = "genomic coordinate",
         y = "panel") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(outF,height = 18, width = 12)
  
  outF <-  paste0(outDir,"/genomic_positions_of_panels_for_gene_",gene,"_with_counts.pdf")
  ggplot(densityPlotDf, aes(x = coordinate, y = SEQ_ASSAY_ID, shape = pos, size=size,color=type)) +
    geom_point(alpha=.2,size=1,color="black")+
    geom_point(alpha=.6)+
    #facet_grid(type ~ ., scales="free_y")+
    labs(title = paste0(gene, "\nGenomic coordiantes tested by panel"),
         x = "genomic coordinate",
         y = "panel") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(outF,height = 18, width = 12)
  
  # create position-level evaluations
  lvl1PosCnts$temp_key <- 1
  geneSubset$temp_key <- 1
  
  # Use merge to create the Cartesian product
  posTargetMerged <- merge(lvl1PosCnts[,!colnames(lvl1PosCnts) %in% c("SEQ_ASSAY_ID")], 
                  geneSubset[,!colnames(geneSubset) %in% c("type","size","")], by = "temp_key")
  
}

################################
### detailed tumor type info ###
################################

c_type_count <- perSubjVar %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(number_of_subjects=n()) %>%
  dplyr::arrange(desc(number_of_subjects))

c_type_detailed <- perSubjVar %>%
  dplyr::group_by(CANCER_TYPE,CANCER_TYPE_DETAILED) %>%
  dplyr::summarise(number_of_subjects=n()) %>%
  dplyr::arrange(desc(number_of_subjects))
outF <- paste0(outDir,"/cancer_type_detailed_counts.txt")
write.table(c_type_detailed,outF,row.names=F,quote=F,sep="\t")

iTop <- c_type_detailed$CANCER_TYPE %in% c_type_count$CANCER_TYPE[1:10]
c_type_subset <- c_type_detailed[iTop,]

c_type_subset$CANCER_TYPE_DETAILED <- factor(c_type_subset$CANCER_TYPE_DETAILED,levels=c_type_subset$CANCER_TYPE_DETAILED)
outF <-  paste0(outDir,"/cancer_type_detailed_counts.pdf")
ggplot(c_type_detailed[1:20,], aes(x = CANCER_TYPE_DETAILED, y = number_of_subjects)) +
  geom_bar(stat = "identity") +
  labs(title = "Sample counts for each detailed cancer type categorizaiton (most common cancers)",
       x = "gene",
       y = "number of assays") +
  facet_wrap(CANCER_TYPE ~ ., scales="free_x") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 14, width = 10)

