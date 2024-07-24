library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)
#library(ggrepel)
library(ggplot2)
#library(ggalluvial)

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 2) {
  stop("Two arguments must be supplied", call. = FALSE)
}

# Assign Arguments to Variables
baseDir <- args[1] # should be the relative path to "balder/code" where this file is located
dbName <- args[2]

utilFile1 <- paste0(baseDir,"/../code/R/snv_indel_annotation.R")
source(utilFile1)
utilFile2 <- paste0(baseDir,"/../code/R/utils.R")
source(utilFile1)

## load key data products
#svCompiled
#aa.genome.representative.oncokb


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

aa.genome.representative.oncokb <- rbind(aa.genome.representative,oncokbOnly)
RSQLite::dbWriteTable(harmonizedDb, "representativeClinicalAnnotatedPatientVariantsWithOncoKB", aa.genome.representative.oncokb, overwrite=T)
### Disconnect from SQL db ###
RSQLite::dbDisconnect(harmonizedDb)

# what proportion of variants are in oncokb database? 
table(svCompiled$ONCOKB_VARIANT_IN_ONCOKB)

# 2x2 of in oncokb and in exiting matching approach
table(svCompiled$ONCOKB_VARIANT_IN_ONCOKB, svCompiled$balderVariantID %in% aa.genome.representative$balderVariantID)

# how many variants are not matched via oncokb, but are matched by other annotation sources? 
table(aa.genome.representative.oncokb$ONCOKB_VARIANT_IN_ONCOKB)

# breakdown of oncokb "highest_effect". How does it match to "actionability.summary"
table(aa.genome.representative.oncokb$ONCOKB_HIGHEST_LEVEL,aa.genome.representative.oncokb$clinical.evidence.summary)
# repeat for cancer type match
ctMatchRep <- aa.genome.representative.oncokb[aa.genome.representative.oncokb$cancerTypeMatch==T,]
table(ctMatchRep$ONCOKB_HIGHEST_LEVEL,ctMatchRep$clinical.evidence.summary)

### concordance of oncokb and CIVIC/MOA results - representative
iCtypeMatch <- (aa.genome.representative.oncokb$cancerTypeMatch == T)
iActionable <- (aa.genome.representative.oncokb$ONCOKB_IS_SENS_DX_PX_OR_RESISTANCE==T) |
  (!is.na(aa.genome.representative.oncokb$actionability.summary) & iCtypeMatch)
okbVariantSummary <- aa.genome.representative.oncokb[iActionable,] %>%
  dplyr::group_by(Hugo_Symbol,HGVSp_Short,ONCOTREE_CODE) %>%
  dplyr::summarise(n=dplyr::n(),
                   sensitivity_entry_cnt_oncokb=sum(ONCOKB_HAS_SENSITIVE_ENTRY),
                   sensitivity_entry_cnt_CIVIC_MOA=sum(actionability.summary=="Predictive therapy actionability"),
                   sensitivity_cnt_match=sensitivity_entry_cnt_oncokb==sensitivity_entry_cnt_CIVIC_MOA,
                   resistance_entry_cnt_oncokb=sum(ONCOKB_HAS_RESISTANCE_ENTRY),
                   resistance_entry_cnt_CIVIC_MOA=sum(actionability.summary=="Therapy resistance"),
                   resistance_cnt_match=resistance_entry_cnt_oncokb==resistance_entry_cnt_CIVIC_MOA,
                   Px_entry_cnt_oncokb=sum(ONCOKB_HAS_PX_ENTRY),
                   Px_entry_cnt_CIVIC_MOA=sum(actionability.summary=="Prognostic"),
                   Px_cnt_match=Px_entry_cnt_oncokb==Px_entry_cnt_CIVIC_MOA,
                   DX_entry_cnt_oncokb=sum(ONCOKB_HAS_DX_ENTRY),
                   other_cnt_MOA_CIVIC=sum(actionability.summary=="Other")) %>%
  dplyr::arrange(desc(n))
outF <- paste0(outDir,"/oncokb_MOA_CIVIC_variant_summary_table_v2.txt")
write.table(okbVariantSummary,outF,row.names=F,quote=F,sep="\t")

### concordance of oncokb and CIVIC/MOA results - exhaustive
iCtypeMatch <- aa.genome.exahustive$cancerTypeMatch == T
iActionable <- (aa.genome.exahustive$ONCOKB_IS_SENS_DX_PX_OR_RESISTANCE==T) |
  (!is.na(aa.genome.exahustive$actionability.summary) & iCtypeMatch)
okbVariantExhaustiveSummary <- aa.genome.exahustive[iActionable,] %>%
  dplyr::group_by(Hugo_Symbol,HGVSp_Short,ONCOTREE_CODE) %>%
  dplyr::summarise(n=dplyr::n_distinct(Tumor_Sample_Barcode),
                   sensitivity_entry_cnt_oncokb=sum(ONCOKB_HAS_SENSITIVE_ENTRY),
                   sensitivity_entry_cnt_CIVIC_MOA=sum(actionability.summary=="Predictive therapy actionability"),
                   sensitivity_cnt_match=sensitivity_entry_cnt_oncokb==sensitivity_entry_cnt_CIVIC_MOA,
                   resistance_entry_cnt_oncokb=sum(ONCOKB_HAS_RESISTANCE_ENTRY),
                   resistance_entry_cnt_CIVIC_MOA=sum(actionability.summary=="Therapy resistance"),
                   resistance_cnt_match=resistance_entry_cnt_oncokb==resistance_entry_cnt_CIVIC_MOA,
                   Px_entry_cnt_oncokb=sum(ONCOKB_HAS_PX_ENTRY),
                   Px_entry_cnt_CIVIC_MOA=sum(actionability.summary=="Prognostic"),
                   Px_cnt_match=Px_entry_cnt_oncokb==Px_entry_cnt_CIVIC_MOA,
                   DX_entry_cnt_oncokb=sum(ONCOKB_HAS_DX_ENTRY),
                   other_cnt_MOA_CIVIC=sum(actionability.summary=="Other")) %>%
  dplyr::arrange(desc(n))
outF <- paste0(outDir,"/oncokb_MOA_CIVIC_variant_summary_table_exhaustive.txt")
write.table(okbVariantExhaustiveSummary,outF,row.names=F,quote=F,sep="\t")

## manual check
# iGene <- aa.genome.representative.oncokb$Hugo_Symbol == "TERT"
# iOT <- aa.genome.representative.oncokb$ONCOTREE_CODE == "BLCA"
# iHGVSP <- aa.genome.representative.oncokb$HGVSp_Short == ""
# #
# iGene <- aa.genome.representative.oncokb$Hugo_Symbol == "EGFR"
# iOT <- aa.genome.representative.oncokb$ONCOTREE_CODE == "LUAD"
# iHGVSP <- aa.genome.representative.oncokb$HGVSp_Short == "p.E746_A750del"
# 
# head(data.frame(aa.genome.representative.oncokb[iGene & iOT & iHGVSP,]))

### build summary of oncokb annotations alone
iDx <- grepl("LEVEL_Dx",aa.genome.representative.oncokb$ONCOKB_HIGHEST_LEVEL_SUMMARY)
dxAnnotated <- aa.genome.representative.oncokb[iDx,]
table(dxAnnotated$ONCOKB_HIGHEST_LEVEL_SUMMARY,dxAnnotated$actionability.summary)

#iPx <- grepl("LEVEL_Px",aa.genome.representative.oncokb$ONCOKB_HIGHEST_LEVEL_SUMMARY) | 
#  grepl("Prognostic",aa.genome.representative.oncokb$actionability.summary)
iPx <- grepl("LEVEL_Px",aa.genome.representative.oncokb$ONCOKB_HIGHEST_LEVEL_SUMMARY)
pxAnnotated <- aa.genome.representative.oncokb[iPx,]
table(pxAnnotated$ONCOKB_HIGHEST_LEVEL_SUMMARY,pxAnnotated$actionability.summary,exclude=NULL)

######################################
### Build primary reporting tables ###
######################################

### non-annotated, background/control variants - 
non.annotated.variants <- svCompiled %>%
  dplyr::filter(!balderVariantID %in% aa.genome.representative$balderVariantID) %>%
  dplyr::mutate(clinical_annotation_status="not annotated")
annotated.and.non.annotated <- rbind(aa.genome.representative,non.annotated.variants) # same as svCompiled, but with annotation columns

### per-patient summary - from exhaustive table
exhaustive.matched.patient.summary <- aa.genome.exahustive %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::summarize(n.variants.matched=dplyr::n_distinct(balderVariantID),
                   anyOtCodeCancerTypeMatch=TRUE %in% cancerTypeMatch,
                   exactOtCodeCancerTypeMatch=TRUE %in% cancerTypeMatchPrimary,
                   actionSummary=paste0(unique(actionability.summary),collapse=";"),
                   predTherapyAction="Predictive therapy actionability" %in% actionability.summary,
                   prognosticAction="Prognostic" %in% actionability.summary,
                   resistanceAction="Therapy resistance" %in% actionability.summary,
                   otherAction="Other" %in% actionability.summary,
                   clinEvidienceSummary=paste0(unique(clinical.evidence.summary),collapse=";"),
                   highConfidenceMatch="High Confidence" %in% clinical.evidence.summary,
                   otCode=paste0(unique(ONCOTREE_CODE),collapse=";"),
                   disease=paste0(unique(Disease),collapse=";"),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.variants.matched))

### per-patient summary - from representative table
# exhaustive.matched.patient.summary <- aa.genome.representative %>%
#   dplyr::group_by(Tumor_Sample_Barcode) %>%
#   dplyr::summarize(n.variants.matched=dplyr::n_distinct(balderVariantID),
#                    anyCancerTypeMatch=TRUE %in% cancerTypeMatch,
#                    actionSummary=paste0(unique(actionability.summary),collapse=";"),
#                    predTherapyAction="Predictive therapy actionability" %in% actionability.summary,
#                    prognosticAction="Prognostic" %in% actionability.summary,
#                    resistanceAction="Therapy resistance" %in% actionability.summary,
#                    otherAction="Other" %in% actionability.summary,
#                    clinEvidienceSummary=paste0(unique(clinical.evidence.summary),collapse=";"),
#                    highConfidenceMatch="High Confidence" %in% clinical.evidence.summary,
#                    otCode=paste0(unique(ONCOTREE_CODE),collapse=";"),
#                    disease=paste0(unique(Disease),collapse=";"),
#                    source=paste0(unique(source),collapse=";"),
#                    Drugs=paste0(unique(Drugs),collapse=";"),
#                    ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
#   dplyr::arrange(desc(n.variants.matched))


### summary table
# Combine the pairs into a matrix, then convert to a data frame
pairs_matrix <- rbind(c("number of starting subjects",length(unique(svCompiled$Tumor_Sample_Barcode))), 
                      c("number of starting variants",dim(svCompiled)[[1]]),
                      c("number of exhaustive annotation matches (one-to-many mapping)",dim(aa.genome.exahustive)[[1]]),
                      c("number of subjects with at least one match",length(unique(aa.genome.exahustive$Tumor_Sample_Barcode))),
                      c("number of unique variants with annotation match",print(dim(aa.genome.representative))[[1]]),
                      c("number of subjects with an exact OncoTree cancer type annotation match",sum(exhaustive.matched.patient.summary$exactOtCodeCancerTypeMatch==T)),
                      c("number of subjects with any OncoTree hierarchy cancer type annotation match",sum(exhaustive.matched.patient.summary$anyOtCodeCancerTypeMatch==T)))
reportTable <- data.frame(pairs_matrix)
write.csv(reportTable,"../../output/actionability_db_curration_20231220/actionability_exhaustive_representative_counts.csv")


### create an incidence table by cancer type
cType.sample.counts <- svCompiled %>%
  dplyr::group_by(ONCOTREE_CODE) %>%
  dplyr::summarise(n.patients.total=dplyr::n_distinct(Tumor_Sample_Barcode),
                   CANCER_TYPE=paste0(unique(CANCER_TYPE),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.total))
head(cType.sample.counts)

cType.unrestricted.match.summary <- aa.genome.exahustive %>%
  dplyr::group_by(ONCOTREE_CODE) %>%
  dplyr::summarize(n.variants.unrestricted.matched=dplyr::n_distinct(balderVariantID),
                   n.patients.unrestricted.matched=dplyr::n_distinct(Tumor_Sample_Barcode),
                   #n.patients.CancerTypeMatch=sum(cancerTypeMatch==TRUE & !is.na(cancerTypeMatch)),
                   actionSummary=paste0(unique(actionability.summary),collapse=";"),
                   predTherapyAction="Predictive therapy actionability" %in% actionability.summary,
                   prognosticAction="Prognostic" %in% actionability.summary,
                   resistanceAction="Therapy resistance" %in% actionability.summary,
                   otherAction="Other" %in% actionability.summary,
                   clinEvidienceSummary=paste0(unique(clinical.evidence.summary),collapse=";"),
                   highConfidenceMatch="High Confidence" %in% clinical.evidence.summary,
                   disease=paste0(unique(Disease),collapse=";"),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.unrestricted.matched))

cType.restricted.match.summary <- aa.genome.exahustive %>%
  dplyr::filter(cancerTypeMatch==TRUE & !is.na(cancerTypeMatch)) %>%
  dplyr::group_by(ONCOTREE_CODE) %>%
  dplyr::summarize(n.variants.restricted.matched=dplyr::n_distinct(balderVariantID),
                   n.patients.restricted.matched=dplyr::n_distinct(Tumor_Sample_Barcode),
                   #n.patients.CancerTypeMatch=sum(cancerTypeMatch==TRUE & !is.na(cancerTypeMatch)),
                   actionSummary=paste0(unique(actionability.summary),collapse=";"),
                   predTherapyAction="Predictive therapy actionability" %in% actionability.summary,
                   prognosticAction="Prognostic" %in% actionability.summary,
                   resistanceAction="Therapy resistance" %in% actionability.summary,
                   otherAction="Other" %in% actionability.summary,
                   clinEvidienceSummary=paste0(unique(clinical.evidence.summary),collapse=";"),
                   highConfidenceMatch="High Confidence" %in% clinical.evidence.summary,
                   disease=paste0(unique(Disease),collapse=";"),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.restricted.matched))

cTypeSummary <- cType.sample.counts %>%
  dplyr::left_join(cType.unrestricted.match.summary[,c("n.variants.unrestricted.matched",
                                                       "n.patients.unrestricted.matched",
                                                       "ONCOTREE_CODE")],by="ONCOTREE_CODE") %>%
  dplyr::left_join(cType.restricted.match.summary[,c("n.variants.restricted.matched",
                                                     "n.patients.restricted.matched",
                                                     "highConfidenceMatch",
                                                     "disease",
                                                     "ONCOTREE_CODE")],by="ONCOTREE_CODE") %>%
  dplyr::arrange(desc(n.patients.total)) %>%
  dplyr::mutate(percPatientsUnrestircted=100*round(n.patients.unrestricted.matched/n.patients.total,3),
                percPatientsRestircted=100*round(n.patients.restricted.matched/n.patients.total,3),
                patientsUnrestricted=paste0(n.patients.unrestricted.matched," (",percPatientsUnrestircted,"%)"),
                patientsRestricted=paste0(n.patients.restricted.matched," (",percPatientsRestircted,"%)")) %>%
  dplyr::select(CANCER_TYPE, ONCOTREE_CODE, n.patients.total, patientsUnrestricted, patientsRestricted,
                everything()) #re-order columns
write.csv(cTypeSummary,"../../output/actionability_db_curration_20231220/test_cancer_type_match_summary_patient_sorted.csv")



#Columns of summary table:
# - cancer type
# - number of patients mutated
# - av. number mutations per patient
# - proportion with cancer type match
# - prop with high-confidence match
# - prop with action type: therapy, prognostic, resistance, or other

studyCntsTbl <- svCompiled %>%
  dplyr::group_by(SourceStudy) %>%
  dplyr::summarise(nSubjects=dplyr::n_distinct(Tumor_Sample_Barcode))

##############################
### reporting and analysis ###
##############################

### Cancer Type match summary plots - all samples

cTypeBar <- cTypeSummary %>%
  dplyr::filter(n.patients.total > 1000) %>%
  dplyr::arrange(desc(percPatientsUnrestircted)) %>%
  dplyr::mutate(RestrictedCancerType = replace(percPatientsRestircted, is.na(percPatientsRestircted), 0),
                UnrestrictedCancerType = percPatientsUnrestircted-RestrictedCancerType) %>%
  dplyr::select(c("ONCOTREE_CODE","RestrictedCancerType","UnrestrictedCancerType")) %>%
  tidyr::pivot_longer(!ONCOTREE_CODE,names_to = "match_type",values_to = "value")
cTypeBar$ONCOTREE_CODE <- factor(cTypeBar$ONCOTREE_CODE,levels=unique(cTypeBar$ONCOTREE_CODE))

outF <-  paste0(outDir,"/cancer_type_clinical_annotation_match_summary.pdf")
ggplot(cTypeBar,aes(x=ONCOTREE_CODE,y=value,fill=match_type))+
  geom_bar(stat = "identity", position = "stack",alpha=.6) +
  theme_bw()+
  #scale_fill_manual(values = c("X1" = "blue", "X2" = "red")) + 
  #facet_grid(SourceStudy~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("Percent of patients")+
  ggtitle(paste0("Clinical annotation rates for the 30 most profiled cancer types"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)

### restricted - primary vs. met 

cType.primary.met <- aa.genome.representative %>%
  dplyr::filter(!is.na(cancerTypeMatch)) %>%
  dplyr::group_by(SAMPLE_TYPE,ONCOTREE_CODE,cancerTypeMatch) %>%
  dplyr::summarize(n.patients.matched=dplyr::n_distinct(Tumor_Sample_Barcode)) #%>%

cType.sample.counts.primary.met <- svCompiled %>%
  dplyr::group_by(SAMPLE_TYPE,ONCOTREE_CODE) %>%
  dplyr::summarise(n.patients=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients)) %>%
  data.frame()

cTypeSummaryPrimaryMet <- cType.sample.counts.primary.met %>%
  dplyr::left_join(cType.primary.met,by=c("ONCOTREE_CODE","SAMPLE_TYPE")) %>%
  dplyr::filter(SAMPLE_TYPE %in% c("Primary","Metastasis")) %>%
  dplyr::mutate(percPatientsMatched=100*round(n.patients.matched/n.patients,3),
                patientsRestricted=paste0(n.patients.matched," (",percPatientsMatched,"%)")) %>%
  data.frame() #%>%

iPrimary <- cType.sample.counts.primary.met$SAMPLE_TYPE=="Primary"
iHighCount <- cType.sample.counts.primary.met$n.patients > 500
TopPrimTissue <- cType.sample.counts.primary.met[iPrimary & iHighCount,"ONCOTREE_CODE"]

cTypeBarPrimMet <- cTypeSummaryPrimaryMet %>%
  dplyr::rename(patient_count=n.patients,
                match_count=n.patients.matched) %>%
  tidyr::pivot_longer(c("patient_count","match_count"),names_to = "count_type",values_to = "value") %>%
  dplyr::filter(ONCOTREE_CODE %in% TopPrimTissue)

codeLevels <- c(cType.sample.counts.primary.met[cType.sample.counts.primary.met$SAMPLE_TYPE=="Primary","ONCOTREE_CODE"])
codeLevelsMod <- codeLevels[codeLevels %in% c(cTypeBarPrimMet$ONCOTREE_CODE)]
cTypeBarPrimMet$ONCOTREE_CODE <- factor(cTypeBarPrimMet$ONCOTREE_CODE,
                                        levels=codeLevels)

outF <-  paste0(outDir,"/cancer_type_clinical_annotation_match_summary_primary_met.pdf")
ggplot(cTypeBarPrimMet,aes(x=ONCOTREE_CODE,y=value,color=count_type))+
  geom_point(alpha=.4)+
  #geom_bar(stat = "identity", position = "stack",alpha=.6) +
  theme_bw()+
  #scale_fill_manual(values = c("X1" = "blue", "X2" = "red")) +
  facet_grid(SAMPLE_TYPE~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("number of patients")+
  #ggtitle(paste0("Clinical annotation rates for the 30 most profiled cancer types"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)

barData <- cTypeSummaryPrimaryMet[cTypeSummaryPrimaryMet$ONCOTREE_CODE %in% TopPrimTissue,] %>%
  dplyr::arrange(desc(percPatientsMatched))
barData$ONCOTREE_CODE <- factor(barData$ONCOTREE_CODE,levels=unique(barData$ONCOTREE_CODE))

outF <-  paste0(outDir,"/cancer_type_clinical_annotation_match_summary_primary_met_perc.pdf")
ggplot(barData,aes(x=ONCOTREE_CODE,y=percPatientsMatched,fill=SAMPLE_TYPE))+
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  theme_bw()+
  facet_grid(cancerTypeMatch~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("percent of patients")+
  ylim(c(0,99))+
  #ggtitle(paste0("Clinical annotation rates for the 30 most profiled cancer types"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)

cTypeSummaryPrimaryMetWide <- cTypeSummaryPrimaryMet[,c("SAMPLE_TYPE",
                                                        "ONCOTREE_CODE",
                                                        "n.patients",
                                                        "cancerTypeMatch",
                                                        "n.patients.matched")] %>%
  dplyr::filter(!is.na(cancerTypeMatch)) %>%
  pivot_wider(names_from = cancerTypeMatch,values_from=n.patients.matched) %>%
  dplyr::rename(restrictedMatchCount=`TRUE`,
                UnrestrictedTmpCount=`FALSE`) %>%
  dplyr::mutate(restrictedMatchCount=replace(restrictedMatchCount, is.na(restrictedMatchCount), 0),
                UnrestrictedTmpCount=replace(UnrestrictedTmpCount, is.na(UnrestrictedTmpCount), 0),
                UnrestrictedMatchCount=UnrestrictedTmpCount+restrictedMatchCount,
                PercMatchRestricted=100*(restrictedMatchCount/n.patients),
                PercMatchUnestricted=100*(UnrestrictedMatchCount/n.patients))


### perform chi-squared test
chiSqTbl <- cTypeSummaryPrimaryMet[,c("SAMPLE_TYPE",
                                      "ONCOTREE_CODE",
                                      "n.patients",
                                      "cancerTypeMatch",
                                      "n.patients.matched")] %>%
  dplyr::filter(!is.na(cancerTypeMatch)) %>%
  pivot_wider(names_from = SAMPLE_TYPE,values_from=c(n.patients.matched,n.patients)) %>%
  dplyr::rename(n_matched_Primary=n.patients.matched_Primary,
                n_matched_Metastasis=n.patients.matched_Metastasis,
                n_patients_Primary=n.patients_Primary,
                n_patients_Metastasis=n.patients_Metastasis) %>%
  dplyr::mutate(n_matched_Primary=replace(n_matched_Primary, is.na(n_matched_Primary), 0),
                n_matched_Metastasis=replace(n_matched_Metastasis, is.na(n_matched_Metastasis), 0),
                n_patients_Primary=replace(n_patients_Primary, is.na(n_patients_Primary), 0),
                n_patients_Metastasis=replace(n_patients_Metastasis, is.na(n_patients_Metastasis), 0),
                #
                n_not_matched_Primary=n_patients_Primary-n_matched_Primary,
                n_not_matched_Metastasis=n_patients_Metastasis-n_matched_Metastasis) %>%
  dplyr::filter(ONCOTREE_CODE %in% TopPrimTissue)

# Function to perform chi-squared test on a single row
perform_chi_squared_test <- function(row) {
  # Construct the contingency table for the current row
  M <- as.table(rbind(c(row["n_matched_Primary"][[1]], row["n_not_matched_Primary"][[1]]), 
                      c(row["n_matched_Metastasis"][[1]], row["n_not_matched_Metastasis"][[1]])))
  dimnames(M) <- list(c("Primary", "Metastatic"),
                      c("Matched", "Not Matched"))
  
  # Perform chi-square test
  test_result <- chisq.test(M)
  
  # Return the p-value (or any other statistic of interest)
  return(test_result$p.value)
}

# Apply the function to each row of the dataframe and collect results
cntCols <- c("n_matched_Primary","n_not_matched_Primary","n_matched_Metastasis","n_not_matched_Metastasis")
chiSqTbl$p_value <- apply(data.frame(chiSqTbl[, cntCols]), 1, function(row) perform_chi_squared_test(as.list(row)))
chiSqTbl <- chiSqTbl %>%
  dplyr::arrange(p_value) %>%
  dplyr::mutate(significance_str = if_else(p_value < 0.05, "*", "")) %>%
  dplyr::left_join(cType.sample.counts,
                   by="ONCOTREE_CODE")

outF <- paste0(outDir,"/chi_sq_test_of_prim_met_match_rate.csv")
write.csv(chiSqTbl,outF)

### Plot chiSqTbl results





ddBar <- cTypeSummaryPrimaryMetWide %>%
  dplyr::rename(restricted_cancer_type_match=PercMatchRestricted,
                unrestricted_cancer_type_match=PercMatchUnestricted) %>%
  tidyr::pivot_longer(c("restricted_cancer_type_match","unrestricted_cancer_type_match"),names_to = "match_type",values_to = "value") %>%
  dplyr::filter(ONCOTREE_CODE %in% TopPrimTissue) #%>%
ddBar$ONCOTREE_CODE <- factor(ddBar$ONCOTREE_CODE, levels=unique(chiSqTbl$ONCOTREE_CODE))


outF <-  paste0(outDir,"/cancer_type_clinical_annotation_match_summary_primary_met_perc2.pdf")
ggplot(ddBar,aes(x=ONCOTREE_CODE,y=value,fill=SAMPLE_TYPE))+
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  theme_bw()+
  facet_grid(match_type~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("percent of patients")+
  ylim(c(0,99))+
  #ggtitle(paste0("Clinical annotation rates for the 30 most profiled cancer types"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)


# restricted cancer types only
ddBar2 <- ddBar %>%
  dplyr::filter(match_type=="restricted_cancer_type_match",
                value>0) %>%
  dplyr::left_join(chiSqTbl[,c("ONCOTREE_CODE","significance_str")],
                   by="ONCOTREE_CODE")
ddBar2$ONCOTREE_CODE <- factor(ddBar2$ONCOTREE_CODE, levels=unique(chiSqTbl$ONCOTREE_CODE))
ddBar2[ddBar2$SAMPLE_TYPE=="Primary","significance_str"] <- ""

outF <-  paste0(outDir,"/cancer_type_clinical_annotation_match_summary_primary_met_perc3.pdf")
ggplot(ddBar2,aes(x=ONCOTREE_CODE,y=value,fill=SAMPLE_TYPE))+
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  theme_bw()+
  geom_text(aes(label = significance_str, y = value + 5))+
  #facet_grid(match_type~.,scale="free_y",)+
  xlab("Oncotree Code")+
  ylab("percent of patients")+
  ylim(c(0,99))+
  #ggtitle(paste0("Restricted clinical annotation"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 6,width = 6)


### look at resistance markers in primary vs. met

cType.primary.met.resistance <- aa.genome.exahustive %>%
  dplyr::filter(actionability.summary=="Therapy resistance",
                cancerTypeMatch == TRUE,
                !is.na(cancerTypeMatch),
                clinical.evidence.summary=="High Confidence") %>%
  dplyr::group_by(SAMPLE_TYPE,ONCOTREE_CODE) %>%
  dplyr::summarize(n.patients.matched=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients.matched))

# cType.sample.counts.primary.met <- svCompiled %>%
#   dplyr::group_by(SAMPLE_TYPE,ONCOTREE_CODE) %>%
#   dplyr::summarise(n.patients=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
#   dplyr::arrange(desc(n.patients)) %>%
#   data.frame()

highCntsSamples <- cType.sample.counts[cType.sample.counts$n.patients.total>100,]

cTypeSummaryPrimaryMetResistance <- cType.primary.met.resistance %>%
  dplyr::filter(ONCOTREE_CODE %in% highCntsSamples$ONCOTREE_CODE) %>%
  dplyr::left_join(cType.sample.counts.primary.met,by=c("ONCOTREE_CODE","SAMPLE_TYPE")) %>%
  dplyr::left_join(ot_code_full[,c("ot_code","ot_name")],by=c("ONCOTREE_CODE"="ot_code")) %>%
  dplyr::filter(SAMPLE_TYPE %in% c("Primary","Metastasis")) %>%
  dplyr::mutate(percPatientsMatched=100*round(n.patients.matched/n.patients,3),
                patientsRestricted=paste0(n.patients.matched," (",percPatientsMatched,"%)")) %>%
  dplyr::arrange(desc(percPatientsMatched)) %>%
  data.frame() #%>%
cTypeSummaryPrimaryMetResistance$ot_name <- factor(cTypeSummaryPrimaryMetResistance$ot_name,
                                                   levels=unique(cTypeSummaryPrimaryMetResistance$ot_name))
#cTypeSummaryPrimaryMetResistance$ONCOTREE_CODE <- factor(cTypeSummaryPrimaryMetResistance$ONCOTREE_CODE,
#                                                         levels=unique(cTypeSummaryPrimaryMetResistance$ONCOTREE_CODE))


outF <-  paste0(outDir,"/cancer_type_resistance_primary_met_perc.pdf")
ggplot(cTypeSummaryPrimaryMetResistance,aes(y=ot_name,x=percPatientsMatched,fill=SAMPLE_TYPE))+
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  theme_bw()+
  #geom_text(aes(label = significance_str, y = value + 5))+
  #facet_grid(match_type~.,scale="free_y",)+
  xlab("Oncotree Code")+
  xlab("percent of patients")+
  xlim(c(0,60))+
  #ggtitle(paste0("Restricted clinical annotation"))+
  scale_fill_brewer(palette="Set2",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 8)


### resistance markers primary vs. met

AAChange.primary.met.resistance <- aa.genome.exahustive %>%
  dplyr::filter(actionability.summary=="Therapy resistance",
                cancerTypeMatch == TRUE,
                !is.na(cancerTypeMatch),
                clinical.evidence.summary=="High Confidence") %>%
  dplyr::group_by(SAMPLE_TYPE,Hugo_Symbol,AAChange) %>%
  dplyr::summarize(n.patients.matched=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients.matched))

sampleTypeCnt <- aa.genome.exahustive %>%
  dplyr::group_by(SAMPLE_TYPE) %>%
  dplyr::summarize(n.patients=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n.patients))

AAChangePrimaryMetResistance <- AAChange.primary.met.resistance %>%
  #dplyr::filter(ONCOTREE_CODE %in% highCntsSamples$ONCOTREE_CODE) %>%
  dplyr::filter(SAMPLE_TYPE %in% c("Primary","Metastasis")) %>%
  dplyr::left_join(sampleTypeCnt,by=c("SAMPLE_TYPE")) %>%
  dplyr::mutate(aaStr=paste0(Hugo_Symbol,"-",AAChange),
                percPatientsMatched=100*round(n.patients.matched/n.patients,3),
                patientsRestricted=paste0(n.patients.matched," (",percPatientsMatched,"%)")) %>%
  dplyr::arrange(desc(percPatientsMatched)) %>%
  data.frame() #%>%
AAChangePrimaryMetResistance$aaStr <- factor(AAChangePrimaryMetResistance$aaStr,
                                             levels=unique(AAChangePrimaryMetResistance$aaStr))

outF <-  paste0(outDir,"/resistance_primary_met_perc_AA_change.pdf")
ggplot(AAChangePrimaryMetResistance,aes(y=aaStr,x=percPatientsMatched,fill=SAMPLE_TYPE))+
  geom_bar(stat = "identity", position = "dodge",alpha=.6) +
  theme_bw()+
  #geom_text(aes(label = significance_str, y = value + 5))+
  xlab("Resistance marker")+
  xlab("percent of patients")+
  #ggtitle(paste0("Restricted clinical annotation"))+
  scale_fill_brewer(palette="Set4",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 6,width = 4)

#################################
### Secondary reporting table ###
#################################

# Computing actionable variant counts for different match approaches:
# 1) Unrestricted where no cancer type matching was imposed
# 2) Restricted with raw cancer type string matching

###
### join with observed variant table with clinically annotated - by genome position
genome.matched.unrestricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome_annotation",
                                 "Start_Position"="pos",
                                 "Reference_Allele"="ref",
                                 "Tumor_Seq_Allele2"="alt")) %>%
  dplyr::group_by(Chromosome,Start_Position,Tumor_Seq_Allele2) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(genome.matched.unrestricted))
print(table(genome.matched.unrestricted$source,exclude = NULL))

### Unristricted match events - AA Change
aa.matched.unrestricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Hugo_Symbol"="gene",
                                 "AAChangeObserved"="AAChange")) %>%
  dplyr::group_by(Hugo_Symbol,AAChangeObserved) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(aa.matched.unrestricted))
print(table(aa.matched.unrestricted$source,exclude = NULL))

### impose cancer type matches (genome) - oncotree
genome.matched.ot.restricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome_annotation",
                                 "Start_Position"="pos",
                                 "Reference_Allele"="ref",
                                 "Tumor_Seq_Allele2"="alt",
                                 "ONCOTREE_CODE"="oncotree_code_annotation")) %>%
  dplyr::group_by(Chromosome,Start_Position,Tumor_Seq_Allele2,ONCOTREE_CODE) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   disease=paste0(unique(Disease),collapse=";"),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(genome.matched.unrestricted))
print(table(genome.matched.unrestricted$source,exclude = NULL))

### impose cancer type matches (protein) - oncotree
aa.matched.ot.restricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Hugo_Symbol"="gene",
                                 "AAChangeObserved"="AAChange",
                                 "ONCOTREE_CODE"="oncotree_code_annotation")) %>%
  dplyr::group_by(Hugo_Symbol,AAChangeObserved,ONCOTREE_CODE) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   disease=paste0(unique(Disease),collapse=";"),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(aa.matched.unrestricted))
print(table(aa.matched.unrestricted$source,exclude = NULL))

### Build reporting table
n.patients <- length(unique(svCompiled$Tumor_Sample_Barcode))
n.variants <- length(unique(dbAlterations$AAChangeObserved))

c1 <- c("Restricted","Restricted","Unrestricted","Unrestricted")
c2 <- c("genome","AA change", "genome", "AA change")
c3 <- c(dim(genome.matched.ot.restricted)[[1]], 
        dim(aa.matched.ot.restricted)[[1]], 
        dim(genome.matched.unrestricted)[[1]], 
        dim(aa.matched.unrestricted)[[1]])
#c4 <- paste0(round(100*(c3/n.variants),1),"%")
c5 <- c(sum(genome.matched.ot.restricted$n.patients.mutated), 
        sum(aa.matched.ot.restricted$n.patients.mutated), 
        sum(genome.matched.unrestricted$n.patients.mutated), 
        sum(aa.matched.unrestricted$n.patients.mutated))
c6 <- paste0(round(100*(c5/n.patients),1),"%") 
# Create a data frame
df <- data.frame(c1=c1, c2=c2, c3=c3, c5=c5, c6=c6)
colnames(df) <- c("Cancer type match restriction", "Cancer type matching approach", "number of variants","number of patients","percent of patients")
# Print the data frame
#print(df)
knitr::kable(df, )
write.csv(df,"../../output/actionability_db_curration_20231220/actionability_restricted_unrestricted_counts.csv")
