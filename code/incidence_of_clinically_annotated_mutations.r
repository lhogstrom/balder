library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)
library(ggrepel)
library(ggplot2)
library(ggalluvial)

### connect to result DB and get variants
outDir <- "../../output/clinical_annotation_matching_20240412"
bDir <- "../../data/processed/balderResultsDb"
dbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v20240412.sqlite")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), dbName)

# load clinical annotations
dbRules <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM MoaCiVICRuleEntries') %>%
  dplyr::rename(oncotree_code_annotation=oncotree_code,
                oncotree_term_annotation=oncotree_term,
                chromosome_annotation=chromosome)
dbRules$balderRuleID <- seq(1,dim(dbRules)[[1]])
dbRules$matchFlag <- "Y"
dbAlterations <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM actionableSNVsByAAChange')
dbGenome <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM actionableSNVsByGenomicCoordinate')
ot_code_full <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM OncoTreeCodeHierarchy')

# load patient observed variants
poVariants<- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM patientObservedVariantTable') %>%
  dplyr::mutate(AAChangeObserved=stringr::str_replace(HGVSp_Short,"p.", ""),
                MAF=100*(t_alt_count / (t_ref_count+t_alt_count)))
poVariants$balderVariantID <- seq(1,dim(poVariants)[[1]])
# replace MAF with TVAF for hartwig variants
iHartwig <- poVariants$SourceStudy == "Hartwig-data"
poVariants[iHartwig,"MAF"] <- poVariants[iHartwig,"TVAF"] 

print(dim(poVariants))

# load sample info for observed variants
sampleInfoCompiled <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM harmonizedSampleInfo')
sampleInfoCompiled[is.na(sampleInfoCompiled$SAMPLE_TYPE),"sampleInfoCompiled"] <- "Unspecified"
#compile hartwig and genie seq or report to death
sampleInfoCompiled$INT_Seq_or_Biopsy_To_Death <- sampleInfoCompiled$INT_Biopsy_To_Death
iNonNullSeqDeath <- !is.na(sampleInfoCompiled[,"INT_SEQ_TO_DEATH"])
sampleInfoCompiled[iNonNullSeqDeath,"INT_Seq_or_Biopsy_To_Death"] <- sampleInfoCompiled[iNonNullSeqDeath,"INT_SEQ_TO_DEATH"]


otCodeCols <- c("ot_code", "ot_name", "Highest_Non_Null_Level", "oncotree_level",
                "level_1_disease", "ot_code_level_1",
                "level_2_disease", "ot_code_level_2",
                "level_3_disease", "ot_code_level_3",
                "level_4_disease", "ot_code_level_4",
                "level_5_disease", "ot_code_level_5",
                "level_6_disease", "ot_code_level_6")

### join sample info to variant data
svCompiled <- poVariants %>%
  dplyr::inner_join(sampleInfoCompiled[,!colnames(sampleInfoCompiled)=="SourceStudy"],
                   by=c("Tumor_Sample_Barcode"="SAMPLE_ID")) %>% # exclude samples that do not have any sample information
  dplyr::left_join(ot_code_full[,otCodeCols],by=c("ONCOTREE_CODE"="ot_code")) %>% # join ot_codes to annotation
  dplyr::filter(!SAMPLE_TYPE == "Cell line")
print(dim(svCompiled))

######################################
### Build primary reporting tables ###
######################################

aa.match.exhaustive <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Hugo_Symbol"="gene",
                                 "AAChangeObserved"="AAChange"),
                    relationship = "many-to-many") %>%
  dplyr::mutate(matchAndVarID=paste0(balderVariantID,"-",balderRuleID)) %>%
  dplyr::mutate(AAChange=AAChangeObserved,
                gene=Hugo_Symbol) # create dummy variable for combining
  
genome.match.exhaustive <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome_annotation",
                                 "Start_Position"="pos",
                                 "Reference_Allele"="ref",
                                 "Tumor_Seq_Allele2"="alt"),
                    relationship = "many-to-many") %>%
  dplyr::mutate(matchAndVarID=paste0(balderVariantID,"-",balderRuleID)) %>%
  dplyr::mutate(chromosome_annotation=Chromosome,
                pos=Start_Position,
                ref=Reference_Allele,
                alt=Tumor_Seq_Allele2) # create dummy variable for combining

aa.genome.exahustive <- rbind(aa.match.exhaustive,genome.match.exhaustive) %>% 
  dplyr::group_by(balderVariantID,balderRuleID) %>%
  dplyr::filter(dplyr::row_number()==1) %>%
  dplyr::mutate(cancerTypeMatchPrimary=oncotree_code_annotation==ONCOTREE_CODE & !is.na(ONCOTREE_CODE), # ot code from variants (CAPS) & ot code from annotation (not caps)
                cancerTypeMatchLevel1=oncotree_code_annotation==ot_code_level_1 & !is.na(ot_code_level_1),
                cancerTypeMatchLevel2=oncotree_code_annotation==ot_code_level_2 & !is.na(ot_code_level_2),
                cancerTypeMatchLevel3=oncotree_code_annotation==ot_code_level_3 & !is.na(ot_code_level_3),
                cancerTypeMatchLevel4=oncotree_code_annotation==ot_code_level_4 & !is.na(ot_code_level_4),
                cancerTypeMatchLevel5=oncotree_code_annotation==ot_code_level_5 & !is.na(ot_code_level_5),
                cancerTypeMatchLevel6=oncotree_code_annotation==ot_code_level_6 & !is.na(ot_code_level_6),
                cancerTypeMatch = cancerTypeMatchPrimary | cancerTypeMatchLevel1 | cancerTypeMatchLevel2 | cancerTypeMatchLevel3 | cancerTypeMatchLevel4 | cancerTypeMatchLevel5 | cancerTypeMatchLevel6)

#%>% # filter entries if there is a AA match and genomic match
  #dplyr::mutate(annotationMatchGenomicCoord=matchAndVarID %in% genome.match.exhaustive$matchAndVarID,
                #annotationMatchAAChange=matchAndVarID %in% aa.match.exhaustive$matchAndVarID)
aa.genome.exahustive$annotationMatchGenomicCoord <- aa.genome.exahustive$matchAndVarID %in% genome.match.exhaustive$matchAndVarID
aa.genome.exahustive$annotationMatchAAChange <- aa.genome.exahustive$matchAndVarID %in% aa.match.exhaustive$matchAndVarID

# SQL write
RSQLite::dbWriteTable(harmonizedDb, "exhaustiveClinicalAnnotatedPatientVariants", aa.genome.exahustive,overwrite=T)

# count variants matched on genomic coord and/ or AA change
table(aa.genome.exahustive$annotationMatchGenomicCoord,aa.genome.exahustive$annotationMatchAAChange,exclude=NULL)

### AA mismatch events between observed and annotation names? 
# how often is there a mismatch b/t AA change calls for genomic-based targets? 
table(genome.match.exhaustive$AAChangeObserved == genome.match.exhaustive$AAChange)

# pick a single representative variant per subject
aa.genome.representative <- aa.genome.exahustive %>%
  dplyr::arrange(desc(cancerTypeMatch),
                 clinical.evidence.summary) %>% # prioritize restricted cancer type matches that are high-confidence
  dplyr::group_by(Tumor_Sample_Barcode,balderVariantID) %>%
  dplyr::filter(dplyr::row_number()==1) %>%
  dplyr::mutate(clinical_annotation_status="clinically annotated")
write.csv(aa.genome.representative,
          "../../output/actionability_db_curration_20231220/actionability_representative_variant_table.csv",
          row.names = F)
# SQL write
RSQLite::dbWriteTable(harmonizedDb, "representativeClinicalAnnotatedPatientVariants", aa.genome.representative,overwrite=T)


# count before and after representative variant selection
print(dim(aa.genome.exahustive))
print(dim(aa.genome.representative))
table(aa.genome.representative$cancerTypeMatch)

length(unique(genome.match.exhaustive$balderVariantID))
length(unique(genome.match.exhaustive$Tumor_Sample_Barcode))


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
                      c("number of subjects with at least one match",length(unique(genome.match.exhaustive$Tumor_Sample_Barcode))),
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


### Disconnect from SQL db ###
RSQLite::dbDisconnect(harmonizedDb)

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



##############################################
### mutant allele frequency (MAF) analysis ###
##############################################


iMafOver100 <- annotated.and.non.annotated$MAF > 100
outF <-  paste0(outDir,"/MAF_distribution_clinically_actionable_vs_non_actionable.pdf")
ggplot(annotated.and.non.annotated[!iMafOver100,],aes(x=MAF,fill=clinical_annotation_status))+
  #geom_histogram(alpha=.4)+
  geom_density(alpha=.4)+
  theme_bw()+
  facet_grid(SourceStudy~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of unrestricted annotated and non-annotated variants \n in source studies"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)

### Examine MAF distributions ###

### select cancer type unrestricted data for evaluation
#data_for_evaluation <- aa.genome.representative
#MafOutDir <- paste0(outDir,"/MAF_analysis_unrestricted_cancer_types")

### select cancer type restricted data for evaluation
data_for_evaluation <- aa.genome.representative %>%
 dplyr::filter(cancerTypeMatch==T)
MafOutDir <- paste0(outDir,"/MAF_analysis_restricted_cancer_types")


outF <-  paste0(MafOutDir,"/MAF_distribution_clinically_actionable.pdf")
ggplot(data_for_evaluation,aes(x=MAF,fill=actionability.summary))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  facet_grid(SourceStudy~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n in source studies"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 10,width = 10)


outF <-  paste0(MafOutDir,"/MAF_distribution_clinically_actionable_type.pdf")
ggplot(data_for_evaluation,aes(x=MAF,fill=actionability.summary))+
  geom_density(alpha=.4)+
  theme_bw()+
  facet_grid(SourceStudy~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n by type"))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 10,width = 10)


isMetPrim <- data_for_evaluation$SAMPLE_TYPE %in% c("Metastasis","Primary")
isMC3 <- data_for_evaluation$SourceStudy == "TCGA-MC3-data"
met.or.primary <- data_for_evaluation[isMetPrim & !isMC3,]

### density plots of actionability type,met/primary, and data source
outF <-  paste0(MafOutDir,"/MAF_distribution_clinically_actionable_type_primary_met.pdf")
ggplot(met.or.primary,aes(x=MAF,fill=actionability.summary))+
  geom_density(alpha=.4)+
  theme_bw()+
  #facet_grid(SourceStudy~SAMPLE_TYPE,scale="free_y",)+
  facet_grid(SAMPLE_TYPE~SourceStudy,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n by type"))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 10,width = 10)

outF <-  paste0(MafOutDir,"/MAF_distribution_met_vs_primary.pdf")
ggplot(met.or.primary,aes(x=MAF,fill=SAMPLE_TYPE))+
  geom_density(alpha=.4)+
  theme_bw()+
  #facet_grid(SourceStudy~SAMPLE_TYPE,scale="free_y",)+
  facet_grid(.~SourceStudy,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n by tissue type"))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 6,width = 10)

### create comparisons of distribution types
testResDf <- data.frame()

srcStudies <- unique(data_for_evaluation$SourceStudy)
#sampleTypes <- unique(data_for_evaluation$SAMPLE_TYPE)
sampleTypes <- c("Metastasis","Primary")
actionSummTypes <- unique(data_for_evaluation$actionability.summary)
for (testSource in srcStudies) {
  print(testSource)
  ixSrc <- data_for_evaluation$SourceStudy == testSource
  for (xMetPrim in sampleTypes) {
    ixMetPrim <- data_for_evaluation$SAMPLE_TYPE == xMetPrim
    for (i in 1:(length(actionSummTypes) - 1)) {
      for (j in (i + 1):length(actionSummTypes)) { 
        
        xAS <- actionSummTypes[i]
        yAS <- actionSummTypes[j]
        
        # conditions of test
        ixActionSumm <- data_for_evaluation$actionability.summary == xAS
        iyActionSumm <- data_for_evaluation$actionability.summary == yAS
        
        xVarSubset <- data_for_evaluation[ixSrc & ixMetPrim & ixActionSumm,]
        yVarSubset <- data_for_evaluation[ixSrc & ixMetPrim & iyActionSumm,]
        print(dim(xVarSubset))
        print(dim(yVarSubset))
        
        xSize <- dim(xVarSubset)[[1]]
        ySize <- dim(yVarSubset)[[1]]
        
        # skip statsitical test if there are fewer than N exemplars in a given group
        if ((xSize < 100) | (ySize < 100)) {
          next
        }
        
        # create vectors of MAF values
        xMAFVecTmp <- xVarSubset$MAF
        yMAFVecTmp <- yVarSubset$MAF
 
        ixNa <- is.na(xMAFVecTmp)
        iyNa <- is.na(yMAFVecTmp)
        ixInf <- is.infinite(xMAFVecTmp)
        iyInf <- is.infinite(yMAFVecTmp)
        
        # remove nulls and Inf from vectors
        xMAFVec <- xMAFVecTmp[!ixNa & !ixInf]
        yMAFVec <- yMAFVecTmp[!iyNa & !iyInf]
        
        # Organize data for density plot
        mafDf <- data.frame(value=c(xMAFVec,yMAFVec),
                            condition=xAS)
        iC2 <- seq((length(xMAFVec)+1),dim(mafDf)[[1]])
        mafDf[iC2,"condition"] <- yAS
        #mafDf$condition <- factor(mafDf$condition)
        
        ecdf1_data <- data.frame(x = unique(xMAFVec), y = ecdf(xMAFVec)(unique(xMAFVec)))
        ecdf2_data <- data.frame(x = unique(yMAFVecTmp), y = ecdf(yMAFVecTmp)(unique(yMAFVecTmp)))
        
        # Label the data
        ecdf1_data$Condition <- xAS
        ecdf2_data$Condition <- yAS
        
        # Combine the data frames
        ecdf_data <- rbind(ecdf1_data, ecdf2_data)
        
        # Calculate the KS statistic and the point of maximum difference
        ks_test <- ks.test(xMAFVec, yMAFVec)
        ks_stat <- ks_test$statistic
        xx <- seq(0,100,.1)
        diffs <- abs(ecdf(xMAFVec)(xx) - ecdf(yMAFVec)(xx))
        ks_value <- xx[which.max(diffs)]
        
        # Plot the empirical CDFs using ggplot2
        outF <-  paste0(MafOutDir,"/MAF_distributions_actionable_",xAS,"_",yAS,"_",xMetPrim,"_",testSource,"_ecdf.pdf")
        ggplot(ecdf_data, aes(x = x, y = y, color = Condition)) + 
          geom_line() +
          geom_vline(xintercept = ks_value, linetype = "dotted", color = "black", size = 1) +
          annotate("text", x = ks_value, y = 0.5, label = paste("KS Statistic = ", round(ks_stat, digits=4)), hjust = 1.2) +
          #scale_color_manual(values = c(xAS = 'blue', yAS = 'red')) +
          labs(title = "Empirical CDFs with KS Statistic", x = "Value", y = "ECDF") +
          theme_minimal()
        ggsave(outF,height = 6,width = 6)
        
        ## density plot
        outF <-  paste0(MafOutDir,"/MAF_distributions_actionable_",xAS,"_",yAS,"_",xMetPrim,"_",testSource,".pdf")
        ggplot(mafDf,aes(x=value,fill=condition))+
          geom_density(alpha=.4)+
          theme_bw()+
          geom_vline(xintercept = ks_value, linetype = "dotted", color = "black", size = 1)+ 
          xlab("MAF (%)")+
          ggtitle(paste0("MAF values in ",xMetPrim," tissue \n",testSource))+
          theme(plot.title = element_text(hjust = 0.5))
        ggsave(outF,height = 6,width = 6)
        
        
        # perform wilcox test
        wrs.twosided <- wilcox.test(xMAFVec,
                                    yMAFVec,
                                    alternative = "two.sided",
                                    exact=FALSE,
                                    conf.int=TRUE,
                                    paired=FALSE)
        
        ## create outputs
        c1_temp_summary <- summary(xMAFVec)
        c2_temp_summary <- summary(yMAFVec)
        
        # Create a row to append to the dataframe
        new_row <- data.frame(
          # condition info
          source_study = testSource,
          SampleType = xMetPrim,
          c1_actionability_summary = xAS,
          c2_actionability_summary = yAS,
          c1_variant_count = xSize,
          c2_variant_count = ySize,
          # summary stats condition 1
          c1_min = c1_temp_summary[[1]],
          c1_first_quartile = c1_temp_summary[[2]],
          c1_median = c1_temp_summary[[3]],
          c1_mean = c1_temp_summary[[4]],
          c1_third_quartile = c1_temp_summary[[5]],
          c1_max = c1_temp_summary[[6]],
          # summary stats condition 2
          c2_min = c2_temp_summary[[1]],
          c2_first_quartile = c2_temp_summary[[2]],
          c2_median = c2_temp_summary[[3]],
          c2_mean = c2_temp_summary[[4]],
          c2_third_quartile = c2_temp_summary[[5]],
          c2_max = c2_temp_summary[[6]],
          # wilcox results
          wilcox_p_value = wrs.twosided$p.value,
          wilcox_statistic = wrs.twosided$statistic[[1]],
          wilcox_estimate = wrs.twosided$estimate[[1]]
        )
        
        # Append the new row to the summary dataframe
        testResDf <- rbind(testResDf, new_row)
        
      }
    }
  }
}

# Apply FDR correction
testResDf$p_adjusted_wilcox <- p.adjust(testResDf$wilcox_p_value, method = "fdr")

testResDf <- testResDf %>%
  dplyr::arrange(wilcox_p_value) %>%
  dplyr::mutate(significant = p_adjusted_wilcox < 0.05,
                c1_c2_median_MAF_difference = c1_median - c2_median)
outF <- paste0(MafOutDir,"/actionability_type_MAF_wilcox_test.csv")
write.csv(testResDf,outF)

### plot directionality of relationships and significance
outF <-  paste0(MafOutDir,"/MAF_distribution_wilcox_test_results.pdf")
ggplot(testResDf,aes(x=source_study,y=wilcox_estimate,shape=significant,color=SampleType))+
  geom_point(alpha=.8)+
  theme_bw()+
  facet_grid(c1_actionability_summary~c2_actionability_summary,scale="free_y",)+
  ylim(-8,8)+
  geom_hline(yintercept = 0,linetype="dotted")+
  ggtitle(paste0("Wilcox tests restults of MAF distributions \n by actionability type"))+
  #scale_color_brewer(palette="Set1",drop=FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)

### plot wilcox estimate and median MAF difference 
outF <-  paste0(MafOutDir,"/MAF_wilcox_estimate_and_median_diff.pdf")
ggplot(testResDf,aes(x=wilcox_estimate,y=c1_c2_median_MAF_difference,shape=significant,color=SampleType))+
  geom_point(alpha=.8)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 5,width = 5)

### top variant summary
atGroupVariantSummary <- data_for_evaluation %>%
  dplyr::group_by(actionability.summary,Hugo_Symbol,HGVSp_Short) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::arrange(desc(n))

atGroupVariantSummaryTop <- atGroupVariantSummary %>%
  dplyr::group_by(actionability.summary) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::filter(dplyr::row_number() < 16) %>%
  dplyr::mutate(variantName=paste0(Hugo_Symbol," ",HGVSp_Short))

outF <-  paste0(MafOutDir,"/top_variants_by_actionability_type_unristricted.pdf")
ggplot(atGroupVariantSummaryTop,aes(x=actionability.summary,y=n))+
  geom_point(alpha=.8)+
  ggtitle(paste0("Top 15 variants by actionability type \n unrestricted cancer type"))+
  #geom_text(aes(label = variantName), nudge_x = 0.0, nudge_y = 0.2)+
  geom_text_repel(aes(label = variantName))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 9,width = 9)

outF <-  paste0(MafOutDir,"/top_variants_by_actionability_type_unristricted_log.pdf")
ggplot(atGroupVariantSummaryTop[atGroupVariantSummaryTop$n>1,],aes(x=actionability.summary,y=n))+
  geom_point(alpha=.8)+
  #geom_jitter(alpha=.8,width = .1)+
  ggtitle(paste0("Top 15 variants by actionability type"))+
  #geom_text(aes(label = variantName), nudge_x = 0.0, nudge_y = 0.2)+
  geom_text_repel(aes(label = variantName))+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 9,width = 12)

table(data_for_evaluation$actionability.summary,data_for_evaluation$cancerTypeMatch)

#####################################
### sequencing or report to death ###
#####################################

# compare GENIE to Hartwig 
outF <-  paste0(outDir,"/hartwig_vs_aacr_time_to_death.pdf")
ggplot(sampleInfoCompiled,aes(x=INT_Seq_or_Biopsy_To_Death,fill=SourceStudy))+
  geom_density(alpha=.4)+
  theme_bw()+
  xlim(-500,2000)+
  ggtitle(paste0("Infered days between biopsy or sequencing \nreport and death"))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 6,width = 6)

outF <-  paste0(outDir,"/hartwig_vs_aacr_time_to_death_sample_type.pdf")
iMetPrim = (sampleInfoCompiled$SAMPLE_TYPE=="Primary" | sampleInfoCompiled$SAMPLE_TYPE=="Metastasis") & !is.na(sampleInfoCompiled$SAMPLE_TYPE)
ggplot(sampleInfoCompiled[iMetPrim,],aes(x=INT_Seq_or_Biopsy_To_Death,fill=SourceStudy))+
  geom_density(alpha=.4)+
  facet_grid(SAMPLE_TYPE~.)+
  theme_bw()+
  xlim(-500,2000)+
  ggtitle(paste0("Infered days between biopsy or sequencing \nreport and death"))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 6,width = 6)

####################################
### aggregate prognostic markers ###
####################################

# which samples have a prognostic marker: high conf, low, conf
prognosticRepCnts <- aa.genome.representative %>% 
  dplyr::filter(cancerTypeMatch==T) %>%
  dplyr::group_by(Tumor_Sample_Barcode) %>%
  dplyr::summarise(nVariants=n(),
                   nPrognosticVars=sum(actionability.summary == "Prognostic"),
                   nMarkerBetterOutcome=sum(clinical_significance == "Better Outcome"),
                   nMarkerPoorOutcome=sum(clinical_significance == "Poor Outcome"),
                   nPrognosticHighConfVars=sum((actionability.summary == "Prognostic") & (clinical.evidence.summary=="High Confidence")),
                   nPrognosticLowConfVarssum=sum((actionability.summary == "Prognostic") & (clinical.evidence.summary=="Lower Confidence")))


# some samples have both poor outcome and better outcome prognostic markers
table(prognosticRepCnts$nMarkerBetterOutcome,prognosticRepCnts$nMarkerPoorOutcome)

# make prognostic marker assignments by sample
iProgMarkers <- prognosticRepCnts$nPrognosticVars > 0
sampleInfoCompiled$hasPrognosticMarker <- sampleInfoCompiled$SAMPLE_ID %in% prognosticRepCnts[iProgMarkers,"Tumor_Sample_Barcode"][[1]]

sampleInfoCompiled$prognosticMarkerType <- "Prognostic Marker Unspecified"
iBetterProgMarkers <- prognosticRepCnts$nMarkerBetterOutcome > 0
iBetter <- sampleInfoCompiled$SAMPLE_ID %in% prognosticRepCnts[iBetterProgMarkers,"Tumor_Sample_Barcode"][[1]]
sampleInfoCompiled[iBetter,"prognosticMarkerType"] <- "Better Outcome Prognostic Marker"
iPoorProgMarkers <- prognosticRepCnts$nMarkerPoorOutcome > 0
iPoor <- sampleInfoCompiled$SAMPLE_ID %in% prognosticRepCnts[iPoorProgMarkers,"Tumor_Sample_Barcode"][[1]]
sampleInfoCompiled[iPoor,"prognosticMarkerType"] <- "Poor Outcome Prognostic Marker"
sampleInfoCompiled[!sampleInfoCompiled$hasPrognosticMarker,"prognosticMarkerType"] <- "No Prognostic Marker"
  
outF <-  paste0(outDir,"/hartwig_vs_aacr_time_to_death_sample_type.pdf")
iMetPrim = (sampleInfoCompiled$SAMPLE_TYPE=="Primary" | sampleInfoCompiled$SAMPLE_TYPE=="Metastasis") & !is.na(sampleInfoCompiled$SAMPLE_TYPE) & !is.na(sampleInfoCompiled$INT_Seq_or_Biopsy_To_Death)
ggplot(sampleInfoCompiled[iMetPrim,],aes(x=INT_Seq_or_Biopsy_To_Death,fill=prognosticMarkerType))+
  geom_density(alpha=.4)+
  facet_grid(SourceStudy~.)+
  theme_bw()+
  xlim(-500,2000)+
  ggtitle(paste0("Infered days between biopsy or sequencing \nreport and death"))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 6,width = 8)

####################################
### specific prognostic markers ###
####################################

prognosticTypeMutCnts <- aa.genome.representative %>%
  dplyr::filter(cancerTypeMatch==T,
                actionability.summary == "Prognostic") %>%
  dplyr::group_by(ONCOTREE_CODE,CANCER_TYPE,gene,AAChangeObserved,clinical_significance) %>%
  dplyr::summarise(n=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n))
outF <- paste0(ProgOutDir,"/prognostic_variant_frequency_rank.csv")
write.csv(prognosticTypeMutCnts,outF)

prognosticTypeCnts <- aa.genome.representative %>%
  dplyr::filter(cancerTypeMatch==T,
                actionability.summary == "Prognostic") %>%
  dplyr::group_by(CANCER_TYPE,gene,clinical_significance) %>%
  dplyr::summarise(n=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(n))

ProgOutDir <- "../../output/clinical_annotation_matching_20240412/prognostic_analysis_restricted_cancer_types"
### create comparisons of distribution types
testResProgDf <- data.frame()

#otCodes <- unique(prognosticTypeCnts$ONCOTREE_CODE)
ctTypes <- unique(prognosticTypeCnts$CANCER_TYPE)
progGenes <- unique(prognosticTypeCnts$gene)
progSig <- unique(prognosticTypeCnts$clinical_significance)

for (iX in rownames(prognosticTypeCnts)) {
  print(iX)
  ctType <- prognosticTypeCnts[iX,"CANCER_TYPE"][[1]]
  progGene <- prognosticTypeCnts[iX,"gene"][[1]]
  progSigDir <- prognosticTypeCnts[iX,"clinical_significance"][[1]]
  print(paste(ctType,"-",progGene,"-",progSigDir))
  
  # select variants with the prognostic markers
  xVarSubset <- aa.genome.representative %>%
    dplyr::filter(cancerTypeMatch==T,
                  actionability.summary == "Prognostic",
                  CANCER_TYPE == ctType,
                  gene == progGene,
                  clinical_significance == progSigDir)
  # all other variants for the cancer type
  yVarSubset <- aa.genome.representative %>%
    dplyr::filter(CANCER_TYPE == ctType,
                  !Tumor_Sample_Barcode %in% xVarSubset$Tumor_Sample_Barcode)
  
  # select the samples with the variants (and the metric)
  xSubSetSamples <- xVarSubset %>%
    dplyr::group_by(Tumor_Sample_Barcode) %>%
    dplyr::filter(dplyr::row_number()==1)
  
  ySubSetSamples <- yVarSubset %>%
    dplyr::group_by(Tumor_Sample_Barcode) %>%
    dplyr::filter(dplyr::row_number()==1)
  
  xSize <- dim(xSubSetSamples)[[1]]
  ySize <- dim(ySubSetSamples)[[1]]
  
  # skip statsitical test if there are fewer than N exemplars in a given group
  if ((xSize < 30) | (ySize < 30)) {
    next
  }
  
  # create vectors of MAF values
  xSeqToDeathTmp <- xSubSetSamples$INT_Seq_or_Biopsy_To_Death
  ySeqToDeathTmp <- ySubSetSamples$INT_Seq_or_Biopsy_To_Death
  
  ixNa <- is.na(xSeqToDeathTmp)
  iyNa <- is.na(ySeqToDeathTmp)
  ixInf <- is.infinite(xSeqToDeathTmp)
  iyInf <- is.infinite(ySeqToDeathTmp)
  
  # remove nulls and Inf from vectors
  xSeqToDeathVec <- xSeqToDeathTmp[!ixNa & !ixInf]
  ySeqToDeathVec <- ySeqToDeathTmp[!iyNa & !iyInf]
  
  # Organize data for density plot
  mafDf <- data.frame(value=c(xSeqToDeathVec,ySeqToDeathVec),
                      condition="sample has prognostic marker")
  iC2 <- seq((length(xSeqToDeathVec)+1),dim(mafDf)[[1]])
  mafDf[iC2,"condition"] <- "no prognostic marker in sample"
  
  # Calculate the KS statistic and the point of maximum difference
  ks_test <- ks.test(xSeqToDeathVec, ySeqToDeathVec)
  ks_stat <- ks_test$statistic
  xx <- seq(0,100,.1)
  diffs <- abs(ecdf(xSeqToDeathVec)(xx) - ecdf(ySeqToDeathVec)(xx))
  ks_value <- xx[which.max(diffs)]
  
  outF <-  paste0(ProgOutDir,"/days_from_seq_to_death_",ctType,"_",progGene,"_",progSigDir,".pdf")
  ggplot(mafDf,aes(x=value,fill=condition))+
    geom_density(alpha=.4)+
    theme_bw()+
    geom_vline(xintercept = ks_value, linetype = "dotted", color = "black", size = 1)+
    xlab("Days between sequencing and death")+
    ggtitle(paste0("Death outcome for ",ctType,"\n with prognostic ",progGene," mutation \n expected direction: ",progSigDir))+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(outF,height = 6,width = 6)
  
  # perform wilcox test
  wrs.twosided <- wilcox.test(xSeqToDeathVec,
                              ySeqToDeathVec,
                              alternative = "two.sided",
                              exact=FALSE,
                              conf.int=TRUE,
                              paired=FALSE)
  
  ## create outputs
  c1_temp_summary <- summary(xSeqToDeathVec)
  c2_temp_summary <- summary(ySeqToDeathVec)
  
  # Create a row to append to the dataframe
  new_row <- data.frame(
    # condition info
    cancer_type = ctType,
    prognostic_marker_gene=progGene,
    expected_significance=progSigDir,
    # summary stats condition 1
    c1_min = c1_temp_summary[[1]],
    c1_first_quartile = c1_temp_summary[[2]],
    c1_median = c1_temp_summary[[3]],
    c1_mean = c1_temp_summary[[4]],
    c1_third_quartile = c1_temp_summary[[5]],
    c1_max = c1_temp_summary[[6]],
    # summary stats condition 2
    c2_min = c2_temp_summary[[1]],
    c2_first_quartile = c2_temp_summary[[2]],
    c2_median = c2_temp_summary[[3]],
    c2_mean = c2_temp_summary[[4]],
    c2_third_quartile = c2_temp_summary[[5]],
    c2_max = c2_temp_summary[[6]],
    # wilcox results
    wilcox_p_value = wrs.twosided$p.value,
    wilcox_statistic = wrs.twosided$statistic[[1]],
    wilcox_estimate = wrs.twosided$estimate[[1]]
  )
  
  # Append the new row to the summary dataframe
  testResProgDf <- rbind(testResProgDf, new_row)
}

# Apply FDR correction
testResProgDf$p_adjusted_wilcox <- p.adjust(testResProgDf$wilcox_p_value, method = "fdr")

testResProgDf <- testResProgDf %>%
  dplyr::arrange(wilcox_p_value) %>%
  dplyr::mutate(significant = p_adjusted_wilcox < 0.05,
                c1_c2_median_MAF_difference = c1_median - c2_median)
outF <- paste0(ProgOutDir,"/prognostic_markers_days_from_seq_to_death_wilcox_test.csv")
write.csv(testResProgDf,outF)




#######################################################
### Sankey diagram of biomarker to drug assignments ###
#######################################################

# data <- data.frame(
#   id = 1:4,
#   stage1 = c("A", "A", "B", "B"),
#   stage2 = c("C", "D", "D", "C"),
#   stage3 = c("E", "E", "F", "F"),
#   weight = c(1, 2, 1, 1)
# )
# 
# ggplot(data = data,
#        aes(axis1 = stage1, axis2 = stage2, axis3 = stage3, y = weight)) +
#   geom_alluvium(aes(fill = stage1)) +  # Fill color based on the first stage
#   geom_stratum() +  # Add stratum to show stages
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +  # Add text labels
#   theme_minimal() +  # Use a minimal theme
#   ggtitle("Sankey Diagram with Three Stages")  # Add a title

### Sankey diagram of cancer type --> biomarker --> prognostic assignment
