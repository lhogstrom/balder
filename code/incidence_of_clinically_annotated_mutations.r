library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)
library(ggrepel)
library(ggplot2)

### connect to result DB and get variants
bDir <- "../../data/processed/balderResultsDb"
dbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v20240311.sqlite")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), dbName)

# load clinical annotations
dbRules <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM MoaCiVICRuleEntries')
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
print(dim(poVariants))

# load sample info for observed variants
sampleInfoCompiled <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM harmonizedSampleInfo')
sampleInfoCompiled[is.na(sampleInfoCompiled$SAMPLE_TYPE),"sampleInfoCompiled"] <- "Unspecified"
print(dim(sampleInfoCompiled))
RSQLite::dbDisconnect(harmonizedDb)


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
                   by=c("Tumor_Sample_Barcode"="SAMPLE_ID")) #%>%
  #dplyr::left_join(ot_code_full[,otCodeCols],by=c("oncotree_code"="ot_code")) # join ot_codes to annotation entry
  #dplyr::left_join(sampleInfoCompiled[,!colnames(sampleInfoCompiled)=="SourceStudy"],
  #                 by=c("Tumor_Sample_Barcode"="SAMPLE_ID"))
print(dim(svCompiled))
table(svCompiled$SAMPLE_TYPE,exclude=NULL)


# Computing actionable variant counts for different match approaches
# 
# Different actionablity matching approaches were computed:
#   
# Unrestricted where no cancer type matching was imposed
# 
# Restricted with raw cancer type string matching

dbRules <- dbRules %>%
  dplyr::left_join(ot_code_full[,otCodeCols],by=c("oncotree_code"="ot_code")) # join ot_codes to annotation entry

###
### join with observed variant table with clinically annotated - by genome position
genome.matched.unrestricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome",
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
  dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome",
                                 "Start_Position"="pos",
                                 "Reference_Allele"="ref",
                                 "Tumor_Seq_Allele2"="alt",
                                 "ONCOTREE_CODE"="oncotree_code")) %>%
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
                                 "ONCOTREE_CODE"="oncotree_code")) %>%
  dplyr::group_by(Hugo_Symbol,AAChangeObserved,ONCOTREE_CODE) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   disease=paste0(unique(Disease),collapse=";"),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(aa.matched.unrestricted))
print(table(aa.matched.unrestricted$source,exclude = NULL))

### Reporting table

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
  dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome",
                                 "Start_Position"="pos",
                                 "Reference_Allele"="ref",
                                 "Tumor_Seq_Allele2"="alt"),
                    relationship = "many-to-many") %>%
  dplyr::mutate(matchAndVarID=paste0(balderVariantID,"-",balderRuleID)) %>%
  dplyr::mutate(chromosome=Chromosome,
                pos=Start_Position,
                ref=Reference_Allele,
                alt=Tumor_Seq_Allele2) # create dummy variable for combining

aa.genome.exahustive <- rbind(aa.match.exhaustive,genome.match.exhaustive) %>% 
  dplyr::group_by(balderVariantID,balderRuleID) %>%
  dplyr::filter(dplyr::row_number()==1) %>%
  dplyr::mutate(cancerTypeMatchPrimary=oncotree_code==ONCOTREE_CODE & !is.na(ONCOTREE_CODE), # ot code from variants (CAPS) & ot code from annotation (not caps)
                cancerTypeMatchLevel1=ONCOTREE_CODE==ot_code_level_1 & !is.na(ot_code_level_1),
                cancerTypeMatchLevel2=ONCOTREE_CODE==ot_code_level_2 & !is.na(ot_code_level_2),
                cancerTypeMatchLevel3=ONCOTREE_CODE==ot_code_level_3 & !is.na(ot_code_level_3),
                cancerTypeMatchLevel4=ONCOTREE_CODE==ot_code_level_4 & !is.na(ot_code_level_4),
                cancerTypeMatchLevel5=ONCOTREE_CODE==ot_code_level_5 & !is.na(ot_code_level_5),
                cancerTypeMatchLevel6=ONCOTREE_CODE==ot_code_level_6 & !is.na(ot_code_level_6),
                cancerTypeMatch = cancerTypeMatchPrimary | cancerTypeMatchLevel1 | cancerTypeMatchLevel2 | cancerTypeMatchLevel3 | cancerTypeMatchLevel4 | cancerTypeMatchLevel5 | cancerTypeMatchLevel6)

#%>% # filter entries if there is a AA match and genomic match
  #dplyr::mutate(annotationMatchGenomicCoord=matchAndVarID %in% genome.match.exhaustive$matchAndVarID,
                #annotationMatchAAChange=matchAndVarID %in% aa.match.exhaustive$matchAndVarID)
aa.genome.exahustive$annotationMatchGenomicCoord <- aa.genome.exahustive$matchAndVarID %in% genome.match.exhaustive$matchAndVarID
aa.genome.exahustive$annotationMatchAAChange <- aa.genome.exahustive$matchAndVarID %in% aa.match.exhaustive$matchAndVarID

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
  dplyr::arrange(desc(n.patients))
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

outDir <- "../../output/actionability_db_curration_20231220"


iMafOver100 <- annotated.and.non.annotated$MAF > 100
outF <-  paste0(outDir,"/MAF_distribution_clinically_actionable_vs_non_actionable.pdf")
ggplot(annotated.and.non.annotated[!iMafOver100,],aes(x=MAF,fill=clinical_annotation_status))+
  #geom_histogram(alpha=.4)+
  geom_density(alpha=.4)+
  theme_bw()+
  facet_grid(SourceStudy~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of annotated and non-annotated variants \n in source studies"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 8,width = 7)


### Examine MAF distributions

outF <-  paste0(outDir,"/MAF_distribution_clinically_actionable.pdf")
ggplot(aa.genome.representative,aes(x=MAF,fill=actionability.summary))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  facet_grid(SourceStudy~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n in source studies"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 10,width = 10)


outF <-  paste0(outDir,"/MAF_distribution_clinically_actionable_type.pdf")
ggplot(aa.genome.representative,aes(x=MAF,fill=actionability.summary))+
  geom_density(alpha=.4)+
  theme_bw()+
  facet_grid(SourceStudy~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n by type"))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 10,width = 10)


isMetPrim <- aa.genome.representative$SAMPLE_TYPE %in% c("Metastasis","Primary")
isMC3 <- aa.genome.representative$SourceStudy == "TCGA-MC3-data"
met.or.primary <- aa.genome.representative[isMetPrim & !isMC3,]

### density plots of actionability type,met/primary, and data source
outF <-  paste0(outDir,"/MAF_distribution_clinically_actionable_type_primary_met.pdf")
ggplot(met.or.primary,aes(x=MAF,fill=actionability.summary))+
  geom_density(alpha=.4)+
  theme_bw()+
  #facet_grid(SourceStudy~SAMPLE_TYPE,scale="free_y",)+
  facet_grid(SAMPLE_TYPE~SourceStudy,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n by type"))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 10,width = 10)

outF <-  paste0(outDir,"/MAF_distribution_met_vs_primary.pdf")
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

srcStudies <- unique(aa.genome.representative$SourceStudy)
#sampleTypes <- unique(aa.genome.representative$SAMPLE_TYPE)
sampleTypes <- c("Metastasis","Primary")
actionSummTypes <- unique(aa.genome.representative$actionability.summary)
for (testSource in srcStudies) {
  print(testSource)
  ixSrc <- aa.genome.representative$SourceStudy == testSource
  for (xMetPrim in sampleTypes) {
    ixMetPrim <- aa.genome.representative$SAMPLE_TYPE == xMetPrim
    for (i in 1:(length(actionSummTypes) - 1)) {
      for (j in (i + 1):length(actionSummTypes)) { 
        
        xAS <- actionSummTypes[i]
        yAS <- actionSummTypes[j]
        
        # conditions of test
        ixActionSumm <- aa.genome.representative$actionability.summary == xAS
        iyActionSumm <- aa.genome.representative$actionability.summary == yAS
        
        xVarSubset <- aa.genome.representative[ixSrc & ixMetPrim & ixActionSumm,]
        yVarSubset <- aa.genome.representative[ixSrc & ixMetPrim & iyActionSumm,]
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
        outF <-  paste0(outDir,"/MAF_distributions_actionable_",xAS,"_",yAS,"_",xMetPrim,"_",testSource,"_ecdf.pdf")
        ggplot(ecdf_data, aes(x = x, y = y, color = Condition)) + 
          geom_line() +
          geom_vline(xintercept = ks_value, linetype = "dotted", color = "black", size = 1) +
          annotate("text", x = ks_value, y = 0.5, label = paste("KS Statistic = ", round(ks_stat, digits=4)), hjust = 1.2) +
          #scale_color_manual(values = c(xAS = 'blue', yAS = 'red')) +
          labs(title = "Empirical CDFs with KS Statistic", x = "Value", y = "ECDF") +
          theme_minimal()
        ggsave(outF,height = 6,width = 6)
        
        ## density plot
        outF <-  paste0(outDir,"/MAF_distributions_actionable_",xAS,"_",yAS,"_",xMetPrim,"_",testSource,".pdf")
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

write.csv(testResDf,"../../output/actionability_db_curration_20231220/actionability_type_MAF_wilcox_test.csv")

### plot directionality of relationships and significance
outF <-  paste0(outDir,"/MAF_distribution_wilcox_test_results.pdf")
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
outF <-  paste0(outDir,"/MAF_wilcox_estimate_and_median_diff.pdf")
ggplot(testResDf,aes(x=wilcox_estimate,y=c1_c2_median_MAF_difference,shape=significant,color=SampleType))+
  geom_point(alpha=.8)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 5,width = 5)

### top variant summary
atGroupVariantSummary <- aa.genome.representative %>%
  dplyr::group_by(actionability.summary,Hugo_Symbol,HGVSp_Short) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::arrange(desc(n))

atGroupVariantSummaryTop <- atGroupVariantSummary %>%
  dplyr::group_by(actionability.summary) %>%
  dplyr::arrange(desc(n)) %>%
  dplyr::filter(dplyr::row_number() < 16) %>%
  dplyr::mutate(variantName=paste0(Hugo_Symbol," ",HGVSp_Short))

outF <-  paste0(outDir,"/top_variants_by_actionability_type_unristricted.pdf")
ggplot(atGroupVariantSummaryTop,aes(x=actionability.summary,y=n))+
  geom_point(alpha=.8)+
  ggtitle(paste0("Top 15 variants by actionability type \n unrestricted cancer type"))+
  #geom_text(aes(label = variantName), nudge_x = 0.0, nudge_y = 0.2)+
  geom_text_repel(aes(label = variantName))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 9,width = 9)

table(aa.genome.representative$actionability.summary,aa.genome.representative$cancerTypeMatch)



