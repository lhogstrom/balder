#library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

# load variant DB

figDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/mutation_incidence_20220913"

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-VariantSummaries.tsv"
# load and order by variant score
df <- read.csv(inFile,sep="\t") %>%
  dplyr::arrange(desc(civic_variant_evidence_score))

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-AssertionSummaries.tsv"
assertions <- read.csv(inFile,sep="\t")

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"
clinical <- read.csv(inFile,sep="\t")
clinical$chromosomeNum <- as.character(clinical$chromosome)
clinical$chr <- paste0("chr",clinical$chromosomeNum)

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-VariantGroupSummaries.tsv"
varGroups <- read.csv(inFile,sep="\t")

# COMMAND ----------
### load MOAlmanac which is similar to CIViC
inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/MOAlmanac_43018_2021_243_MOESM2_ESM.txt"
moa <- read.csv(inFile,sep="\t")

print(table(moa$gene))
print(table(moa$predictive_implication))
print(table(moa$source_type))
print(table(moa$feature_type))
print(table(moa$variant_annotation))
print(table(moa$cosmic_signature_number,exclude=NULL))

signatures.moa <- moa[moa$feature_type=="Mutational signature",]

### make plot of feature type 
cnts.moa.type <- moa %>%
  dplyr::group_by(feature_type) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::arrange(desc(n))
cnts.moa.type$feature_type <- factor(cnts.moa.type$feature_type,levels=cnts.moa.type$feature_type)

ggplot(cnts.moa.type,aes(x=feature_type,y=n))+
  #geom_bar(stat_count="identity")+
  geom_point(alpha=.8,color='darkblue')+
  theme_bw()+
  ylab("nuber of MOA entries")+
  ggtitle(paste0("MOAlmanac aberrations by type"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### load MSK-IMPACT
inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/impact_2017_annotated_per_variant.tsv"
msk <- read.csv(inFile,sep="\t")
#inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/msk_impact_data_clinical_sample2.txt"
#patient.msk <- read.csv(inFile,sep="\t")
inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/msk_impact_2017/data_clinical_sample.txt"
patient.msk <- read.csv(inFile,sep="\t",skip = 4)
msk$Chromosome.str <- paste0("chr",msk$Chromosome)
msk$MAF <- 100*(msk$t_alt_count / (msk$t_ref_count+msk$t_alt_count))

### Is there a primary and secondary alt listed? 
#table(is.na(msk$n_alt_count),exclude=NULL)
#table(msk$Reference_Allele == msk$Tumor_Seq_Allele1,exclude=NULL)

# join clinical info onto variant info
all(msk$Tumor_Sample_Barcode %in% patient.msk$SAMPLE_ID)

msk <- msk %>% 
  dplyr::left_join(patient.msk,by=c("Tumor_Sample_Barcode"="SAMPLE_ID"))

# COMMAND ----------
### create actionability type assignemnts 
## MOA
moa$actionability.summary <- "Other"
isPrognostic <- moa$therapy_type=="" & !is.na(moa$favorable_prognosis)
moa[isPrognostic,"actionability.summary"] <- "Prognostic"
#
isPredTherapyAct <- is.na(moa$therapy_resistance) & !moa$therapy_type==""
moa[isPredTherapyAct,"actionability.summary"] <- "Predictive therapy actionability"
#
isTherapyResistance <- !is.na(moa$therapy_resistance) & !moa$therapy_type==""
moa[isTherapyResistance,"actionability.summary"] <- "Therapy resistance"
table(moa$actionability.summary,exclude=NULL)

## CIViC
clinical$actionability.summary <- "Other"
cSigPrognostic <- c("Better Outcome","N/A","Negative","Positive","Poor Outcome")
isPrognostic <- clinical$evidence_type=="Prognostic" & (clinical$clinical_significance %in% cSigPrognostic)
clinical[isPrognostic,"actionability.summary"] <- "Prognostic"
#
isPredTherapyAct <- clinical$evidence_type=="Predictive" | (clinical$clinical_significance=="Sensitivity/Response")
clinical[isPredTherapyAct,"actionability.summary"] <- "Predictive therapy actionability"
#
isTherapyResistance <- clinical$clinical_significance=="Resistance"
clinical[isTherapyResistance,"actionability.summary"] <- "Therapy resistance"
table(clinical$actionability.summary,exclude=NULL)

### make evidence assignments 
moa$clinical.evidence.summary <- "Lower Confidence"
moa[moa$source_type == "FDA" | moa$source_type == "Guideline","clinical.evidence.summary"] <- "High Confidence"
table(moa$clinical.evidence.summary)

clinical$clinical.evidence.summary <- "Lower Confidence"
clinical[clinical$evidence_level == "A" | clinical$evidence_level == "B","clinical.evidence.summary"] <- "High Confidence"
table(clinical$clinical.evidence.summary)




# COMMAND ----------
### add matching columns
moa$source <- "MOAlmanac"
moa$FDAApproved <- NA
moa$Indicated <- moa$predictive_implication
moa$EvidenceText <- moa$description
moa$ReferenceOrTrialID <- moa$citation
moa$Disease <- moa$disease
moa$TrialID <- NA
moa$Phase <- NA
moa$Drugs <- moa$therapy_name
moa$AberrationsLabelsMod <- NA
moa <- moa %>%
  tidyr::separate(protein_change,sep="\\.",into=c("nullP","AAChange"),remove=F)
moa$chr <- paste0("chr",moa$chromosome)
moa$pos <- moa$start_position
moa$ref <- moa$reference_allele
moa$alt <- moa$alternate_allele

clinical$source <- "CIViC-clinical"
clinical$FDAApproved <- NA
clinical$Indicated <- clinical$evidence_type
clinical$EvidenceText <- clinical$evidence_statement
clinical$ReferenceOrTrialID <- clinical$gene_civic_url
clinical$Disease <- clinical$disease
clinical$TrialID <- NA
clinical$Phase <- NA
clinical$Drugs <- clinical$drugs
clinical$AberrationsLabelsMod <- NA
clinical$AAChange <- clinical$variant
clinical$pos <- clinical$start
clinical$ref <- clinical$reference_bases
clinical$alt <- clinical$variant_bases


### combine MOA and CIViC data

### create cancer type maps CIVIC-->MSK and MOA-->MSK
civic.cType.cnt <- clinical %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(number.of.subjects=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.subjects))

moa.cType.cnt <- moa %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(number.of.subjects=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.subjects))

### which are the most common cancer types in msk, moa, and civic
msk.cType.cnt <- msk %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(number.of.subjects=dplyr::n_distinct(PATIENT_ID)) %>%
  dplyr::arrange(desc(number.of.subjects))
outF <-  paste0(figDir,"/msk_cancer_types.txt")
write.table(msk.cType.cnt,outF,row.names=F,quote=F,sep="\t")

## create draft cancer type maps
iCMmatch <- civic.cType.cnt$Disease %in% msk.cType.cnt$CANCER_TYPE
CivicMSKCancerTypeMap <- data.frame(CivicCancerType=civic.cType.cnt$Disease)
rownames(CivicMSKCancerTypeMap) <- CivicMSKCancerTypeMap$CivicCancerType
CivicMSKCancerTypeMap[iCMmatch,"MskCancerType"] <- as.character(civic.cType.cnt$Disease[iCMmatch])
mapF <-  paste0(figDir,"/civic_cancer_type_map.txt")
write.table(CivicMSKCancerTypeMap,mapF,row.names=F,quote=F,sep="\t")

### load refined map
#mapF <- paste0(figDir,"/civic_to_msk_cancer_type_map.tsv")
mapF <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Pancan_biomarker_incidence_v2/civic_cancer_type_map.txt"
CivicMSKCancerTypeMap <- read.csv(mapF,sep="\t")
clinical <- clinical %>%
  dplyr::left_join(CivicMSKCancerTypeMap,by=c("disease"="CivicCancerType"))
table(clinical$MskCancerType,exclude=NULL)
# what fraction of entries are still missing MSK cancer type assignment
print(sum(is.na(clinical$MskCancerType))/dim(clinical)[[1]])

### create draft map 
iCMmatch <- moa.cType.cnt$disease %in% msk.cType.cnt$CANCER_TYPE
MoaMSKCancerTypeMap <- data.frame(MoaCancerType=moa.cType.cnt$disease)
rownames(MoaMSKCancerTypeMap) <- MoaMSKCancerTypeMap$CivicCancerType
MoaMSKCancerTypeMap[iCMmatch,"MskCancerType"] <- as.character(moa.cType.cnt$disease[iCMmatch])
mapF <-  paste0(figDir,"/moa_cancer_type_map.txt")
write.table(MoaMSKCancerTypeMap,mapF,row.names=F,quote=F,sep="\t")

### load refined map
#mapF <-  paste0(figDir,"/moa_cancer_type_map.tsv")
mapF <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Pancan_biomarker_incidence_v2/moa_cancer_type_map.tsv"
MoaMSKCancerTypeMap <- read.csv(mapF,sep="\t")
MoaMSKCancerTypeMap$MoaCancerType <- as.character(MoaMSKCancerTypeMap$MoaCancerType)
moa <- moa %>%
  dplyr::left_join(MoaMSKCancerTypeMap,by=c("disease"="MoaCancerType"))
table(moa$MskCancerType,exclude=NULL)
# what fraction of entries are still missing MSK cancer type assignment
print(sum(is.na(moa$MskCancerType))/dim(moa)[[1]])

# now combine MOA and CIVIC
matchCols <- c("source","gene","AAChange","Drugs","FDAApproved","ReferenceOrTrialID","EvidenceText","Phase",
               "Indicated","Disease","chr","pos","ref","alt","MskCancerType","actionability.summary","clinical.evidence.summary")
dbRules <- rbind(moa[,matchCols],clinical[,matchCols]) %>%
  dplyr::mutate(EvidenceText=gsub("\t","-",EvidenceText)) %>%
  dplyr::mutate(ReferenceOrTrialID=gsub("\t","-",ReferenceOrTrialID))
outF <-  paste0(figDir,"/civic_MOA_clinically_actionable_list.txt")
write.table(dbRules,outF,row.names=F,quote=F,sep="\t")

### summarize
dbGenome <- dbRules %>%
  dplyr::group_by(chr,pos,ref,alt) %>%
  dplyr::summarize(n.db.entries=dplyr::n(),
                   source=paste0(unique(source),collapse=";")) %>% 
  dplyr::arrange(desc(n.db.entries))
dbGenome <- dbGenome[!dbGenome$ref=="" & !dbGenome$alt=="",]
outF <-  paste0(figDir,"/civic_MOA_clinically_actionable_genomic_pos.txt")
write.table(dbGenome,outF,row.names=F,quote=F,sep="\t")

dbAlteration <- dbRules %>%
  dplyr::group_by(gene,AAChange) %>%
  dplyr::summarize(n.db.entries=dplyr::n(),
                   source=paste0(unique(source),collapse=";"),
                   diseases=paste0(unique(Disease),collapse=";"),
                   actionability.summary=paste0(unique(actionability.summary),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";")) %>% 
  dplyr::arrange(desc(n.db.entries))
dbAlteration <- dbAlteration[!is.na(dbAlteration$AAChange),]
print(table(dbAlteration$source))

# find explicit AA changes
iAACoord <- nchar(dbAlteration$AAChange) <= 8 
excludeList <- c("Loss", "LOH", "LOSS", "ITD", "VIII", "BCR-ABL", "DELETION", "EML4-ALK", "Fusion", "NPM-ALK", "Mutation", "MUTATION","AGK-BRAF")
iAACoord <- iAACoord & !(dbAlteration$AAChange %in% excludeList)
dbProtein <- dbAlteration[iAACoord,]
dbOtherAlt <- dbAlteration[!iAACoord,]
# output non-match events
outF <-  paste0(figDir,"/civic_MOA_non_AA_change_alterations.txt")
write.table(dbOtherAlt,outF,row.names=F,quote=F,sep="\t")

table(dbProtein[,"source"]) # except 
table(dbOtherAlt[,"source"])
