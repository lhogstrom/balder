library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)

# load variant DB

figDir <- "../../output/actionability_db_curration_20231220"

### load clinical evidence entries
inFile <- "../../data/CIViC/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"
clinical <- read.csv(inFile,sep="\t")
clinical$chromosomeNum <- as.character(clinical$chromosome)
clinical$chr <- paste0("chr",clinical$chromosomeNum)

### Paths to other CIViC data (not currently used)
# inFile <- "../../data/CIViC/CIViC-01-Dec-2021-VariantSummaries.tsv"
# varSummaries <- read.csv(inFile,sep="\t") %>%
#   dplyr::arrange(desc(civic_variant_evidence_score)) # load and order by variant score
# 
# inFile <- "../../data/CIViC//CIViC-01-Dec-2021-AssertionSummaries.tsv"
# assertions <- read.csv(inFile,sep="\t")
# 
# inFile <- "../../data/CIViC//CIViC-01-Dec-2021-VariantGroupSummaries.tsv"
# varGroups <- read.csv(inFile,sep="\t")

# COMMAND ----------
### load MOAlmanac which is similar to CIViC
inFile <- "../../data/MOA/MOAlmanac_43018_2021_243_MOESM2_ESM.txt"
moa <- read.csv(inFile,sep="\t")

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

### which are the most common cancer types in msk, moa, and civic
civic.cType.cnt <- clinical %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(number.of.subjects=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.subjects))

moa.cType.cnt <- moa %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(number.of.subjects=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.subjects))

### Use MSK-IMPACT as cancer type anchor 

inFile <- "../../data/MSK_IMPACT/impact_2017_annotated_per_variant.tsv"
msk <- read.csv(inFile,sep="\t")
inFile <- "../../data/MSK_IMPACT/msk_impact_2017/data_clinical_sample.txt"
patient.msk <- read.csv(inFile,sep="\t",skip = 4)
msk <- msk %>% 
  dplyr::left_join(patient.msk,by=c("Tumor_Sample_Barcode"="SAMPLE_ID"))

msk.cType.cnt <- msk %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(number.of.subjects=dplyr::n_distinct(PATIENT_ID)) %>%
  dplyr::arrange(desc(number.of.subjects))
outF <-  paste0(figDir,"/msk_cancer_types.txt")
write.table(msk.cType.cnt,outF,row.names=F,quote=F,sep="\t")

### create cancer type maps CIVIC-->MSK and MOA-->MSK

## create draft cancer type maps
iCMmatch <- civic.cType.cnt$Disease %in% msk.cType.cnt$CANCER_TYPE
CivicMSKCancerTypeMap <- data.frame(CivicCancerType=civic.cType.cnt$Disease)
rownames(CivicMSKCancerTypeMap) <- CivicMSKCancerTypeMap$CivicCancerType
CivicMSKCancerTypeMap[iCMmatch,"MskCancerType"] <- as.character(civic.cType.cnt$Disease[iCMmatch])
mapF <-  paste0(figDir,"/civic_cancer_type_map.txt")
write.table(CivicMSKCancerTypeMap,mapF,row.names=F,quote=F,sep="\t")

### load refined map
#mapF <- paste0(figDir,"/civic_to_msk_cancer_type_map.tsv")
mapF <- "../../data/curration/civic_cancer_type_map.txt"
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
mapF <- "../../data/curration/moa_cancer_type_map.tsv"
MoaMSKCancerTypeMap <- read.csv(mapF,sep="\t")
MoaMSKCancerTypeMap$MoaCancerType <- as.character(MoaMSKCancerTypeMap$MoaCancerType)
moa <- moa %>%
  dplyr::left_join(MoaMSKCancerTypeMap,by=c("disease"="MoaCancerType"))
table(moa$MskCancerType,exclude=NULL)
# what fraction of entries are still missing MSK cancer type assignment
print(sum(is.na(moa$MskCancerType))/dim(moa)[[1]])

#####################
### MOA and CIVIC ###
#####################

matchCols <- c("source","gene","AAChange","Drugs","FDAApproved","ReferenceOrTrialID","EvidenceText","Phase",
               "Indicated","Disease","chr","pos","ref","alt","MskCancerType","actionability.summary","clinical.evidence.summary")
dbRules <- rbind(moa[,matchCols],clinical[,matchCols]) %>%
  dplyr::mutate(EvidenceText=gsub("\t","-",EvidenceText)) %>%
  dplyr::mutate(ReferenceOrTrialID=gsub("\t","-",ReferenceOrTrialID))
outF <-  paste0(figDir,"/civic_MOA_clinically_actionable_list.txt")
write.table(dbRules,outF,row.names=F,quote=F,sep="\t")

### group by genomic position and AA change
dbGenome <- dbRules %>%
  dplyr::group_by(chr,pos,ref,alt) %>%
  dplyr::summarize(n.db.entries=dplyr::n(),
                   source=paste0(unique(source),collapse=";")) %>% 
  dplyr::arrange(desc(n.db.entries))
dbGenome <- dbGenome[!dbGenome$ref=="" & !dbGenome$alt=="",]
outF <-  paste0(figDir,"/civic_MOA_clinically_actionable_genomic_pos.txt")
write.table(dbGenome,outF,row.names=F,quote=F,sep="\t")

## group by AA change
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
dbAlteration$isHGVSpChange <- iAACoord

# output non-match events
outF <-  paste0(figDir,"/civic_MOA_non_AA_change_alterations.txt")
write.table(dbOtherAlt,outF,row.names=F,quote=F,sep="\t")

########################
### write results db ###
########################

bDir <- "../../data/processed/balderResultsDb"
timestamp <- format(Sys.time(), "%Y%m%d")

# output db for harmonized data
outDbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v",timestamp,".sqlite")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), outDbName)
# output db for compiled raw data
outDbName2 <- paste0(bDir,"/balder-compiled-raw-data-v",timestamp,".sqlite")
rawDataDb <- DBI::dbConnect(RSQLite::SQLite(), outDbName2)

RSQLite::dbWriteTable(harmonizedDb, "MoaCiVICRuleEntries", dbRules)
RSQLite::dbWriteTable(harmonizedDb, "actionableSNVsByGenomicCoordinate", dbGenome)
RSQLite::dbWriteTable(harmonizedDb, "actionableSNVsByAAChange", dbAlteration)
RSQLite::dbWriteTable(harmonizedDb, "MoaCiVICOtherAlterations", dbOtherAlt)

######################################################
### patient and sample tables from various studies ###
######################################################

### load hartwig data
baseDir <- "../../data/Hartwig/data" 
inFile <- paste0(baseDir,"/samples.txt")
sampleList <- read.csv(inFile,sep=",", header = F)

inFile <- paste0(baseDir,"/metadata.tsv")
hMeta <- read.csv(inFile,sep="\t")

inFile <- paste0(baseDir,"/pre_biopsy_drugs.tsv")
hTreatment <- read.csv(inFile,sep="\t")

RSQLite::dbWriteTable(rawDataDb, "hartwigSampleInfo", sampleList)
RSQLite::dbWriteTable(rawDataDb, "hartwigMetadata", hMeta)
RSQLite::dbWriteTable(rawDataDb, "hartwigPreBiopsyDrugs", hTreatment)

# sample columns needed: cancer type, met/primary
hSampleSubset <- hMeta[,c("sampleId","primaryTumorLocation","primaryTumorType")]
colnames(hSampleSubset) <- c("SAMPLE_ID","CANCER_TYPE","CANCER_TYPE_DETAILED")
hSampleSubset$ONCOTREE_CODE <- NA # no code listed by authors
hSampleSubset$SAMPLE_TYPE <- "Metastasis" # assume all samples are met from Hartwig
hSampleSubset$SourceStudy <- "Hartwig-data"
RSQLite::dbWriteTable(harmonizedDb, "harmonizedSampleInfo", hSampleSubset,overwrite=T)

###  MSK-IMPACT
inFile <- "../../data/MSK_IMPACT/msk_impact_data_clinical_sample2.txt"
patient.msk <- read.csv(inFile,sep="\t")
RSQLite::dbWriteTable(rawDataDb, "mskMetadata", patient.msk)

mskSampleSubset <- patient.msk[,c("SAMPLE_ID","CANCER_TYPE","CANCER_TYPE_DETAILED","ONCOTREE_CODE","SAMPLE_TYPE")]
mskSampleSubset$SourceStudy <- "MSK-IMPACT-2017-data"
RSQLite::dbWriteTable(harmonizedDb, "harmonizedSampleInfo", mskSampleSubset,append=T)

### PANCAN
inFile <- "../../data/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_clinical_sample.txt"
pancanSampInfo <- read.csv(inFile,sep="\t",skip = 4)
RSQLite::dbWriteTable(rawDataDb, "pancanMetadata", pancanSampInfo)

pancanSampleSubset <- pancanSampInfo[,c("SAMPLE_ID","CANCER_TYPE","CANCER_TYPE_DETAILED","ONCOTREE_CODE","SAMPLE_TYPE")]
pancanSampleSubset$SourceStudy <- "PANCAN-WGS-data"
RSQLite::dbWriteTable(harmonizedDb, "harmonizedSampleInfo", pancanSampleSubset,append=T)

### AACR GENIE
inFile <- "../../data/AACR_Project_GENIE/Release_14p1_public/data_clinical_patient.txt"
patient.genie <- read.csv(inFile,sep="\t",skip=4)
RSQLite::dbWriteTable(rawDataDb, "GeniePatientData", patient.genie)

### detailed column definitions for patient data
# gpatientCols <- c("Patient.Identifier", 
#                   "Sex", 
#                   "Primary.Race", 
#                   "Ethnicity.Category", 
#                   "Center", 
#                   "Interval.in.days.from.DOB.to.date.of.last.contact",
#                   "Interval.in.days.from.DOB.to.DOD", 
#                   "Year.of.last.contact",                
#                   "Vital.Status",                            
#                   "Year.of.death")

inFile <- "../../data/AACR_Project_GENIE/Release_14p1_public/data_clinical_sample.txt"
sample.genie <- read.csv(inFile,sep="\t",skip=4)
RSQLite::dbWriteTable(rawDataDb, "GenieClinicalSampleData", sample.genie)

genieSampleSubset <- sample.genie[,c("SAMPLE_ID","CANCER_TYPE","CANCER_TYPE_DETAILED","ONCOTREE_CODE","SAMPLE_TYPE")]
genieSampleSubset$SourceStudy <- "AACR-GENIE-data"
RSQLite::dbWriteTable(harmonizedDb, "harmonizedSampleInfo", genieSampleSubset,append=T)

# to-do: compile minimal patient columns across data sets: age, gender, race, ethnicity 

### MC3 TCGA
# TCGA site codes
inFile <- "../../data/curration/TCGA_tissue_source_site_codes.csv"
tss <- read.csv(inFile,sep=",")
RSQLite::dbWriteTable(harmonizedDb, "Mc3TCGASiteCodes", tss)
# barcode naming convention: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
# tissue site codes: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tissue-source-site-codes

## Review harmonized output
#sRes <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM harmonizedSampleInfo')
#print(dim(sRes))

###########################################
### short variants - SNV and indel data ###
###########################################

# for each dataset: 
#     1) write the full table (*PatientVarients)
#     2) a subset of columns across all datasets (patientObservedVariantTable)

### pancan
#inFile <- "../../data/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_cna.txt"
#cna.full <- read.csv(inFile,sep="\t")
#
#inFile <- "../../data/ICGC_TCGA_WGS_2020/CNA_Genes.txt"
#cyto <- read.csv(inFile,sep="\t")
#
inFile <- "../../data/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_mutations.txt"
svPancan <- read.csv(inFile,sep="\t",skip = 2)
RSQLite::dbWriteTable(harmonizedDb, "PancanPatientVarients", svPancan)

#subset variant columns for variant table
varTable <- svPancan[,c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")]
colnames(varTable) <- c("chrom","pos","ref","alt","sample")
varTable$SourceStudy <- "PANCAN-WGS-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTable)

### MSK-IMPACT - see above for initial data loading
msk$Chromosome.str <- paste0("chr",msk$Chromosome)
msk$MAF <- 100*(msk$t_alt_count / (msk$t_ref_count+msk$t_alt_count))
RSQLite::dbWriteTable(harmonizedDb, "Msk2017PatientVarients", msk)

varTableMsk <- msk[,c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")]
colnames(varTableMsk) <- c("chrom","pos","ref","alt","sample")
varTableMsk$SourceStudy <- "MSK-IMPACT-2017-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableMsk,append=T)

### MC3 TCGA
inFile <- "../../data/mc3_tcga/scratch.sample.mc3.maf"
#inFile <- "../../data/mc3_tcga/mc3.v0.2.8.PUBLIC.maf"
mc3 <- read.csv(inFile,sep="\t") %>%
  dplyr::rename(StrandPlusMinus=STRAND)
RSQLite::dbWriteTable(harmonizedDb, "Mc3PatientVarients", mc3)

### consequence df
# selectColM3 <- c("Hugo_Symbol",
#                  "Variant_Classification",
#                  "Tumor_Sample_Barcode",
#                  "HGVSp_Short",
#                  "Transcript_ID")
#conseq.mc3 <- mc3[,selectColM3]

varTableMC3 <- mc3[,c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")]
colnames(varTableMC3) <- c("chrom","pos","ref","alt","sample")
varTableMC3$SourceStudy <- "TCGA-MC3-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableMC3,append=T)

### to-do: check to see if previous MSK-IMPACT data set is a subset of the newest GENIE dataset
### add label for assay used - what does each panel capture? 

### AACR GENIE
inFile <- "../../data/AACR_Project_GENIE/Release_14p1_public/data_mutations_extended.txt"
genie <- read.csv(inFile,sep="\t")
RSQLite::dbWriteTable(harmonizedDb, "GeniePatientVarients", genie)

varTableGenie <- genie[,c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode")]
colnames(varTableGenie) <- c("chrom","pos","ref","alt","sample")
varTableGenie$SourceStudy <- "AACR-GENIE-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableGenie,append=T)
# vRes <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM patientObservedVariantTable')
# print(dim(vRes))

### Disconnect from SQL db
RSQLite::dbDisconnect(harmonizedDb)

########
outRFile <- paste0(figDir,"/actionDbCurration20240227.RData")
save.image(file = outRFile)

### to-do
# columns needed: Hugo_Symbol, AA change, cancer type
# sample columns needed: cancer type, met/primary
# patient columns needed: age, gender, race, ethnicity 

