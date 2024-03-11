library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)

outDir <- "../../output/actionability_db_curration_20231220"

################################
### establish db connection ###
################################

bDir <- "../../data/processed/balderResultsDb"
timestamp <- format(Sys.time(), "%Y%m%d")

# output db for harmonized data
outDbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v",timestamp,".sqlite")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), outDbName)
# output db for compiled raw data
outDbName2 <- paste0(bDir,"/balder-compiled-raw-data-v",timestamp,".sqlite")
rawDataDb <- DBI::dbConnect(RSQLite::SQLite(), outDbName2)

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

################################################
### short variants - load SNV and indel data ###
################################################

### pancan
inFile <- "../../data/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_mutations.txt"
svPancan <- read.csv(inFile,sep="\t",skip = 2)
RSQLite::dbWriteTable(rawDataDb, "PancanPatientVarients", svPancan)
#inFile <- "../../data/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_cna.txt"
#cna.full <- read.csv(inFile,sep="\t")
#
#inFile <- "../../data/ICGC_TCGA_WGS_2020/CNA_Genes.txt"
#cyto <- read.csv(inFile,sep="\t")
#

### MSK-IMPACT 
inFile <- "../../data/MSK_IMPACT/impact_2017_annotated_per_variant.tsv"
msk <- read.csv(inFile,sep="\t")
msk$MAF <- 100*(msk$t_alt_count / (msk$t_ref_count+msk$t_alt_count))
RSQLite::dbWriteTable(rawDataDb, "Msk2017PatientVarients", msk)

### MC3 TCGA
inFile <- "../../data/mc3_tcga/scratch.sample.mc3.maf"
#inFile <- "../../data/mc3_tcga/mc3.v0.2.8.PUBLIC.maf"
mc3 <- read.csv(inFile,sep="\t") %>%
  dplyr::rename(StrandPlusMinus=STRAND)
RSQLite::dbWriteTable(rawDataDb, "Mc3PatientVarients", mc3)

### to-do: check to see if previous MSK-IMPACT data set is a subset of the newest GENIE dataset
### add label for assay used - what does each panel capture? 

### AACR GENIE
inFile <- "../../data/AACR_Project_GENIE/Release_14p1_public/data_mutations_extended.txt"
genie <- read.csv(inFile,sep="\t")
RSQLite::dbWriteTable(rawDataDb, "GeniePatientVarients", genie)

###############################################################
### write common variant columns for each dataset to SQLite ###
###############################################################

# indentify intersection of column names across all variant sets
pancanCols <- colnames(svPancan)
genieCols <- colnames(genie)
mc3Cols <- colnames(mc3)
mskCols <- colnames(msk)
uCols <- base::intersect(base::intersect(pancanCols,genieCols),base::intersect(mc3Cols,mskCols))

varTable <- svPancan[,uCols]
varTable$SourceStudy <- "PANCAN-WGS-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTable)

varTableMsk <- msk[,uCols]
varTableMsk$SourceStudy <- "MSK-IMPACT-2017-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableMsk,append=T)

varTableMC3 <- mc3[,uCols]
varTableMC3$SourceStudy <- "TCGA-MC3-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableMC3,append=T)

varTableGenie <- genie[,uCols]
varTableGenie$SourceStudy <- "AACR-GENIE-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableGenie,append=T)

# vRes <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM patientObservedVariantTable')
# print(dim(vRes))

##############################
### Disconnect from SQL db ###
##############################

RSQLite::dbDisconnect(harmonizedDb)
RSQLite::dbDisconnect(rawDataDb)

########
outRFile <- paste0(outDir,"/genomicStudyiesDbCurration-",timestamp,".RData")
save.image(file = outRFile)
