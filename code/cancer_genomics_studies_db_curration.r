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
# define output columns for harmonized sample info
harmSampleCols <- c("SAMPLE_ID",
                    "CANCER_TYPE",
                    "CANCER_TYPE_DETAILED",
                    "ONCOTREE_CODE",
                    "SAMPLE_TYPE",
                    "SAMPLE_TYPE_DETAILED",
                    "SEQ_ASSAY_ID",
                    "STAGE",
                    "TUMOR_PURITY",
                    "INT_Biopsy_To_Death",
                    "AGE_AT_SEQ_REPORT",
                    "DEAD",
                    "INT_SEQ_TO_DEATH",
                    "INT_SEQ_TO_CONTACT",
                    "SourceStudy")

### load hartwig data
baseDir <- "../../data/Hartwig/data" 
inFile <- paste0(baseDir,"/samples.txt")
sampleList <- read.csv(inFile,sep=",", header = F)

inFile <- paste0(baseDir,"/metadata.tsv")
hMeta <- read.csv(inFile,sep="\t")

inFile <- paste0(baseDir,"/pre_biopsy_drugs.tsv")
hTreatment <- read.csv(inFile,sep="\t")

RSQLite::dbWriteTable(rawDataDb, "hartwigSampleInfo", sampleList, overwrite = T)
RSQLite::dbWriteTable(rawDataDb, "hartwigMetadata", hMeta, overwrite = T)
RSQLite::dbWriteTable(rawDataDb, "hartwigPreBiopsyDrugs", hTreatment, overwrite = T)

### calculate difference of age and sample type
biopsyDateNum <- as.Date(hMeta$biopsyDate,"%Y-%m-%d")
deathDateNum <- as.Date(hMeta$deathDate,"%Y-%m-%d")
hMeta$INT_Biopsy_To_Death <- deathDateNum - biopsyDateNum

# sample columns needed: cancer type, met/primary
hSampleSubset <- hMeta[,c("sampleId","primaryTumorLocation","primaryTumorType","tumorPurity","INT_Biopsy_To_Death")]
colnames(hSampleSubset) <- c("SAMPLE_ID","CANCER_TYPE","CANCER_TYPE_DETAILED","tumorPurity","INT_Biopsy_To_Death")
hSampleSubset$TUMOR_PURITY <- hSampleSubset$tumorPurity*0.01
hSampleSubset$SAMPLE_TYPE <- "Metastasis" # assume all samples are met from Hartwig
hSampleSubset$SourceStudy <- "Hartwig-data"
hSampleSubset$ONCOTREE_CODE <- NA # no code listed by authors
hSampleSubset$SAMPLE_TYPE_DETAILED <- NA
hSampleSubset$SEQ_ASSAY_ID <- NA
hSampleSubset$AGE_AT_SEQ_REPORT <- NA
hSampleSubset$DEAD <- NA
hSampleSubset$STAGE <- NA
hSampleSubset$INT_SEQ_TO_DEATH <- NA
hSampleSubset$INT_SEQ_TO_CONTACT <- NA

RSQLite::dbWriteTable(harmonizedDb, "harmonizedSampleInfo", hSampleSubset[,harmSampleCols],overwrite=T)

###  MSK-IMPACT
inFile <- "../../data/MSK_IMPACT/msk_impact_data_clinical_sample2.txt"
patient.msk <- read.csv(inFile,sep="\t")
RSQLite::dbWriteTable(rawDataDb, "mskMetadata", patient.msk, overwrite=T)

mskSampleSubset <- patient.msk[,c("SAMPLE_ID","CANCER_TYPE","CANCER_TYPE_DETAILED","ONCOTREE_CODE","SAMPLE_TYPE","TUMOR_PURITY")]
mskSampleSubset$SourceStudy <- "MSK-IMPACT-2017-data"
mskSampleSubset$TUMOR_PURITY <- mskSampleSubset$TUMOR_PURITY*0.01
mskSampleSubset$SAMPLE_TYPE_DETAILED <- NA
mskSampleSubset$SEQ_ASSAY_ID <- NA
mskSampleSubset$AGE_AT_SEQ_REPORT <- NA
mskSampleSubset$STAGE <- NA
mskSampleSubset$DEAD <- NA
mskSampleSubset$INT_Biopsy_To_Death <- NA
mskSampleSubset$INT_SEQ_TO_DEATH <- NA
mskSampleSubset$INT_SEQ_TO_CONTACT <- NA

RSQLite::dbWriteTable(harmonizedDb, "harmonizedSampleInfo", mskSampleSubset[,harmSampleCols],append=T)

### PANCAN
inFile <- "../../data/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_clinical_sample.txt"
pancanSampInfo <- read.csv(inFile,sep="\t",skip = 4)
RSQLite::dbWriteTable(rawDataDb, "pancanMetadata", pancanSampInfo, overwrite=T)

pancanSampleSubset <- pancanSampInfo[,c("SAMPLE_ID","CANCER_TYPE","CANCER_TYPE_DETAILED","ONCOTREE_CODE","SAMPLE_TYPE","STAGE","PURITY")]
pancanSampleSubset$SourceStudy <- "PANCAN-WGS-data"
pancanSampleSubset$TUMOR_PURITY <- pancanSampleSubset$PURITY
pancanSampleSubset$SAMPLE_TYPE_DETAILED <- NA
pancanSampleSubset$SEQ_ASSAY_ID <- NA
pancanSampleSubset$AGE_AT_SEQ_REPORT <- NA
pancanSampleSubset$DEAD <- NA
pancanSampleSubset$INT_Biopsy_To_Death <- NA
pancanSampleSubset$INT_SEQ_TO_DEATH <- NA
pancanSampleSubset$INT_SEQ_TO_CONTACT <- NA

RSQLite::dbWriteTable(harmonizedDb, "harmonizedSampleInfo", pancanSampleSubset[,harmSampleCols],append=T)

### AACR GENIE
inFile <- "../../data/AACR_Project_GENIE/Release_14p1_public/data_clinical_patient.txt"
patient.genie <- read.csv(inFile,sep="\t",skip=4)
RSQLite::dbWriteTable(rawDataDb, "GeniePatientData", patient.genie, overwrite=T)

inFile <- "../../data/AACR_Project_GENIE/Release_14p1_public/data_clinical_sample.txt"
sample.genie <- read.csv(inFile,sep="\t",skip=4)
RSQLite::dbWriteTable(rawDataDb, "GenieClinicalSampleData", sample.genie, overwrite=T)

### join sample data to patient data and create inferred variables
sample.genie$AGE_AT_SEQ_REPORT_NUMERIC <- as.numeric(sample.genie$AGE_AT_SEQ_REPORT)
sample.genie$INF_AGE_SEQ_REPORT_DAYS <- (sample.genie$AGE_AT_SEQ_REPORT_NUMERIC * 365) + 183
samplePatient <- sample.genie %>%
  dplyr::left_join(patient.genie,by="PATIENT_ID") %>%
  dplyr::mutate(
    INT_SEQ_TO_DEATH = if_else(
      !is.na(as.numeric(INT_DOD)) & !is.na(as.numeric(INF_AGE_SEQ_REPORT_DAYS)),
      as.numeric(INT_DOD) - INF_AGE_SEQ_REPORT_DAYS,
      NA  # Assign NA of type double if conditions are not met
    ),
    INT_SEQ_TO_CONTACT = if_else(
      !is.na(as.numeric(INT_CONTACT)) & !is.na(as.numeric(INF_AGE_SEQ_REPORT_DAYS)),
      as.numeric(INT_CONTACT) - INF_AGE_SEQ_REPORT_DAYS,
      NA  # Assign NA of type double if conditions are not met
    )
  )

genieSampleSubset <- samplePatient[,c("SAMPLE_ID",
                                      "CANCER_TYPE",
                                      "CANCER_TYPE_DETAILED",
                                      "ONCOTREE_CODE",
                                      "SAMPLE_TYPE",
                                      "SAMPLE_TYPE_DETAILED",
                                      "SEQ_ASSAY_ID",
                                      "AGE_AT_SEQ_REPORT",
                                      "DEAD",
                                      "INT_SEQ_TO_DEATH",
                                      "INT_SEQ_TO_CONTACT")]
genieSampleSubset$SourceStudy <- "AACR-GENIE-data"
genieSampleSubset$STAGE <- NA
genieSampleSubset$TUMOR_PURITY <- NA
genieSampleSubset$INT_Biopsy_To_Death <- NA

RSQLite::dbWriteTable(harmonizedDb, "harmonizedSampleInfo", genieSampleSubset[,harmSampleCols],append=T)

# to-do: compile minimal patient columns across data sets: age, gender, race, ethnicity 

### MC3 TCGA
# TCGA site codes
inFile <- "../../data/curration/TCGA_tissue_source_site_codes.csv"
tss <- read.csv(inFile,sep=",")
RSQLite::dbWriteTable(harmonizedDb, "Mc3TCGASiteCodes", tss, overwrite=T)
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
RSQLite::dbWriteTable(rawDataDb, "PancanPatientVarients", svPancan, overwrite=T)
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
RSQLite::dbWriteTable(rawDataDb, "Msk2017PatientVarients", msk, overwrite=T)

### MC3 TCGA
inFile <- "../../data/mc3_tcga/scratch.sample.mc3.maf"
#inFile <- "../../data/mc3_tcga/mc3.v0.2.8.PUBLIC.maf"
mc3 <- read.csv(inFile,sep="\t") %>%
  dplyr::rename(StrandPlusMinus=STRAND)
RSQLite::dbWriteTable(rawDataDb, "Mc3PatientVarients", mc3, overwrite=T)

### to-do: check to see if previous MSK-IMPACT data set is a subset of the newest GENIE dataset
### add label for assay used - what does each panel capture? 

### AACR GENIE
inFile <- "../../data/AACR_Project_GENIE/Release_14p1_public/data_mutations_extended.txt"
genie <- read.csv(inFile,sep="\t")
RSQLite::dbWriteTable(rawDataDb, "GeniePatientVarients", genie, overwrite=T)

### Hartwig
inFile <- "../../data/Hartwig/data/output_v4/hartwig_clinically_actionable_pcgr_entries.txt"
hartwig.pcgr <- read.csv(inFile,sep="\t")
RSQLite::dbWriteTable(rawDataDb, "HartwigPatientVarients", hartwig.pcgr, overwrite=T)


harwigMod <- hartwig.pcgr %>%
  dplyr::mutate(Hugo_Symbol = SYMBOL,
                Entrez_Gene_Id = ENTREZ_ID,
                Center = NA, 
                NCBI_Build = NA, 
                Chromosome = CHROM, 
                Start_Position = POS, 
                End_Position = NA, 
                Strand = NA, 
                Variant_Classification = NA, 
                Variant_Type = NA, 
                Reference_Allele = REF, 
                Tumor_Seq_Allele1 = NA, 
                Tumor_Seq_Allele2 = ALT, 
                dbSNP_RS = NA, 
                dbSNP_Val_Status = NA, 
                Tumor_Sample_Barcode = sample, 
                Matched_Norm_Sample_Barcode = NA, 
                Match_Norm_Seq_Allele1 = NA, 
                Match_Norm_Seq_Allele2 = NA, 
                Tumor_Validation_Allele1 = NA, 
                Tumor_Validation_Allele2 = NA, 
                Match_Norm_Validation_Allele1 = NA, 
                Match_Norm_Validation_Allele2 = NA, 
                Verification_Status = NA, 
                Validation_Status = NA, 
                Mutation_Status = NA, 
                Sequencing_Phase = NA, 
                Sequence_Source = NA, 
                Validation_Method = NA, 
                Score = NA, 
                BAM_File = NA, 
                Sequencer = NA, 
                t_ref_count = NA, 
                t_alt_count = NA, 
                n_ref_count = NA, 
                n_alt_count = NA, 
                HGVSp_Short = HGVSp_short, 
                Transcript_ID = NA)

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
varTable$TVAF <- NA
varTable$SourceStudy <- "PANCAN-WGS-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTable, overwrite=T)

varTableMsk <- msk[,uCols]
varTableMsk$TVAF <- NA
varTableMsk$SourceStudy <- "MSK-IMPACT-2017-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableMsk,append=T)

varTableMC3 <- mc3[,uCols]
varTableMC3$TVAF <- NA
varTableMC3$SourceStudy <- "TCGA-MC3-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableMC3,append=T)

varTableGenie <- genie[,uCols]
varTableGenie$TVAF <- NA
varTableGenie$SourceStudy <- "AACR-GENIE-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableGenie,append=T)

varTableHartwig <- harwigMod[,c(uCols,"TVAF")]
varTableHartwig$SourceStudy <- "Hartwig-data"
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTable", varTableHartwig,append=T)

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
