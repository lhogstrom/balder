library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)
#source("R/snv_indel_annotation.R")

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 2) {
  stop("Two arguments must be supplied", call. = FALSE)
}

# Assign Arguments to Variables
baseDir <- args[1] # should be the relative path to "balder/code" where this file is located
dbName <- args[2]

utilFile <- paste0(baseDir,"/../code/R/snv_indel_annotation.R")
source(utilFile)

### connect to result DB and get variants
outDir <- paste0(baseDir,"/../../output/actionability_db_curration_20240712")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), dbName)


### connect to result DB and get variants
#outDir <- "../../output/clinical_annotation_matching_20240503"
#bDir <- "../../data/processed/balderResultsDb"
#dbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v20240412.sqlite")
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
                MAF=100*(t_alt_count / (t_ref_count+t_alt_count)),
                t_alt_plus_ref_count = t_ref_count+t_alt_count)
poVariants$balderVariantID <- seq(1,dim(poVariants)[[1]])
# replace MAF with TVAF for hartwig variants
iHartwig <- poVariants$SourceStudy == "Hartwig-data"
poVariants[iHartwig,"MAF"] <- poVariants[iHartwig,"TVAF"]*100

# load sample info for observed variants
sampleInfoCompiled <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM harmonizedSampleInfo')
sampleInfoCompiled[is.na(sampleInfoCompiled$SAMPLE_TYPE),"sampleInfoCompiled"] <- "Unspecified"
#compile hartwig and genie seq or report to death
sampleInfoCompiled$INT_Seq_or_Biopsy_To_Death <- sampleInfoCompiled$INT_Biopsy_To_Death
iNonNullSeqDeath <- !is.na(sampleInfoCompiled[,"INT_SEQ_TO_DEATH"])
sampleInfoCompiled[iNonNullSeqDeath,"INT_Seq_or_Biopsy_To_Death"] <- sampleInfoCompiled[iNonNullSeqDeath,"INT_SEQ_TO_DEATH"]

# add oncoKB results
oncokbRes <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM oncoKBAnnotatedVariantTable')
colnames(oncokbRes) <- paste0("ONCOKB_",colnames(oncokbRes))

svCompiled <- combine_patient_table_and_oncokb(poVariants,
                                               sampleInfoCompiled,
                                               ot_code_full,
                                               oncokbRes)
RSQLite::dbWriteTable(harmonizedDb, "patientObservedVariantTableWithSampleInfoOncokb", svCompiled,overwrite=T)

##########################################################
### exhaustive and representative annotation selection ###
##########################################################

aa.genome.exahustive <- annotation_matching_genomic_coordinate_and_AA_change(svCompiled,dbRules)
RSQLite::dbWriteTable(harmonizedDb, "exhaustiveClinicalAnnotatedPatientVariants", aa.genome.exahustive,overwrite=T)

# pick a single representative variant per subject
aa.genome.representative <- aa.genome.exahustive %>%
  dplyr::arrange(desc(cancerTypeMatch),
                 clinical.evidence.summary) %>% # prioritize restricted cancer type matches that are high-confidence
  dplyr::group_by(Tumor_Sample_Barcode,balderVariantID) %>%
  dplyr::filter(dplyr::row_number()==1) %>%
  dplyr::mutate(clinical_annotation_status="clinically annotated")
write.csv(aa.genome.representative,
          paste0(outDir,"/actionability_representative_variant_table.csv"),
          row.names = F)
# SQL write
RSQLite::dbWriteTable(harmonizedDb, "representativeClinicalAnnotatedPatientVariants", aa.genome.representative,overwrite=T)
