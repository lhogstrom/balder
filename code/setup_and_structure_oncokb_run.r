library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 3) {
  stop("Three arguments must be supplied", call. = FALSE)
}

# Assign Arguments to Variables
baseDir <- args[1] # should be the relative path to "balder/code" where this file is located
dbName <- args[2]
outVarTable <- args[3]
#bDir <- "../../data/processed/balderResultsDb"/
#dbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v20240412.sqlite")

### connect to result DB and get variants
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), dbName)

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

### join sample info to variant data
svCompiledPrelim <- poVariants %>%
  dplyr::inner_join(sampleInfoCompiled[,!colnames(sampleInfoCompiled)=="SourceStudy"],
                    by=c("Tumor_Sample_Barcode"="SAMPLE_ID")) %>% # exclude samples that do not have any sample information
  dplyr::filter(!SAMPLE_TYPE == "Cell line")
print(dim(svCompiledPrelim))

## Load list of oncokb genes
#inFile <- "../../data/oncokb/cancerGeneList.tsv"
#genesOncokb <- read.csv(inFile,sep="\t")

## filter variants to oncoKB gene list
# okbCompiledSubset <- svCompiledPrelim %>%
#   dplyr::filter(Hugo_Symbol %in% genesOncokb$Hugo.Symbol) %>%
#   dplyr::group_by(Hugo_Symbol,
#                   HGVSp,
#                   HGVSp_Short,
#                   ONCOTREE_CODE) %>%
#   dplyr::summarise(n=dplyr::n()) %>%
#   dplyr::filter(!HGVSp == "") %>%
#   dplyr::arrange(desc(n))

oncoKBCols <- c("Hugo_Symbol",
                "Tumor_Sample_Barcode",
                "HGVSp",
                "HGVSp_Short",
                #"HGVSg",
                "ONCOTREE_CODE")
#outVarTable <- paste0(outDir,"/compiled_mutations_column_subset_all_studies.txt")
write.table(svCompiledPrelim[,oncoKBCols],outVarTable,row.names=F,quote=F,sep="\t")

