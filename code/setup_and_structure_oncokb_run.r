library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)
library(ggplot2)

### connect to result DB and get variants
outDir <- "../../output/clinical_annotation_matching_20240503"
bDir <- "../../data/processed/balderResultsDb"
dbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v20240412.sqlite")
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

# Load list of oncokb genes
inFile <- "../../data/oncokb/cancerGeneList.tsv"
genesOncokb <- read.csv(inFile,sep="\t")

# filter variants to oncoKB gene list
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
outF <- paste0(outDir,"/compiled_mutations_column_subset_all_studies.txt")
write.table(svCompiledPrelim[,oncoKBCols],outF,row.names=F,quote=F,sep="\t")

##########################
### load ONCOKB output ###
##########################

# Define the order of levels
level_order <- c("LEVEL_R1", "LEVEL_1", "LEVEL_2", "LEVEL_Dx1", "LEVEL_Px1", 
                 "LEVEL_3A", "LEVEL_3B", "LEVEL_4", "LEVEL_R2", "LEVEL_Dx2", 
                 "LEVEL_Px2", "LEVEL_Dx3", "LEVEL_Px3")

# Create a function to determine the highest level
get_highest_level <- function(levels) {
  levels <- levels[levels %in% level_order]
  if(length(levels) == 0) return(NA)
  return(levels[order(match(levels, level_order))][1])
}

### load OncoKB annotator results - create summary columns
inFile <- paste0(outDir,"/compiled_mutations_column_subset_all_studies.oncokb.txt")
oncokbRes <- read.csv(inFile,sep="\t") %>%
  dplyr::mutate(IS_SENS_DX_PX_OR_RESISTANCE=!(HIGHEST_SENSITIVE_LEVEL == "") | 
                  !(HIGHEST_RESISTANCE_LEVEL == "") | 
                  !(HIGHEST_DX_LEVEL == "") | 
                  !(HIGHEST_PX_LEVEL == ""),
                HAS_SENSITIVE_ENTRY=!(HIGHEST_SENSITIVE_LEVEL == ""),
                HAS_RESISTANCE_ENTRY=!(HIGHEST_RESISTANCE_LEVEL == ""),
                HAS_DX_ENTRY=!(HIGHEST_DX_LEVEL == ""),
                HAS_PX_ENTRY=!(HIGHEST_PX_LEVEL == ""),
                HIGHEST_SENS_DX_PX_OR_RESISTANCE = paste(HIGHEST_SENSITIVE_LEVEL, 
                                                         HIGHEST_DX_LEVEL, 
                                                         HIGHEST_PX_LEVEL, 
                                                         HIGHEST_RESISTANCE_LEVEL, recycle0 = T)) %>%
  rowwise() %>%
  mutate(HIGHEST_LEVEL_SUMMARY = get_highest_level(c(HIGHEST_DX_LEVEL, HIGHEST_PX_LEVEL, HIGHEST_LEVEL))) %>%
  ungroup()

RSQLite::dbWriteTable(harmonizedDb, "oncoKBAnnotatedVariantTable", oncokbRes, overwrite=T)

# confirm oncokb entries match order of annotated variants
table(oncokbRes$Tumor_Sample_Barcode == svCompiledPrelim$Tumor_Sample_Barcode,exclude=NULL)
table(oncokbRes$HGVSp == svCompiledPrelim$HGVSp,exclude=NULL)
table(oncokbRes$ONCOTREE_CODE == svCompiledPrelim$ONCOTREE_CODE,exclude=NULL)

### summary of oncokb results
okbVariantSummary <- oncokbRes[oncokbRes$IS_SENS_DX_PX_OR_RESISTANCE==T,] %>%
  dplyr::group_by(Hugo_Symbol,HGVSp_Short,ONCOTREE_CODE) %>%
  dplyr::summarise(n=dplyr::n(),
                   sensitivity_entry_cnt=sum(HAS_SENSITIVE_ENTRY),
                   resistance_entry_cnt=sum(HAS_RESISTANCE_ENTRY),
                   DX_entry_cnt=sum(HAS_DX_ENTRY),
                   Px_entry_cnt=sum(HAS_PX_ENTRY)) %>%
  dplyr::arrange(desc(n))
outF <- paste0(outDir,"/oncokb_variant_summary_table.txt")
write.table(okbVariantSummary,outF,row.names=F,quote=F,sep="\t")


