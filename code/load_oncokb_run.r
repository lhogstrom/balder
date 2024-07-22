library(dplyr)
library(DBI)
library(RSQLite)

args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 2) {
  stop("Two arguments must be supplied", call. = FALSE)
}

# Assign Arguments to Variables
baseDir <- args[1] # should be the relative path to "balder/code" where this file is located
dbName <- args[2]

### connect to result DB and get variants
outDir <- paste0(baseDir,"/../../output/actionability_db_curration_20240712")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), dbName)

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
#table(oncokbRes$Tumor_Sample_Barcode == svCompiledPrelim$Tumor_Sample_Barcode,exclude=NULL)
#table(oncokbRes$HGVSp == svCompiledPrelim$HGVSp,exclude=NULL)
#table(oncokbRes$ONCOTREE_CODE == svCompiledPrelim$ONCOTREE_CODE,exclude=NULL)

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


