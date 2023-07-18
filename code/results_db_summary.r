library(DBI)
library(RSQLite)
library(phenOncoX)
library(dplyr)

## Goal ##
# The goal of this script is to summarize the contents of the results db

figDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers"

bDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/resultsDb"
mydb <- DBI::dbConnect(RSQLite::SQLite(), paste0(bDir,"/actionable-biomarker-db.sqlite"))
RSQLite::dbListTables(mydb)

# show metadata table
subsetHmeta <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM hartwigMetadata LIMIT 5')
print(head(subsetHmeta))
poVariantTable <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM patientObservedVariantTable LIMIT 5')
print(head(poVariantTable))


# view a head and dim of each table
tableList <- RSQLite::dbListTables(mydb)
for (tbl in tableList) {
  print(tbl)
  sQuery <- paste0('SELECT * FROM ',tbl)
  tblDf <- RSQLite::dbGetQuery(mydb, sQuery)
  print(dim(tblDf))
  nColMax <- min(15,dim(tblDf)[[2]])
  print(head(tblDf[1:5,colnames(tblDf)[1:nColMax]]))
}

RSQLite::dbDisconnect(mydb)

#### test cancer type mappings ###
inF <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/oncotree_tumor_types_2021.txt"
oncotree <- read.csv(inF,sep="\t")
l1Onco <- unique(oncotree$level_1)

## load current MOA maps ##
inF <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Pancan_biomarker_incidence_v2/moa_cancer_type_map.tsv"
moaMap <- read.csv(inF,sep="\t")

iX <-1
xMoa <- moaMap[iX,"MoaCancerType"]

grepl(xMoa,l1Onco)
l1Onco[grepl("Lung",l1Onco)]

### get metadata tables
actionRules <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM MoaCiVICRuleEntries')
table(actionRules$Disease)
actionRules$CANCER_TYPE <- actionRules$Disease
actionCTypes <- actionRules %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(nEntries=n())
actionCTypes$source <- "MoaCiVICRuleEntries"

mskMetadata <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM mskMetadata')
table(mskMetadata$CANCER_TYPE)
mskCTypes <- mskMetadata %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(nEntries=n())
mskCTypes$source <- "mskMetadata"

pancanMetadata <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM pancanMetadata')
table(pancanMetadata$CANCER_TYPE)
pancanCTypes <- pancanMetadata %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(nEntries=n())
pancanCTypes$source <- "pancanMetadata"

hartwigMetadata <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM hartwigMetadata')
table(hartwigMetadata$CANCER_TYPE)
hartwigMetadata$CANCER_TYPE <- hartwigMetadata$primaryTumorLocation
hartwigCTypes <- hartwigMetadata %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(nEntries=n())
hartwigCTypes$source <- "hartwigMetadata"

sourceCancerTypes <- rbind(actionCTypes,mskCTypes,pancanCTypes,hartwigCTypes)
sourceCancerTypes$CANCER_TYPE_lower <- tolower(sourceCancerTypes$CANCER_TYPE)

### compare to phenOncoX
download_dir <- tempdir()
phenTreeDepth2 <- phenOncoX::get_tree(
  cache_dir = download_dir, max_tree_depth = 2)

oncoterms <- phenOncoX::get_terms(
  cache_dir = download_dir)

head(phenTreeDepth2$records)

### what proportion of names are at different levels? 
# sourceCancerTypesEval <- sourceCancerTypes %>%
#   dplyr::mutate(oncoterms_primary_site = CANCER_TYPE %in% oncoterms$records$primary_site,
#                 oncoterms_ot_main_type = CANCER_TYPE %in% oncoterms$records$ot_main_type,
#                 oncoterms_ot_name = CANCER_TYPE %in% oncoterms$records$ot_name,
#                 oncoterms_ot_code = CANCER_TYPE %in% oncoterms$records$ot_code,
#                 oncoterms_cui_name = CANCER_TYPE %in% oncoterms$records$cui_name,
#                 oncoterms_efo_name = CANCER_TYPE %in% oncoterms$records$efo_name,
#                 oncoterms_do_name = CANCER_TYPE %in% oncoterms$records$do_name)
#                 #phenTreeDepth2_ = CANCER_TYPE %in% phenTreeDepth2$records$,

sourceCancerTypesEval <- sourceCancerTypes %>%
  dplyr::mutate(oncoterms_primary_site = CANCER_TYPE_lower %in% tolower(oncoterms$records$primary_site),
                oncoterms_ot_main_type = CANCER_TYPE_lower %in% tolower(oncoterms$records$ot_main_type),
                oncoterms_ot_name = CANCER_TYPE_lower %in% tolower(oncoterms$records$ot_name),
                oncoterms_ot_code = CANCER_TYPE_lower %in% tolower(oncoterms$records$ot_code),
                oncoterms_cui_name = CANCER_TYPE_lower %in% tolower(oncoterms$records$cui_name),
                oncoterms_efo_name = CANCER_TYPE_lower %in% tolower(oncoterms$records$efo_name),
                oncoterms_do_name = CANCER_TYPE_lower %in% tolower(oncoterms$records$do_name))
#phenTreeDepth2_ = CANCER_TYPE %in% phenTreeDepth2$records$,

outFile <- paste0(figDir,"/oncoterm_mapping_booleans.txt")
write.table(sourceCancerTypesEval,outFile,sep="\t",row.names = F)

srcTblSum <- sourceCancerTypesEval %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(srcCnt = dplyr::n_distinct(CANCER_TYPE),
                   oncoterms_primary_site = 100*sum(oncoterms_primary_site)/srcCnt,
                   oncoterms_ot_main_type = 100*sum(oncoterms_ot_main_type)/srcCnt,
                   oncoterms_ot_name = 100*sum(oncoterms_ot_name)/srcCnt,
                   oncoterms_ot_code = 100*sum(oncoterms_ot_code)/srcCnt,
                   oncoterms_cui_name = 100*sum(oncoterms_cui_name)/srcCnt,
                   oncoterms_efo_name = 100*sum(oncoterms_efo_name)/srcCnt,
                   oncoterms_do_name = 100*sum(oncoterms_do_name)/srcCnt)
outFile <- paste0(figDir,"/oncoterm_mapping_proportions.txt")
write.table(srcTblSum,outFile,sep="\t",row.names = F)
