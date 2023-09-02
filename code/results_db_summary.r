library(DBI)
library(RSQLite)
library(phenOncoX)
library(dplyr)
library(UpSetR)

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

#RSQLite::dbDisconnect(mydb)

#### test cancer type mappings ###
inF <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/oncotree_tumor_types_2021.txt"
oncotree <- read.csv(inF,sep="\t")
l1Onco <- unique(oncotree$level_1)

## load current MOA maps ##
inF <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Pancan_biomarker_incidence_v2/moa_cancer_type_map.tsv"
moaMap <- read.csv(inF,sep="\t")

### search level 1 Oncotree terms for one MOA term
#iX <-1
#xMoa <- moaMap[iX,"MoaCancerType"]
#grepl(xMoa,l1Onco)
#l1Onco[grepl("Lung",l1Onco)]

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
sourceCancerTypes$CANCER_TYPE_lower_nc <- gsub(" cancer","",sourceCancerTypes$CANCER_TYPE_lower)

### compare to phenOncoX
download_dir <- tempdir()
phenTreeDepth2 <- phenOncoX::get_tree(
  cache_dir = download_dir, max_tree_depth = 2)
head(phenTreeDepth2$records)

oncoterms <- phenOncoX::get_terms(
  cache_dir = download_dir)
outFile <- paste0(figDir,"/phenoOncoX_records.txt")
write.table(oncoterms$records,outFile,sep="\t",row.names = F)

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
  dplyr::mutate(oncoterms_primary_site = (CANCER_TYPE_lower %in% tolower(oncoterms$records$primary_site)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$primary_site)),
                oncoterms_ot_main_type = (CANCER_TYPE_lower %in% tolower(oncoterms$records$ot_main_type)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$ot_main_type)),
                oncoterms_ot_name = (CANCER_TYPE_lower %in% tolower(oncoterms$records$ot_name)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$ot_name)),
                oncoterms_ot_code = (CANCER_TYPE_lower %in% tolower(oncoterms$records$ot_code)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$ot_code)),
                oncoterms_cui_name = (CANCER_TYPE_lower %in% tolower(oncoterms$records$cui_name)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$cui_name)),
                oncoterms_efo_name = (CANCER_TYPE_lower %in% tolower(oncoterms$records$efo_name)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$efo_name)),
                oncoterms_do_name = (CANCER_TYPE_lower %in% tolower(oncoterms$records$do_name)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$do_name)),
                oncoterm_match_any = oncoterms_primary_site | oncoterms_ot_main_type | oncoterms_ot_name | oncoterms_ot_code | oncoterms_cui_name | oncoterms_efo_name | oncoterms_do_name)

### create exhaustive matching between target name (from source database) and phenoOncoX oncoterms
tmp.records$CANCER_TYPE <- ctLower
match.df <- data.frame()
oncoterms$records$recordID <- seq(1,dim(oncoterms$records)[[1]])
for (iX in rownames(sourceCancerTypes)) {
  ctLower <- as.character(sourceCancerTypes[iX,"CANCER_TYPE_lower"])
  ctLowerNC <- as.character(sourceCancerTypes[iX,"CANCER_TYPE_lower_nc"])
  xSource <- as.character(sourceCancerTypes[iX,"source"])
  nEntries <- as.character(sourceCancerTypes[iX,"nEntries"])
  i_primary_site <- which((ctLower == tolower(oncoterms$records$primary_site)) | (ctLowerNC == tolower(oncoterms$records$primary_site)))
  i_ot_main_type <- which((ctLower == tolower(oncoterms$records$ot_main_type)) | (ctLowerNC == tolower(oncoterms$records$ot_main_type)))
  i_ot_name <- which((ctLower == tolower(oncoterms$records$ot_name)) | (ctLowerNC == tolower(oncoterms$records$ot_name)))
  i_ot_code <- which((ctLower == tolower(oncoterms$records$ot_code)) | (ctLowerNC == tolower(oncoterms$records$ot_code)))
  i_cui_name <- which((ctLower == tolower(oncoterms$records$cui_name)) | (ctLowerNC == tolower(oncoterms$records$cui_name)))
  i_efo_name <- which((ctLower == tolower(oncoterms$records$efo_name)) | (ctLowerNC == tolower(oncoterms$records$efo_name)))
  i_do_name <- which((ctLower == tolower(oncoterms$records$do_name)) | (ctLowerNC == tolower(oncoterms$records$do_name)))
  iRecords <- unique(c(i_primary_site,i_ot_main_type,i_ot_name,i_ot_code,i_cui_name,i_efo_name,i_do_name))
  print(length(iRecords))
  tmp.records <- oncoterms$records[iRecords,]
  if (dim(tmp.records)[[1]]>0) {
    tmp.records$CANCER_TYPE <- ctLower
    tmp.records$source<- xSource
    tmp.records$nEntries <- nEntries
    match.df <- rbind(match.df,tmp.records)
  }
}

match.select <- match.df %>%
  dplyr::group_by(CANCER_TYPE) %>% 
  dplyr::filter(ot_level==min(ot_level))
outFile <- paste0(figDir,"/phenoOncoX_database_select_matches.txt")
write.table(match.select,outFile,sep="\t",row.names = F)


# sourceIndex <- sourceCancerTypes %>%
#   dplyr::mutate(oncoterms_primary_site = toString(which((CANCER_TYPE_lower == tolower(oncoterms$records$primary_site)) | (CANCER_TYPE_lower_nc == tolower(oncoterms$records$primary_site)))),
#                 oncoterms_ot_main_type = toString(which((CANCER_TYPE_lower == tolower(oncoterms$records$ot_main_type)) | (CANCER_TYPE_lower_nc == tolower(oncoterms$records$ot_main_type)))),
#                 oncoterms_ot_name = toString(which((CANCER_TYPE_lower == tolower(oncoterms$records$ot_name)) | (CANCER_TYPE_lower_nc == tolower(oncoterms$records$ot_name)))))#,
                #oncoterms_ot_code = (CANCER_TYPE_lower %in% tolower(oncoterms$records$ot_code)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$ot_code)),
                #oncoterms_cui_name = (CANCER_TYPE_lower %in% tolower(oncoterms$records$cui_name)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$cui_name)),
                #oncoterms_efo_name = (CANCER_TYPE_lower %in% tolower(oncoterms$records$efo_name)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$efo_name)),
                #oncoterms_do_name = (CANCER_TYPE_lower %in% tolower(oncoterms$records$do_name)) | (CANCER_TYPE_lower_nc %in% tolower(oncoterms$records$do_name)),
                #oncoterm_match_any = oncoterms_primary_site | oncoterms_ot_main_type | oncoterms_ot_name | oncoterms_ot_code | oncoterms_cui_name | oncoterms_efo_name | oncoterms_do_name)

outFile <- paste0(figDir,"/oncoterm_mapping_booleans_v2.txt")
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
                   oncoterms_do_name = 100*sum(oncoterms_do_name)/srcCnt,
                   oncoterms_any_term = 100*sum(oncoterm_match_any)/srcCnt)
outFile <- paste0(figDir,"/oncoterm_mapping_proportions_v2.txt")
write.table(srcTblSum,outFile,sep="\t",row.names = F)


# define and modify source term
sourceTermOrig <- "Adrenal Gland Angiosarcoma"
#sourceTermOrig <- "adrenal gland"
sourceTerm <- tolower(sourceTermOrig)
sourceTerm <- gsub(" cancer", "", sourceTerm)
sourceTerm <- gsub(" disease", "", sourceTerm)

### Create matrix of boolean matches

# define name columns
nameCols <- c("primary_site",
              "ot_main_type",
              "ot_name",
              "ot_code",
              "cui_name",
              "efo_name",
              "do_name")
oncotermNames <- oncoterms$records[,nameCols]
queryMatch <- oncotermNames
# for each column convert to lower
for (xCol in nameCols) {
  resRecords <- tolower(oncotermNames[,xCol])
  resRecords <- gsub("_", " ", resRecords)
  # remove phrases
  oncotermNames[,xCol] <- resRecords
  queryMatch[,xCol] <- grepl(sourceTerm,resRecords)
}

queryCols <- paste0("query_match_",nameCols)
colnames(queryMatch) <- queryCols

# provide index
oncoterms$records$recordID <- seq(1,dim(oncoterms$records)[[1]])
queryMatch$recordID <- seq(1,dim(queryMatch)[[1]])

# join terms with query results
oncoQueryRes <- oncoterms$records %>%
  dplyr::left_join(queryMatch,by="recordID")
oncoQueryRes$queryHits <- rowSums(oncoQueryRes[,queryCols])
oncoQueryRes$queryTerm <- sourceTermOrig

outFile <- paste0(figDir,"/oncoterm_query_result_tmp.txt")
write.table(oncoQueryRes,outFile,sep="\t",row.names = F)


## select representative term
# by max query hit count

# by largest oncotree level - with arbitrary (do term selected)




#' Search take an input string and search exhaustivly across oncoterms
#'
#' @param sourceTerm 
#' @param searchLower 
#'
#' @return
#' @export
#'
#' @examples
getOncotermMatches <- function(sourceTerm, searchLower=T) {
  return(x)
}


### edit distance-based matches
#eDistDf <- adist(sourceCancerTypes$CANCER_TYPE,oncoterms$records[,c("cui_name","efo_name")])


### Upset plot: how many cancer types appear in each data set 
# movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
#                    header = T, sep = ";")
# upset(movies, nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2, 
#       mainbar.y.label = "Genre Intersections", sets.x.label = "Movies Per Genre", 
#       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))
# upset(movies, sets = c("Action", "Adventure", "Comedy"), mb.ratio = c(0.55, 0.45), order.by = "freq")




###To do
# 1. check influence of allowing one edit-distance match
# 2. Rank according to most common exception that doesnt match. Is there a string like "cancer" that should be excluded? 
# 3. Review phenOncoX codetools
# 4. what proportion of pancan cancers are in MOA - create exhaustive venn
# 4. Plot influence of exact vs. generalized cancer type matching


