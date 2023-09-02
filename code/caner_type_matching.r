library(DBI)
library(RSQLite)
library(phenOncoX)
library(dplyr)
library(UpSetR)
library(ggplot2)

## Goal ##
# The goal of this script is to summarize the contents of the results db

figDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/cancer_type_matching"

bDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/resultsDb"
mydb <- DBI::dbConnect(RSQLite::SQLite(), paste0(bDir,"/actionable-biomarker-db.sqlite"))
RSQLite::dbListTables(mydb)

#### test cancer type mappings ###
inF <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/oncotree_tumor_types_2021.txt"
oncotree <- read.csv(inF,sep="\t")
l1Onco <- unique(oncotree$level_1)

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

#### perform queries 
# define name columns
nameCols <- c("primary_site",
              "ot_main_type",
              "ot_name",
              "ot_code",
              "cui_name",
              "efo_name",
              "do_name")

dfRepresenative <- data.frame()
for (sourceTermOrig in sourceCancerTypes$CANCER_TYPE) {
  print(sourceTermOrig)
  
  # define and modify source term
  #sourceTermOrig <- "Adrenal Gland Angiosarcoma"
  #sourceTermOrig <- "adrenal gland"
  sourceTerm <- tolower(sourceTermOrig)
  sourceTerm <- gsub(" cancer", "", sourceTerm)
  sourceTerm <- gsub(" disease", "", sourceTerm)
  
  ### Create matrix of boolean matches
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
  
  oncoQueryRes$queryHitRepresentative <- F
  g <- oncoQueryRes %>% 
    dplyr::filter(queryHits == max(queryHits)) %>%
    dplyr::filter(dplyr::row_number()==1) %>%
    mutate(queryHitRepresentative=T)
  dfRepresenative <- rbind(dfRepresenative,g)

  #outFile <- paste0(figDir,"/oncoterm_query_result_tmp.txt")
  #write.table(oncoQueryRes,outFile,sep="\t",row.names = F)
  
}

outFile <- paste0(figDir,"/oncoterm_query_result_tmp2.txt")
write.table(dfRepresenative,outFile,sep="\t",row.names = F)

### plot the number of query hits, plot the oncotree higherarchy number

outF <-  paste0(figDir,"/oncoterm_query_matches.png")
ggplot(dfRepresenative,aes(x=queryHits))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(Hugo_Symbol~direction,scale="free_y",)+
  xlab("PhenoOncoX query mathces")+
  ggtitle(paste0("Largest number of column matches for \n query hits"))+
  #scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggsave(outF,height = 5,width = 5)
