library(DBI)
library(RSQLite)
library(phenOncoX)
library(dplyr)
library(UpSetR)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

## Goal ##
# The goal of this script is to summarize cancer type matching approaches 

figDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/cancer_type_matching"

bDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/resultsDb"
mydb <- DBI::dbConnect(RSQLite::SQLite(), paste0(bDir,"/actionable-biomarker-db.sqlite"))
RSQLite::dbListTables(mydb)

#### test cancer type mappings ###
# inF <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/oncotree_tumor_types_2021.txt"
# oncotree <- read.csv(inF,sep="\t")
# l1Onco <- unique(oncotree$level_1)

########################################
## load data from compiled databases ##
########################################

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
outFile <- paste0(figDir,"/source_cancer_type_terms.txt")
write.table(sourceCancerTypes,outFile,sep="\t",row.names = F)

####################
## phenOncoX load ##
####################

### load phenOncoX
download_dir <- tempdir()
oncoterms <- phenOncoX::get_terms(
  cache_dir = download_dir)
outFile <- paste0(figDir,"/phenoOncoX_records.txt")
write.table(oncoterms$records,outFile,sep="\t",row.names = F)

##########################################################
## get representative phenOncoX ID for each source term ##
##########################################################

#### perform queries 
# define name columns
nameCols <- c("primary_site",
              "ot_main_type",
              "ot_name",
              "ot_code",
              "cui_name",
              "efo_name",
              "do_name")

stringTypeList <- c("origSourceString","modifiedSourceString")
recordMatchList <- c("exactStringMatch","substringMatch")

dfRepresenativeSets <- data.frame()
for (stringType in stringTypeList) {
  for (recordMatch in recordMatchList) {
    for (sourceTermOrig in sourceCancerTypes$CANCER_TYPE) {
      print(sourceTermOrig)
      
      # define and modify source term
      sourceTerm <- tolower(sourceTermOrig)
      sourceTerm <- gsub(" cancer", "", sourceTerm)
      sourceTerm <- gsub(" disease", "", sourceTerm) # other terms to consider: “recurrent” or “malignant” or “adult”
      
      ### Create matrix of boolean matches
      oncotermNames <- oncoterms$records[,nameCols]
      queryMatch <- oncotermNames
      # for each column convert to lower
      for (xCol in nameCols) {
        resRecordsOrig <- oncotermNames[,xCol]
        resRecords <- tolower(resRecordsOrig)
        resRecords <- gsub("_", " ", resRecords) 
        # remove phrases
        oncotermNames[,xCol] <- resRecords
        
        if (stringType == "origSourceString" & recordMatch == "substringMatch") {
          queryMatch[,xCol] <- grepl(sourceTermOrig,resRecordsOrig) # is it a substring
        }
        
        if (stringType == "modifiedSourceString" & recordMatch == "substringMatch") {
          queryMatch[,xCol] <- grepl(sourceTerm,resRecords) # is it a substring
        }
 
        if (stringType == "origSourceString" & recordMatch == "exactStringMatch") {
          queryMatch[,xCol] <- sourceTermOrig == resRecordsOrig # require exact match
        }
        
        if (stringType == "modifiedSourceString" & recordMatch == "exactStringMatch") {
          queryMatch[,xCol] <- sourceTerm == resRecords # require exact match
        }
        
        # NOTE: what is the highest resoultion you can have with each database? List average/distirbuion OT level for each DB.
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
      
      #outFile <- paste0(figDir,"/oncoterm_query_result_tmp.txt")
      #write.table(oncoQueryRes,outFile,sep="\t",row.names = F)
      
      oncoQueryRes$queryHitRepresentative <- F
      repEntryMaxHit <- oncoQueryRes %>% 
        dplyr::filter(queryHits > 0) %>%
        dplyr::filter(queryHits == max(queryHits)) %>%
        dplyr::filter(dplyr::row_number()==1) %>%
        mutate(queryHitRepresentative=T,
               searchMatchType="MaxHit",
               stringTypeApproach=stringType,
               recordMatchApproach=recordMatch)
      dfRepresenativeSets <- rbind(dfRepresenativeSets,repEntryMaxHit)
      
      repEntryMaxLevel <- oncoQueryRes %>% 
        dplyr::filter(queryHits > 0) %>%
        dplyr::arrange(desc(ot_level)) %>%
        dplyr::filter(dplyr::row_number()==1) %>%
        mutate(queryHitRepresentative=T,
               searchMatchType="MaxLevel",
               stringTypeApproach=stringType,
               recordMatchApproach=recordMatch)
      dfRepresenativeSets <- rbind(dfRepresenativeSets,repEntryMaxLevel)
    }
  }
}

outFile <- paste0(figDir,"/oncoterm_representative_query_results_full.txt")
write.table(dfRepresenativeSets,outFile,sep="\t",row.names = F)

table(dfRepresenativeSets$searchMatchType,dfRepresenativeSets$stringTypeApproach,dfRepresenativeSets$recordMatchApproach)

###
dfRepresenativeMaxHit <- dfRepresenativeSets[dfRepresenativeSets$searchMatchType=="MaxHit" & 
                                             dfRepresenativeSets$stringTypeApproach=="modifiedSourceString" &
                                             dfRepresenativeSets$recordMatchApproach=="substringMatch",]

dfRepresenativeMaxLevel <- dfRepresenativeSets[dfRepresenativeSets$searchMatchType=="MaxLevel" & 
                                                 dfRepresenativeSets$stringTypeApproach=="modifiedSourceString" &
                                                 dfRepresenativeSets$recordMatchApproach=="substringMatch",]

outFile <- paste0(figDir,"/oncoterm_representative_query_result_max_hits.txt")
write.table(dfRepresenativeMaxHit,outFile,sep="\t",row.names = F)

outFile <- paste0(figDir,"/oncoterm_representative_query_result_max_oncotree_level.txt")
write.table(dfRepresenativeMaxLevel,outFile,sep="\t",row.names = F)

# which terms had no matches? 
noMatchTerms <- sourceCancerTypes[!sourceCancerTypes$CANCER_TYPE %in% dfRepresenativeMaxLevel$queryTerm,]
outFile <- paste0(figDir,"/cancer_type_term_no_phenoncox_hit.txt")
write.table(noMatchTerms,outFile,sep="\t",row.names = F)

### summary stats for Max Hit vs Max Level
table(dfRepresenativeMaxHit$ot_level,exclude=NULL)
table(dfRepresenativeMaxLevel$ot_level,exclude=NULL)

table(dfRepresenativeMaxHit$queryHits,exclude=NULL)
table(dfRepresenativeMaxLevel$queryHits,exclude=NULL)

### plot the number of query hits, plot the oncotree higherarchy number
outF <-  paste0(figDir,"/oncoterm_query_matches_max_hit.png")
ggplot(dfRepresenativeMaxHit,aes(x=queryHits))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(Hugo_Symbol~direction,scale="free_y",)+
  xlab("PhenoOncoX query mathces")+
  ggtitle(paste0("Largest number of column matches for \n query hits"))+
  #scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggsave(outF,height = 5,width = 5)

outF <-  paste0(figDir,"/oncoterm_query_matches_max_level.png")
ggplot(dfRepresenativeMaxLevel,aes(x=queryHits))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(Hugo_Symbol~direction,scale="free_y",)+
  xlab("PhenoOncoX query mathces")+
  ggtitle(paste0("Largest number of column matches for \n query hits"))+
  #scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
#theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggsave(outF,height = 5,width = 5)

### summarise which columns drove the hits
dCols <- colnames(dfRepresenativeMaxLevel)
querMatchCols <- grepl("query_",dCols)
d1 <- data.frame(colSums(dfRepresenativeMaxLevel[,querMatchCols]))
d1$repType <- "MaxLevel"
d1$ColumnMatch <- rownames(d1)
colnames(d1) <- c("count","repType","ColumnMatch")

d2 <- data.frame(colSums(dfRepresenativeMaxHit[,querMatchCols]))
d2$repType <- "MaxHit"
d2$ColumnMatch <- rownames(d1)
colnames(d2) <- c("count","repType","ColumnMatch")

hitCnts <- rbind(d1,d2)

outF <-  paste0(figDir,"/oncoterm_query_matches_rep_compare.png")
ggplot(hitCnts,aes(x=ColumnMatch,y=count,color=repType))+
  geom_point()+
  #geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(Hugo_Symbol~direction,scale="free_y",)+
  xlab("PhenoOncoX query mathces")+
  ggtitle(paste0("Sum of matches across source terms"))+
  scale_color_brewer(palette="Accent",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(outF,height = 5,width = 5)

# what proportion of source and match terms were unique
repDf <- dfRepresenativeMaxHit
#repDf <- dfRepresenativeMaxLevel
nHits <- dim(repDf)[[1]]
nUniqueHits <- length(unique(repDf$recordID))
nQueryTerms <- length(unique(repDf$queryTerm))
print("Proportion of unique query terms:")
print(nQueryTerms/nHits)
print("Proportion of unique record terms (after match):")
print(nUniqueHits/nHits)
print("Proportion of increased matches:")
print((nQueryTerms-nUniqueHits)/nQueryTerms)

#### of those match hits, how many were in each respective source DB?
uniqueMatchTerms <- repDf %>%
  dplyr::group_by(queryTerm) %>%
  dplyr::filter(dplyr::row_number() == 1)
sourceMatch <- sourceCancerTypes %>%
  dplyr::left_join(uniqueMatchTerms,by=c("CANCER_TYPE"="queryTerm")) %>%
  data.frame()

sourceMatch$ot_level_num <- as.numeric(sourceMatch$ot_level)
sourceMatch[is.na(sourceMatch$ot_level_num),"ot_level_num"] <- 0

### output source matches
outFile <- paste0(figDir,"/phenoOncoX_record_ID_matches_all_sources_v2.csv")
write.table(sourceMatch,outFile,sep=",",row.names = F)

############################################
## compare terms from source DBs pairwise ##
############################################

## Make pairwise comparisons of term overlap proprtions in source dbs
dbSources <- unique(sourceCancerTypes$source)
pairwiseMatchDf <- data.frame()
for (xSource in dbSources) {
  for (ySource in dbSources) {
    pairString <- paste0(xSource, " & ", ySource)
    print(pairString)
    pairwiseMatchDf[pairString,"source1"] <- xSource
    pairwiseMatchDf[pairString,"source2"] <- ySource
    # dfs for each source
    xI <- sourceMatch$source.x == xSource
    yI <- sourceMatch$source.x == ySource
    pairwiseMatchDf[pairString,"source1_count"] <- sum(xI)
    pairwiseMatchDf[pairString,"source2_count"] <- sum(yI)
    
    # overlap in source terms
    xSourceTerm <- sourceMatch[xI,"CANCER_TYPE"]
    ySourceTerm <- sourceMatch[yI,"CANCER_TYPE"]
    pairwiseMatchDf[pairString,"orig_terms_from_source1_in_source2"] <- sum(xSourceTerm %in% ySourceTerm)
    print(table(xSourceTerm %in% ySourceTerm))
    
    # overlap on record IDs
    xRecordID <- c(sourceMatch[xI,"recordID"])
    yRecordID <- c(sourceMatch[yI,"recordID"])
    print(table(xRecordID %in% yRecordID))
    pairwiseMatchDf[pairString,"phenoOncoX_records_from_source1_in_source2"] <- sum(xRecordID %in% yRecordID)
    
    # no overlap in pheno onco records
    NoMatch <- xRecordID[!(xRecordID %in% yRecordID)]
    iNoMatch <- sourceMatch$recordID %in% NoMatch
    noMatchDf <- sourceMatch[xI & iNoMatch,]
    outFile <- paste0(figDir,"/tmp_out_phenoOncoX_record_ID_",xSource,"_not_found_in_",ySource,".csv")
    write.table(noMatchDf,outFile,sep=",",row.names = F)
    
    # perc of total
    pairwiseMatchDf[pairString,"perc_orig_terms_from_source1_in_source2"] <- 100*(sum(xSourceTerm %in% ySourceTerm) / sum(xI))
    pairwiseMatchDf[pairString,"perc_phenoOncoX_records_from_source1_in_source2"] <- 100*(sum(xRecordID %in% yRecordID) / sum(xI))
    
    pairwiseMatchDf[pairString,"perc_match_improvement"] <- pairwiseMatchDf[pairString,"perc_phenoOncoX_records_from_source1_in_source2"] - pairwiseMatchDf[pairString,"perc_orig_terms_from_source1_in_source2"]
  }
}

outFile <- paste0(figDir,"/phenoOncoX_record_ID_term_match_data.csv")
write.table(pairwiseMatchDf,outFile,sep=",",row.names = F)

sourceTermDfMtrx <- matrix(pairwiseMatchDf$perc_match_improvement, nrow = 4, ncol = 4, byrow = TRUE,
       dimnames = list(dbSources,
                       dbSources))

heatmapFile=paste0(figDir,"/phenoOncoX_record_ID_term_match_heatmap.pdf")
gcdfOut <- pheatmap(sourceTermDfMtrx,#as.matrix(t(sourceTermDfMtrx)),#annotation = msk.protein.cnt.df[,c("Hugo_Symbol","protein_change")],#annotation_colors =ann_colors,
                    show_rownames = T,#labels_row = rownames(raw.prot),
                    show_colnames = T,#labels_col = meta[gnums,"name"],
                    color = rev(brewer.pal(8,"GnBu")),
                    display_numbers = T,
                    fontsize_col = 9,
                    fontsize_row = 9,
                    cluster_cols=F,cluster_rows=F,
                    filename = heatmapFile,width = 8, height = 8)

## how often are the query terms different for the same matching records? 
uniqueQueryTerms <- repDf %>%
  dplyr::group_by(recordID) %>%
  dplyr::summarise(nQueryTerms = dplyr::n_distinct(queryTerm))
print(table(uniqueQueryTerms$nQueryTerms))

## show distribution of 
#dbSummary <- sourceMatch %>%
#  dplyr::group_by(source.x) %>%

outF <-  paste0(figDir,"/oncoterm_oncotree_level_by_db.png")
ggplot(sourceMatch,aes(x=ot_level_num))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  facet_grid(source.x~.)+
  #facet_grid(source.x~.,scale="free_y",)+
  xlab("Oncotree level")+
  ggtitle(paste0("Cancer type ontology level"))+
  #scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
#theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggsave(outF,height = 9,width = 5)
## to-do: covert to fraction of counts (not raw counts). Split out CIViC and MOA. Include CGI actionability dataset.



### Parse “ot_code_path” field to get parent connections of the oncotree levels
# t <- "LUNG-LNET-ALUCA"
# ts <- strsplit(t,"-")
# ts[[1]][1]
# ts[[1]][2]
# ts[[1]][3]


### to-do for sigve:
# 2) list out intersection of terms between src1 and src2
# 3) What is the set of pheonooncox matches for "skin Melanoma" before picking representative entry


#### 
### load MSK-IMPACT
inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/impact_2017_annotated_per_variant.tsv"
msk <- read.csv(inFile,sep="\t")
#inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/msk_impact_data_clinical_sample2.txt"
#patient.msk <- read.csv(inFile,sep="\t")
inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/msk_impact_2017/data_clinical_sample.txt"
patient.msk <- read.csv(inFile,sep="\t",skip = 4)
msk$Chromosome.str <- paste0("chr",msk$Chromosome)
msk$MAF <- 100*(msk$t_alt_count / (msk$t_ref_count+msk$t_alt_count))

### Is there a primary and secondary alt listed? 
#table(is.na(msk$n_alt_count),exclude=NULL)
#table(msk$Reference_Allele == msk$Tumor_Seq_Allele1,exclude=NULL)

# join clinical info onto variant info
all(msk$Tumor_Sample_Barcode %in% patient.msk$SAMPLE_ID)

# join patient info and represenative phenOncoX record
msk <- msk %>% 
  dplyr::left_join(patient.msk,by=c("Tumor_Sample_Barcode"="SAMPLE_ID")) %>%
  dplyr::left_join(dfRepresenativeMaxLevel[,c("queryTerm","recordID")],by=c("CANCER_TYPE"="queryTerm"))

table(is.na(msk$recordID))

####




# load tables:
#pcgrAnnotation <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM pcgrAnnotation') # issues loading this data product


# pcgrAnnotation
# actionableSNVsByGenomicCoordinate
# patientObservedVariantTable

dbAlteration <- actionRules %>%
  dplyr::left_join(dfRepresenativeMaxLevel[,c("queryTerm","recordID")],by=c("CANCER_TYPE"="queryTerm"))
table(is.na(dbAlteration$recordID))

### join with MOA - by protein change
dbAlteration$MOAset <- "Y"
# msk.protein <- msk %>%
#   dplyr::left_join(dbAlteration,by=c("Hugo_Symbol"="gene",
#                                      "protein_change"="AAChange"))
# print(table(msk.protein$MOAset,msk.protein$source,exclude = NULL))
# print(table(msk.protein$feature_type,exclude = NULL))

msk.protein.cnt <- msk %>%
  dplyr::inner_join(dbAlteration,by=c("Hugo_Symbol"="gene",
                                      "protein_change"="AAChange")) %>%
  dplyr::group_by(Hugo_Symbol,protein_change) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID),
                   source=paste0(unique(source),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(msk.protein.cnt))
print(table(msk.protein.cnt$source,exclude = NULL))

### impose cancer type matches (protein) - orig string
msk.protein.cancer <- msk %>%
  dplyr::inner_join(dbAlteration,by=c("Hugo_Symbol"="gene",
                                 "protein_change"="AAChange",
                                 "CANCER_TYPE"="MskCancerType"))
msk.protein.cancer.cnt <- msk.protein.cancer %>%
  dplyr::group_by(Hugo_Symbol,protein_change,CANCER_TYPE) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(msk.protein.cancer.cnt))
print(table(msk.protein.cancer.cnt$source,exclude = NULL))
print(sum(msk.protein.cancer.cnt$n.patients.mutated))


### impose cancer type matches (protein) - phenoOncoX recrod
msk.protein.cancer <- msk %>%
  dplyr::inner_join(dbAlteration,by=c("Hugo_Symbol"="gene",
                                      "protein_change"="AAChange",
                                      "recordID"="recordID"))
msk.protein.cancer.cnt <- msk.protein.cancer %>%
  dplyr::group_by(Hugo_Symbol,protein_change,recordID) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(recordID),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(msk.protein.cancer.cnt))
print(table(msk.protein.cancer.cnt$source,exclude = NULL))
print(sum(msk.protein.cancer.cnt$n.patients.mutated))


########
outRFile <- paste0(figDir,"/cancerTypeMatching20231113.RData")
save.image(file = outRFile)

