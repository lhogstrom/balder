library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)
library(phenOncoX)

outDir <- "../../output/actionability_db_curration_20231220"

### load clinical evidence entries
inFile <- "../../data/CIViC/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"
clinical <- read.csv(inFile,sep="\t")
clinical$chromosomeNum <- as.character(clinical$chromosome)
clinical$chr <- paste0("chr",clinical$chromosomeNum)

### Paths to other CIViC data (not currently used)
# inFile <- "../../data/CIViC/CIViC-01-Dec-2021-VariantSummaries.tsv"
# varSummaries <- read.csv(inFile,sep="\t") %>%
#   dplyr::arrange(desc(civic_variant_evidence_score)) # load and order by variant score
# 
# inFile <- "../../data/CIViC//CIViC-01-Dec-2021-AssertionSummaries.tsv"
# assertions <- read.csv(inFile,sep="\t")
# 
# inFile <- "../../data/CIViC//CIViC-01-Dec-2021-VariantGroupSummaries.tsv"
# varGroups <- read.csv(inFile,sep="\t")

# COMMAND ----------
### load MOAlmanac which is similar to CIViC
inFile <- "../../data/MOA/MOAlmanac_43018_2021_243_MOESM2_ESM.txt"
moa <- read.csv(inFile,sep="\t")

### create actionability type assignemnts 
## MOA
moa$actionability.summary <- "Other"
isPrognostic <- moa$therapy_type=="" & !is.na(moa$favorable_prognosis)
moa[isPrognostic,"actionability.summary"] <- "Prognostic"
#
isPredTherapyAct <- is.na(moa$therapy_resistance) & !moa$therapy_type==""
moa[isPredTherapyAct,"actionability.summary"] <- "Predictive therapy actionability"
#
isTherapyResistance <- !is.na(moa$therapy_resistance) & !moa$therapy_type==""
moa[isTherapyResistance,"actionability.summary"] <- "Therapy resistance"
table(moa$actionability.summary,exclude=NULL)

## CIViC
clinical$actionability.summary <- "Other"
cSigPrognostic <- c("Better Outcome","N/A","Negative","Positive","Poor Outcome")
isPrognostic <- clinical$evidence_type=="Prognostic" & (clinical$clinical_significance %in% cSigPrognostic)
clinical[isPrognostic,"actionability.summary"] <- "Prognostic"
#
isPredTherapyAct <- clinical$evidence_type=="Predictive" | (clinical$clinical_significance=="Sensitivity/Response")
clinical[isPredTherapyAct,"actionability.summary"] <- "Predictive therapy actionability"
#
isTherapyResistance <- clinical$clinical_significance=="Resistance"
clinical[isTherapyResistance,"actionability.summary"] <- "Therapy resistance"
table(clinical$actionability.summary,exclude=NULL)

### make evidence assignments 
moa$clinical.evidence.summary <- "Lower Confidence"
moa[moa$source_type == "FDA" | moa$source_type == "Guideline","clinical.evidence.summary"] <- "High Confidence"
table(moa$clinical.evidence.summary)

clinical$clinical.evidence.summary <- "Lower Confidence"
clinical[clinical$evidence_level == "A" | clinical$evidence_level == "B","clinical.evidence.summary"] <- "High Confidence"
table(clinical$clinical.evidence.summary)

##############################################################
### perform Disease Ontology to Oncotree mapping for CIViC ###
##############################################################


# load OT code hierarchy  
library(purrr)
inFile <- "../../data/oncotree/tumor_types_2021.txt"
ot_code_full <- read.csv(inFile,sep="\t") %>%
  separate(level_1, into = c("level_1_disease", "ot_code_level_1"), sep = " \\(", remove = FALSE) %>%
  separate(level_2, into = c("level_2_disease", "ot_code_level_2"), sep = " \\(", remove = FALSE) %>%
  separate(level_3, into = c("level_3_disease", "ot_code_level_3"), sep = " \\(", remove = FALSE) %>%
  separate(level_4, into = c("level_4_disease", "ot_code_level_4"), sep = " \\(", remove = FALSE) %>%
  separate(level_5, into = c("level_5_disease", "ot_code_level_5"), sep = " \\(", remove = FALSE) %>%
  separate(level_6, into = c("level_6_disease", "ot_code_level_6"), sep = " \\(", remove = FALSE) %>%
  mutate(ot_code_level_1 = gsub("\\)", "", ot_code_level_1),
         ot_code_level_2 = gsub("\\)", "", ot_code_level_2),
         ot_code_level_3 = gsub("\\)", "", ot_code_level_3),
         ot_code_level_4 = gsub("\\)", "", ot_code_level_4),
         ot_code_level_5 = gsub("\\)", "", ot_code_level_5),
         ot_code_level_6 = gsub("\\)", "", ot_code_level_6)) %>%
  mutate(Highest_Non_Null_Level = pmap_chr(select(., starts_with("ot_code_level_")), function(...) {
    # Extract the column names of non-NA values for the current row
    non_na_levels <- names(which(!is.na(c(...))))
    # If there are non-NA levels, return the highest; otherwise, return NA or a placeholder
    if(length(non_na_levels) > 0) {
      last(non_na_levels)
    } else {
      NA_character_
    }
  }))

# loop to get highest OT code and name
for (i in rownames(ot_code_full)) {
  lvl <- ot_code_full[i,"Highest_Non_Null_Level"]
  codeSelect <- ot_code_full[i,lvl]
  numLev <- strsplit(lvl, "_level_")[[1]][2]
  ot_code_full[i,"oncotree_level"] <- as.numeric(numLev)
  ot_code_full[i,"ot_code"] <- codeSelect
  nameCol <- paste0("level_",numLev,"_disease")
  ot_code_full[i,"ot_name"] <- ot_code_full[i,nameCol]
  print(codeSelect)
}

outF <-  paste0(outDir,"/ot_code_map.txt")
write.table(ot_code_full,outF,row.names=F,quote=F,sep="\t")

### load phenOncoX
download_dir <- tempdir()
oncoterms <- phenOncoX::get_terms(
  cache_dir = download_dir)

# select first unique entry with valid DO ID and OncoTree code
nonNullDoEntried <- oncoterms$records[!is.na(oncoterms$records$do_id) & !is.na(oncoterms$records$ot_code),] %>%
  dplyr::group_by(do_id) %>%
  dplyr::filter(dplyr::row_number()==1)
nonNullDoEntried$do_id_num <- as.numeric(sapply(strsplit(nonNullDoEntried$do_id, ":"), "[[", 2))

# which entries do not have a DO match in phenOncoX? 
#table(clinical$doid %in% nonNullDoEntried$do_id_num)
#clinical[!clinical$doid %in% nonNullDoEntried$do_id_num,]

### alterinative mapping tests
# nnOT <- oncoterms$records[!is.na(oncoterms$records$do_id),]
# nnOT$do_id_num <- as.numeric(sapply(strsplit(nnOT$do_id, ":"), "[[", 2))
# nnOT[nnOT$do_id_num==9256,]
# nnOT[nnOT$do_id_num==1324,]
# nnOT[nnOT$do_id_num==1612,]
# nnOT[nnOT$ot_name=="Lung" & !is.na(nnOT$ot_name),]
# nnOT[grepl("Colo",nnOT$ot_name) & !is.na(nnOT$ot_name),]
# nnOT[grepl("east",nnOT$ot_name) & !is.na(nnOT$ot_name),]

mapDf <- nonNullDoEntried[,c("ot_name","ot_code","do_id_num")]
# apply manual mapping for a few entries
mapDf <- rbind(mapDf,c("Breast","BRCA",1612))
mapDf <- rbind(mapDf,c("Colon/Rectum","BOWEL",9256))
mapDf <- rbind(mapDf,c("Lung","LUNG",1324))
mapDf$do_id_num <- as.numeric(mapDf$do_id_num)

clinicalTree <- clinical %>%
  dplyr::left_join(mapDf,
                   by=c("doid"="do_id_num"))

# check results of the join
table(is.na(clinicalTree$ot_code))
head(clinicalTree[,c("disease","doid","ot_name","ot_code")])
clinicalTree[is.na(clinicalTree$ot_code),c("disease","doid","ot_name","ot_code")]

### add matching columns
moa$source <- "MOAlmanac"
moa$FDAApproved <- NA
moa$Indicated <- moa$predictive_implication
moa$EvidenceText <- moa$description
moa$ReferenceOrTrialID <- moa$citation
moa$Disease <- moa$disease
moa$TrialID <- NA
moa$Phase <- NA
moa$clinical_significance <- NA
isFavorable <- (moa$favorable_prognosis == 1) & !is.na(moa$favorable_prognosis)
moa[isFavorable,"clinical_significance"] <- "Better Outcome"
isNotFavorable <- (moa$favorable_prognosis == 0) & !is.na(moa$favorable_prognosis)
moa[isNotFavorable,"clinical_significance"] <- "Poor Outcome"

moa$Drugs <- moa$therapy_name
moa$AberrationsLabelsMod <- NA
moa <- moa %>%
  tidyr::separate(protein_change,sep="\\.",into=c("nullP","AAChange"),remove=F)
moa$chr <- paste0("chr",moa$chromosome)
moa$pos <- moa$start_position
moa$ref <- moa$reference_allele
moa$alt <- moa$alternate_allele

clinical$source <- "CIViC-clinical"
clinical$FDAApproved <- NA
clinical$Indicated <- clinical$evidence_type
clinical$EvidenceText <- clinical$evidence_statement
clinical$ReferenceOrTrialID <- clinical$evidence_civic_url
clinical$Disease <- clinical$disease
clinical$TrialID <- NA
clinical$Phase <- NA
clinical$Drugs <- clinical$drugs
clinical$AberrationsLabelsMod <- NA
clinical$AAChange <- clinical$variant
clinical$pos <- clinical$start
clinical$ref <- clinical$reference_bases
clinical$alt <- clinical$variant_bases
clinical$oncotree_term <- clinicalTree$ot_name
clinical$oncotree_code <- clinicalTree$ot_code

### which are the most common cancer types in msk, moa, and civic
civic.cType.cnt <- clinical %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(number.of.subjects=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.subjects))

moa.cType.cnt <- moa %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(number.of.subjects=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.subjects))

### Use MSK-IMPACT as cancer type anchor 

inFile <- "../../data/MSK_IMPACT/impact_2017_annotated_per_variant.tsv"
msk <- read.csv(inFile,sep="\t")
inFile <- "../../data/MSK_IMPACT/msk_impact_2017/data_clinical_sample.txt"
patient.msk <- read.csv(inFile,sep="\t",skip = 4)
msk <- msk %>% 
  dplyr::left_join(patient.msk,by=c("Tumor_Sample_Barcode"="SAMPLE_ID"))

msk.cType.cnt <- msk %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(number.of.subjects=dplyr::n_distinct(PATIENT_ID)) %>%
  dplyr::arrange(desc(number.of.subjects))
outF <-  paste0(outDir,"/msk_cancer_types.txt")
write.table(msk.cType.cnt,outF,row.names=F,quote=F,sep="\t")

### create cancer type maps CIVIC-->MSK and MOA-->MSK

## create draft cancer type maps
iCMmatch <- civic.cType.cnt$Disease %in% msk.cType.cnt$CANCER_TYPE
CivicMSKCancerTypeMap <- data.frame(CivicCancerType=civic.cType.cnt$Disease)
rownames(CivicMSKCancerTypeMap) <- CivicMSKCancerTypeMap$CivicCancerType
CivicMSKCancerTypeMap[iCMmatch,"MskCancerType"] <- as.character(civic.cType.cnt$Disease[iCMmatch])
mapF <-  paste0(outDir,"/civic_cancer_type_map.txt")
write.table(CivicMSKCancerTypeMap,mapF,row.names=F,quote=F,sep="\t")

### load refined map
#mapF <- paste0(outDir,"/civic_to_msk_cancer_type_map.tsv")
mapF <- "../../data/curration/civic_cancer_type_map.txt"
CivicMSKCancerTypeMap <- read.csv(mapF,sep="\t")
clinical <- clinical %>%
  dplyr::left_join(CivicMSKCancerTypeMap,by=c("disease"="CivicCancerType"))
table(clinical$MskCancerType,exclude=NULL)
# what fraction of entries are still missing MSK cancer type assignment
print(sum(is.na(clinical$MskCancerType))/dim(clinical)[[1]])

### create draft map 
iCMmatch <- moa.cType.cnt$disease %in% msk.cType.cnt$CANCER_TYPE
MoaMSKCancerTypeMap <- data.frame(MoaCancerType=moa.cType.cnt$disease)
rownames(MoaMSKCancerTypeMap) <- MoaMSKCancerTypeMap$CivicCancerType
MoaMSKCancerTypeMap[iCMmatch,"MskCancerType"] <- as.character(moa.cType.cnt$disease[iCMmatch])
mapF <-  paste0(outDir,"/moa_cancer_type_map.txt")
write.table(MoaMSKCancerTypeMap,mapF,row.names=F,quote=F,sep="\t")

### load refined map
#mapF <-  paste0(outDir,"/moa_cancer_type_map.tsv")
mapF <- "../../data/curration/moa_cancer_type_map.tsv"
MoaMSKCancerTypeMap <- read.csv(mapF,sep="\t")
MoaMSKCancerTypeMap$MoaCancerType <- as.character(MoaMSKCancerTypeMap$MoaCancerType)
moa <- moa %>%
  dplyr::left_join(MoaMSKCancerTypeMap,by=c("disease"="MoaCancerType"))
table(moa$MskCancerType,exclude=NULL)
# what fraction of entries are still missing MSK cancer type assignment
print(sum(is.na(moa$MskCancerType))/dim(moa)[[1]])

##########################
### MOA and CIVIC join ###
##########################

# filter out only positive evidence entries
civicSelectedEntries <- clinical %>%
  filter(!evidence_direction=="Does Not Support",
         !clinical_significance=="N/A")

# select columns and perform join
matchCols <- c("source","gene","AAChange","Drugs","FDAApproved","ReferenceOrTrialID","EvidenceText","Phase",
               "Indicated","Disease", "oncotree_term", "oncotree_code",
               "chromosome","pos","ref","alt","MskCancerType","clinical_significance","actionability.summary","clinical.evidence.summary")
dbRules <- rbind(moa[,matchCols],civicSelectedEntries[,matchCols]) %>%
  dplyr::mutate(EvidenceText=gsub("\t","-",EvidenceText)) %>%
  dplyr::mutate(ReferenceOrTrialID=gsub("\t","-",ReferenceOrTrialID))
outF <-  paste0(outDir,"/civic_MOA_clinically_actionable_list.txt")
write.table(dbRules,outF,row.names=F,quote=F,sep="\t")

### group by genomic position and AA change
dbGenome <- dbRules %>%
  dplyr::group_by(chromosome,pos,ref,alt) %>%
  dplyr::summarize(n.db.entries=dplyr::n(),
                   source=paste0(unique(source),collapse=";")) %>% 
  dplyr::arrange(desc(n.db.entries))
dbGenome <- dbGenome[!dbGenome$ref=="" & !dbGenome$alt=="",]
outF <-  paste0(outDir,"/civic_MOA_clinically_actionable_genomic_pos.txt")
write.table(dbGenome,outF,row.names=F,quote=F,sep="\t")

## group by AA change
dbAlteration <- dbRules %>%
  dplyr::group_by(gene,AAChange) %>%
  dplyr::summarize(n.db.entries=dplyr::n(),
                   source=paste0(unique(source),collapse=";"),
                   diseases=paste0(unique(Disease),collapse=";"),
                   actionability.summary=paste0(unique(actionability.summary),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";")) %>% 
  dplyr::arrange(desc(n.db.entries))
dbAlteration <- dbAlteration[!is.na(dbAlteration$AAChange),]
print(table(dbAlteration$source))

# find explicit AA changes
iAACoord <- nchar(dbAlteration$AAChange) <= 8 
excludeList <- c("Loss", "LOH", "LOSS", "ITD", "VIII", "BCR-ABL", "DELETION", "EML4-ALK", "Fusion", "NPM-ALK", "Mutation", "MUTATION","AGK-BRAF")
iAACoord <- iAACoord & !(dbAlteration$AAChange %in% excludeList)
dbProtein <- dbAlteration[iAACoord,]
dbOtherAlt <- dbAlteration[!iAACoord,]
dbAlteration$isHGVSpChange <- iAACoord

# output non-match events
outF <-  paste0(outDir,"/civic_MOA_non_AA_change_alterations.txt")
write.table(dbOtherAlt,outF,row.names=F,quote=F,sep="\t")

########################
### write results db ###
########################

bDir <- "../../data/processed/balderResultsDb"
timestamp <- format(Sys.time(), "%Y%m%d")

# output db for harmonized data
outDbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v",timestamp,".sqlite")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), outDbName)
# output db for compiled raw data
outDbName2 <- paste0(bDir,"/balder-compiled-raw-data-v",timestamp,".sqlite")
#rawDataDb <- DBI::dbConnect(RSQLite::SQLite(), outDbName2)

# Harmonized table
RSQLite::dbWriteTable(harmonizedDb, "MoaCiVICRuleEntries", dbRules, overwrite =T)
RSQLite::dbWriteTable(harmonizedDb, "actionableSNVsByGenomicCoordinate", dbGenome, overwrite =T)
RSQLite::dbWriteTable(harmonizedDb, "actionableSNVsByAAChange", dbAlteration, overwrite =T)
RSQLite::dbWriteTable(harmonizedDb, "MoaCiVICOtherAlterations", dbOtherAlt, overwrite =T)
RSQLite::dbWriteTable(harmonizedDb, "OncoTreeCodeHierarchy", ot_code_full, overwrite =T)

### Disconnect from SQL db
RSQLite::dbDisconnect(harmonizedDb)
#RSQLite::dbDisconnect(rawDataDb)

########
outRFile <- paste0(outDir,"/actionDbCurration-",timestamp,".RData")
save.image(file = outRFile)
