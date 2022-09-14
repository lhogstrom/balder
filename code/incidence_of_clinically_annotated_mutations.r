# Databricks notebook source
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


# COMMAND ----------

# MAGIC %md 
# MAGIC # load CIViC variant knowledge base

# COMMAND ----------

figDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/mutation_incidence_20220913"

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-VariantSummaries.tsv"
# load and order by variant score
df <- read.csv(inFile,sep="\t") %>%
  dplyr::arrange(desc(civic_variant_evidence_score))

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-AssertionSummaries.tsv"
assertions <- read.csv(inFile,sep="\t")

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"
clinical <- read.csv(inFile,sep="\t")
clinical$chromosomeNum <- as.character(clinical$chromosome)
clinical$chr <- paste0("chr",clinical$chromosomeNum)

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-VariantGroupSummaries.tsv"
varGroups <- read.csv(inFile,sep="\t")

# COMMAND ----------
### load MOAlmanac which is similar to CIViC
inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/MOAlmanac_43018_2021_243_MOESM2_ESM.txt"
moa <- read.csv(inFile,sep="\t")

print(table(moa$gene))
print(table(moa$predictive_implication))
print(table(moa$source_type))
print(table(moa$feature_type))
print(table(moa$variant_annotation))
print(table(moa$cosmic_signature_number,exclude=NULL))

signatures.moa <- moa[moa$feature_type=="Mutational signature",]
head(signatures.moa)

### make plot of feature type 
cnts.moa.type <- moa %>%
  dplyr::group_by(feature_type) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::arrange(desc(n))
cnts.moa.type$feature_type <- factor(cnts.moa.type$feature_type,levels=cnts.moa.type$feature_type)

ggplot(cnts.moa.type,aes(x=feature_type,y=n))+
  #geom_bar(stat_count="identity")+
  geom_point(alpha=.8,color='darkblue')+
  theme_bw()+
  ylab("nuber of MOA entries")+
  ggtitle(paste0("MOAlmanac aberrations by type"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### load MSK-IMPACT
inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/impact_2017_annotated_per_variant.tsv"
msk <- read.csv(inFile,sep="\t")
inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/msk_impact_data_clinical_sample2.txt"
patient.msk <- read.csv(inFile,sep="\t")
msk$Chromosome.str <- paste0("chr",msk$Chromosome)
msk$MAF <- 100*(msk$t_alt_count / (msk$t_ref_count+msk$t_alt_count))

# join clinical info onto variant info
all(msk$Tumor_Sample_Barcode %in% patient.msk$SAMPLE_ID)

msk <- msk %>% 
  dplyr::left_join(patient.msk,by=c("Tumor_Sample_Barcode"="SAMPLE_ID"))

# COMMAND ----------
### add matching columns
moa$source <- "MOAlmanac"
moa$FDAApproved <- NA
moa$Indicated <- moa$predictive_implication
moa$EvidenceText <- moa$description
moa$ReferenceOrTrialID <- moa$citation
moa$Disease <- moa$disease
moa$TrialID <- NA
moa$Phase <- NA
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
clinical$ReferenceOrTrialID <- clinical$gene_civic_url
clinical$Disease <- clinical$disease
clinical$TrialID <- NA
clinical$Phase <- NA
clinical$Drugs <- clinical$drugs
clinical$AberrationsLabelsMod <- NA
clinical$AAChange <- clinical$variant
clinical$pos <- clinical$start
clinical$ref <- clinical$reference_bases
clinical$alt <- clinical$variant_bases


### combine MOA and CIViC data

### create cancer type maps CIVIC-->MSK and MOA-->MSK
civic.cType.cnt <- clinical %>%
  dplyr::group_by(Disease) %>%
  dplyr::summarise(number.of.subjects=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.subjects))

moa.cType.cnt <- moa %>%
  dplyr::group_by(disease) %>%
  dplyr::summarise(number.of.subjects=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.subjects))

### which are the most common cancer types in msk, moa, and civic
msk.cType.cnt <- msk %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::summarise(number.of.subjects=dplyr::n_distinct(PATIENT_ID)) %>%
  dplyr::arrange(desc(number.of.subjects))
outF <-  paste0(figDir,"/msk_cancer_types.txt")
write.table(msk.cType.cnt,outF,row.names=F,quote=F,sep="\t")

## create cancer type maps
iCMmatch <- civic.cType.cnt$Disease %in% msk.cType.cnt$CANCER_TYPE
CivicMSKCancerTypeMap <- data.frame(CivicCancerType=civic.cType.cnt$Disease)
rownames(CivicMSKCancerTypeMap) <- CivicMSKCancerTypeMap$CivicCancerType
CivicMSKCancerTypeMap[iCMmatch,"MskCancerType"] <- as.character(civic.cType.cnt$Disease[iCMmatch])
mapF <-  paste0(figDir,"/civic_cancer_type_map.txt")
write.table(CivicMSKCancerTypeMap,mapF,row.names=F,quote=F,sep="\t")
CivicMSKCancerTypeMap <- read.csv(mapF,sep="\t")
clinical <- clinical %>%
  dplyr::left_join(CivicMSKCancerTypeMap,by=c("disease"="CivicCancerType"))
table(clinical$MskCancerType,exclude=NULL)
# what fraction of entries are still missing MSK cancer type assignment
print(sum(is.na(clinical$MskCancerType))/dim(clinical)[[1]])

iCMmatch <- moa.cType.cnt$disease %in% msk.cType.cnt$CANCER_TYPE
MoaMSKCancerTypeMap <- data.frame(MoaCancerType=moa.cType.cnt$disease)
rownames(MoaMSKCancerTypeMap) <- MoaMSKCancerTypeMap$CivicCancerType
MoaMSKCancerTypeMap[iCMmatch,"MskCancerType"] <- as.character(moa.cType.cnt$disease[iCMmatch])
mapF <-  paste0(figDir,"/moa_cancer_type_map.txt")
write.table(MoaMSKCancerTypeMap,mapF,row.names=F,quote=F,sep="\t")
MoaMSKCancerTypeMap <- read.csv(mapF,sep="\t")
MoaMSKCancerTypeMap$MoaCancerType <- as.character(MoaMSKCancerTypeMap$MoaCancerType)
moa <- moa %>%
  dplyr::left_join(MoaMSKCancerTypeMap,by=c("disease"="MoaCancerType"))
table(moa$MskCancerType,exclude=NULL)
# what fraction of entries are still missing MSK cancer type assignment
print(sum(is.na(moa$MskCancerType))/dim(moa)[[1]])

# now combine MOA and CIVIC
matchCols <- c("source","gene","AAChange","Drugs","FDAApproved","ReferenceOrTrialID","EvidenceText","Phase","Indicated","Disease","chr","pos","ref","alt","MskCancerType")
dbRules <- rbind(moa[,matchCols],clinical[,matchCols])
outF <-  paste0(figDir,"/civic_MOA_clinically_actionable_list.txt")
write.table(dbRules,outF,row.names=F,quote=F,sep="\t")

### summarize
dbGenome <- dbRules %>%
  dplyr::group_by(chr,pos,ref,alt) %>%
  dplyr::summarize(n.db.entries=dplyr::n(),
                   source=paste0(unique(source),collapse=";")) %>% 
  dplyr::arrange(desc(n.db.entries))
dbGenome <- dbGenome[!dbGenome$ref=="" & !dbGenome$alt=="",]
outF <-  paste0(figDir,"/civic_MOA_clinically_actionable_genomic_pos.txt")
write.table(dbGenome,outF,row.names=F,quote=F,sep="\t")

dbAlteration <- dbRules %>%
  dplyr::group_by(gene,AAChange) %>%
  dplyr::summarize(n.db.entries=dplyr::n(),
                   source=paste0(unique(source),collapse=";"),
                   diseases=paste0(unique(Disease),collapse=";"),
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
# output non-match events
outF <-  paste0(figDir,"/civic_MOA_non_AA_change_alterations.txt")
write.table(dbOtherAlt,outF,row.names=F,quote=F,sep="\t")


#table(dbProtein[,"AAChange"]) # except 
#table(dbOtherAlt[,"AAChange"])
table(dbProtein[,"source"]) # except 
table(dbOtherAlt[,"source"])

# COMMAND ----------
### join with MOA - by protein change
dbAlteration$MOAset <- "Y"
msk.protein <- msk %>%
  dplyr::left_join(dbAlteration,by=c("Hugo_Symbol"="gene",
                            "protein_change"="AAChange"))
print(table(msk.protein$MOAset,msk.protein$source,exclude = NULL))
print(table(msk.protein$feature_type,exclude = NULL))

msk.protein.cnt <- msk %>%
  dplyr::inner_join(dbAlteration,by=c("Hugo_Symbol"="gene",
                                     "protein_change"="AAChange")) %>%
  dplyr::group_by(Hugo_Symbol,protein_change) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID),
                   source=paste0(unique(source),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(msk.protein.cnt))
print(table(msk.protein.cnt$source,exclude = NULL))

### impose cancer type matches (protein)
msk.protein.cancer <- msk %>%
  dplyr::inner_join(dbRules,by=c("Hugo_Symbol"="gene",
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

# COMMAND ----------
### join with MOA - by genome position
moa$MOAset <- "Y"
msk.full <- msk %>%
  dplyr::left_join(moa,by=c("Chromosome"="chromosome",
                            "Start_Position"="start_position",
                            "Tumor_Seq_Allele2"="alternate_allele"))
print(table(msk.full$MOAset,exclude = NULL))
print(table(msk.full$feature_type,exclude = NULL))
# number of unique variants
print(length(unique(msk.full[msk.full$MOAset=="Y","HGVSc"])))
msk.moa <- msk.full[!is.na(msk.full$MOAset),]

### impose cancer type matches (genome)
msk.genome.cancer <- msk %>%
  dplyr::inner_join(dbRules,by=c("Chromosome.str"="chr",
                            "Start_Position"="pos",
                            "Tumor_Seq_Allele2"="alt",
                            "CANCER_TYPE"="MskCancerType"))
print(table(msk.genome.cancer$MOAset,exclude = NULL))
msk.genome.cancer.cnt <- msk.genome.cancer %>%
  dplyr::group_by(Chromosome,Start_Position,Tumor_Seq_Allele2,CANCER_TYPE) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(msk.genome.cancer.cnt))
print(table(msk.genome.cancer.cnt$source,exclude = NULL))
print(sum(msk.genome.cancer.cnt$n.patients.mutated))


# number of unique variants that are in a db
msk.genome.cnt <- msk %>%
  dplyr::inner_join(dbGenome,by=c("Chromosome.str"="chr",
                                  "Start_Position"="pos",
                                  "Tumor_Seq_Allele2"="alt")) %>%
  dplyr::group_by(Chromosome,Start_Position,Tumor_Seq_Allele2) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID),
                   source=paste0(unique(source),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(msk.genome.cnt))
print(table(msk.genome.cnt$source,exclude = NULL))





## which genes are seen as mutated? 
gene.cnts.moa <- msk.full %>%
  group_by(gene) %>%
  summarize(n=n()) %>%
  arrange(desc(n)) %>%
  data.frame()
gene.cnts.moa <- gene.cnts.moa[!is.na(gene.cnts.moa$gene),]
table(as.character(gene.cnts.moa$gene))

gene.cnts.moa$gene <- factor(gene.cnts.moa$gene,levels=gene.cnts.moa$gene)
ggplot(gene.cnts.moa,aes(x=gene,y=n))+
  #geom_bar(stat_count="identity")+
  geom_point(alpha=.8,color='brown')+
  theme_bw()+
  ylab("nuber of patients with variant in gene")+
  ggtitle(paste0("MOAlmanac variants appearing in the\n MSK-IMPACT dataset"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

msk.full$full_protein_change <- paste0(msk.full$gene,".",msk.full$protein_change.y)
cnts.moa <- msk.full %>%
  group_by(full_protein_change) %>%
  summarise(n=n()) %>% 
  arrange(desc(n))
cnts.moa <- cnts.moa[!cnts.moa$full_protein_change=="NA.NA",]

# also group by cancer type
msk.full <- msk.full %>% 
  group_by(CANCER_TYPE) %>%
  mutate(n_cType_samples=n_distinct(Tumor_Sample_Barcode)) %>%
  data.frame()


# create action_type
msk.full$action_type <- as.character(msk.full$therapy_type)
msk.full[msk.full$action_type == "" & !is.na(msk.full$action_type),"action_type"] <- "Prognostic" # if no therapy type then label as prognostic
msk.full$action_type <- factor(msk.full$action_type)


cnts.moa.cType <- msk.full %>%
  group_by(full_protein_change,CANCER_TYPE,Hugo_Symbol,action_type) %>%
  summarise(n=n_distinct(Tumor_Sample_Barcode),
            n_cType_samples=n_cType_samples,
            perc_cType_mutated=100*(n/n_cType_samples)) %>% 
  unique() %>%
  arrange(desc(perc_cType_mutated))
cnts.moa.cType <- cnts.moa.cType[!cnts.moa.cType$full_protein_change=="NA.NA",]

iExclude1 <- cnts.moa.cType$n_cType_samples < 40
iExclude2 <- grepl("BRAF",cnts.moa.cType$full_protein_change)
iExclude3 <- grepl("EGFR",cnts.moa.cType$full_protein_change)
cnts.moa.cType.subset <- cnts.moa.cType[!iExclude1 & !iExclude2 & !iExclude3,]
outF <-  paste0("/Users/larsonhogstrom/Documents/variant_annotation/MSK-IMPACT_MOAlmanac_cType.txt")
write.table(cnts.moa.cType.subset,outF,row.names=T,quote=F,sep=",")

### make plot for each gene where each cancer type gets its own dot
outF <-  paste0(figDir,"/msk_percent_mutated_with_actionable.png")
ggplot(cnts.moa.cType,aes(x=CANCER_TYPE,y=perc_cType_mutated,color=action_type))+
  geom_point(alpha=.6)+
  theme_bw()+
  #facet_grid(~Hugo_Symbol,scale="free_x", space="free")+
  facet_wrap(~Hugo_Symbol,scale="free_x")+
  ylab("Percent of patients with variant")+
  ggtitle(paste0("Percent of patients with an actionable mutation\n in MSK-IMPACT dataset"))+
  scale_color_brewer(palette="Set2")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#ggsave(outF,height = 4,width = 25)
ggsave(outF,height = 20,width = 20)

### make gene-specifc plots for KRAS and BRAF
cnts.moa.subset <- cnts.moa.cType[cnts.moa.cType$Hugo_Symbol=="KRAS" | cnts.moa.cType$Hugo_Symbol=="BRAF",]
outF <-  paste0(figDir,"/msk_percent_KRAS_or_BRAF_mutated_with_actionable.png")
ggplot(cnts.moa.subset,aes(x=CANCER_TYPE,y=perc_cType_mutated,color=action_type))+
  geom_point(alpha=.6)+
  theme_bw()+
  facet_grid(Hugo_Symbol ~ .)+
  ylab("Percent of patients with variant")+
  ggtitle(paste0("Percent of patients with an actionable mutation\n in MSK-IMPACT dataset"))+
  scale_color_brewer(palette="Set2",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 7,width = 6)

### make gene-specifc plots all but top N
exclude.gene.group <- c("BRAF","EGFR","PIK3CA","KRAS","NRAS")
cnts.moa.subset <- cnts.moa.cType[!cnts.moa.cType$Hugo_Symbol %in% exclude.gene.group,]
outF <-  paste0(figDir,"/msk_percent_mutated_with_actionable_bottom_genes.png")
ggplot(cnts.moa.subset,aes(x=CANCER_TYPE,y=perc_cType_mutated,color=action_type))+
  geom_point(alpha=.6)+
  theme_bw()+
  facet_grid(Hugo_Symbol ~ .)+
  ylab("Percent of patients with variant")+
  ggtitle(paste0("Percent of patients with an actionable mutation\n in MSK-IMPACT dataset"))+
  scale_color_brewer(palette="Set2",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 11,width = 6)

### cancer types with the largest fraction of actionable mutations
cType.cnts <- cnts.moa.cType %>%
  group_by(CANCER_TYPE) %>%
  mutate(sum.perc.mutated=sum(perc_cType_mutated)) %>%
  arrange(desc(sum.perc.mutated))

subsetList <- as.character(unique(cType.cnts$CANCER_TYPE))[1:8]
subsetList <- subsetList[!subsetList=="Leukemia"] # exclude n=1 leukemia entry
cType.cnts.subset <- cType.cnts[cType.cnts$CANCER_TYPE %in% subsetList & cType.cnts$action_type %in% c("Prognostic","Targeted therapy"),]
cType.cnts.subset$CANCER_TYPE <- factor(cType.cnts.subset$CANCER_TYPE,levels=subsetList)
cType.cnts.subset$action_type <- factor(cType.cnts.subset$action_type,levels=c("Prognostic","Targeted therapy"))
cType.cnts.subset$percent_mutated <- cType.cnts.subset$perc_cType_mutated

outF <-  paste0(figDir,"/msk_percent_mutated_by_cancer_type.png")
ggplot(cType.cnts.subset,aes(x=Hugo_Symbol,y=CANCER_TYPE,fill=percent_mutated))+
  geom_tile()+
  theme_bw()+
  facet_grid(action_type ~ .)+
  xlab(" ")+
  ylab(" ")+
  scale_color_brewer(palette="Set2",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 6,width = 8)


### how often do subjects have one or more actionable mutations - rank by cancer type
actionable.counts <- msk.full %>%
  group_by(CANCER_TYPE,Tumor_Sample_Barcode) %>% 
  filter(!full_protein_change == "NA.NA") %>%
  summarise(number.of.actionable.mutations=n_distinct(full_protein_change),
            changes = paste0(full_protein_change, collapse = ";")) %>%
  arrange(desc(number.of.actionable.mutations))
print(table(actionable.counts$number.of.actionable.mutations))
actionable.hist <- as.data.frame(table(actionable.counts$number.of.actionable.mutations))
colnames(actionable.hist) <- c("number.of.actionable.mutations","number.of.subjects")
actionable.hist$number.of.actionable.mutations <- as.numeric(actionable.hist$number.of.actionable.mutations)
num.subjects <- length(unique(msk.full$Tumor_Sample_Barcode))
actionable.hist[5,] <- c(0,num.subjects-dim(actionable.counts)[1])


outF <-  paste0(figDir,"/msk_percent_subjects_with_actionable.png")
ggplot(actionable.hist,aes(x=number.of.actionable.mutations,y=number.of.subjects))+
  geom_col(alpha=.6)+
  theme_bw()+
  ylab("number of subjects")+
  xlab("Number of actionable mutations")+
  ggtitle(paste0("Distribution of actionable mutations \n found in MSK-IMPACT patients"))+
  theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 4,width = 4)


  
# COMMAND ----------
### parse syn/non-syn ratios 
inFile <- paste0(figDir,"/Li_2018_TableS2.txt")
background <- read.csv(inFile,sep="\t")



syn.ratios <- background[,c("Gene","SubType","UniqueMutations","BinnedRole")] %>%
  tidyr::pivot_wider(names_from = SubType, values_from = UniqueMutations)
syn.ratios[is.na(syn.ratios$Silent),"Silent"] <- 0
syn.ratios[is.na(syn.ratios$Nonsense),"Nonsense"] <- 0
syn.ratios[is.na(syn.ratios$Missense),"Missense"] <- 0
syn.ratios$Nonsynonymous <- syn.ratios$Nonsense + syn.ratios$Missense
syn.ratios$SynToNonSynRatio <- (syn.ratios$Silent+1)/(syn.ratios$Nonsynonymous+1)

syn.ratios.subset <- syn.ratios[syn.ratios$Gene %in% msk.full$Hugo_Symbol,]

### get msk mutation counts
msk.gene.cnts <- msk.full %>%
  group_by(Hugo_Symbol) %>%
  summarise(n=n_distinct(Tumor_Sample_Barcode)) 

msk.gene.cnts <- msk.gene.cnts %>%
  left_join(syn.ratios.subset,by=c("Hugo_Symbol"="Gene"))

outF <-  paste0(figDir,"/msk_mut_counts_syn_non_sym_ratio.png")
ggplot(msk.gene.cnts,aes(x=n,y=SynToNonSynRatio,color=BinnedRole))+
  geom_point(alpha=.6)+
  theme_bw()+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  xlab("number of variants in MSK-IMPACT")+
  #ggtitle(paste0("Percent of patients with an actionable mutation\n in MSK-IMPACT dataset"))+
  scale_color_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 5,width = 5)

# COMMAND ----------
### join with CIViC
clinical$civicSet <- "Y"
msk.full <- msk.full %>%
  dplyr::left_join(clinical,by=c("Chromosome"="chromosome",
                                 "Start_Position"="start",
                                 "Tumor_Seq_Allele2"="variant_bases"))
print(table(msk.full$civicSet,exclude = NULL))
print(table(msk.full$feature_type,exclude = NULL))
# number of unique variants
print(length(unique(msk.full[msk.full$civicSet=="Y","HGVSc"])))

### gene counts
gene.cnts.civic <- msk.full %>%
  group_by(gene.y) %>%
  summarize(n=n()) %>%
  arrange(desc(n)) %>%
  data.frame()
gene.cnts.civic <- gene.cnts.civic[!is.na(gene.cnts.civic$gene),]
table(as.character(gene.cnts.civic$gene))

gene.cnts.civic$gene <- factor(gene.cnts.civic$gene,levels=gene.cnts.civic$gene)
ggplot(gene.cnts.civic,aes(x=gene,y=n))+
  #geom_bar(stat_count="identity")+
  geom_point(alpha=.8,color='darkred')+
  theme_bw()+
  ylab("nuber of patients with variant in gene")+
  ggtitle(paste0("CIViC variants appearing in the\n MSK-IMPACT dataset"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )

# COMMAND ----------
### PARP analysis
table(dbRules$gene)
# number of BRCA1 entries 42 and BRCA2 is 40
iBraca <- dbRules$gene=="BRCA1" | dbRules$gene=="BRCA2"
brca.entries <- dbRules[iBraca,]
print(table(brca.entries$AAChange))

iBraca.moa <- moa$gene=="BRCA1" | moa$gene=="BRCA2"
brca.moa <- moa[iBraca.moa,]
table(as.character(brca.moa$Drugs))
table(as.character(brca.moa$disease))

iBraca.civic <- clinical$gene=="BRCA1" | clinical$gene=="BRCA2"
brca.civic <- clinical[iBraca.civic,]
table(as.character(brca.civic$AAChange))

# COMMAND ----------
### MAF results
maf.df.msk <- data.frame(MAF=msk$MAF)
maf.df.msk$source <- "MSK all variants"

# maf of clinically actionable based on protein match
maf.df.protein<- data.frame(MAF=msk.protein.cancer$MAF)
maf.df.protein$source <- "clinically actionable"

maf.df <- rbind(maf.df.msk,maf.df.protein)

outF <-  paste0(figDir,"/MAF_distribution_clinically_actionable.png")
ggplot(maf.df,aes(x=MAF,fill=source))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  facet_grid(source~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n in MSK-IMPACT"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 5,width = 5)

# COMMAND ----------
### combination of mutaitons

# what are the most common pairs of mutations? 
msk.protein.cnt.df <- data.frame(msk.protein.cnt)
msk.protein.cnt.df$var.id <- paste0(msk.protein.cnt.df$Hugo_Symbol,"-",msk.protein.cnt.df$protein_change)
msk.protein.cancer$var.id <- paste0(msk.protein.cancer$Hugo_Symbol,"-",msk.protein.cancer$protein_change)
mut.pairwise <- data.frame()
for (i in rownames(msk.protein.cnt.df)) {
  for (j in rownames(msk.protein.cnt.df)) {
    iID <- msk.protein.cnt.df[i,"var.id"]
    jID <- msk.protein.cnt.df[j,"var.id"]
    #iGene <- msk.protein.cnt.df[i,"Hugo_Symbol"]
    #yGene <- msk.protein.cnt.df[j,"Hugo_Symbol"]
    #iAA <- msk.protein.cnt.df[i,"protein_change"]
    #jAA <- msk.protein.cnt.df[j,"protein_change"]
    
    #what fraction of patients contain both variants? 
    iMatch1 <- msk.protein.cancer$var.id == iID
    iMatch2 <- msk.protein.cancer$var.id == jID
    df.match <- msk.protein.cancer[iMatch1 | iMatch2,]
    subj.cnts <- df.match %>%
      dplyr::group_by(PATIENT_ID) %>%
      dplyr::summarize(n.muts=dplyr::n_distinct(var.id))
    # how many subjects had both mutations
    n.pairs <- sum(subj.cnts$n.muts==2)
    mut.pairwise[i,j] <- n.pairs
  }
}
rownames(mut.pairwise) <- msk.protein.cnt.df$var.id
colnames(mut.pairwise) <- msk.protein.cnt.df$var.id
rownames(msk.protein.cnt.df) <- msk.protein.cnt.df$var.id

# plot raw output
heatmapFile=paste0(figDir,"/mutation_pair_clustering.pdf")
gcdfOut <- pheatmap(as.matrix(t(mut.pairwise)),annotation = msk.protein.cnt.df[,c("Hugo_Symbol","protein_change")],#annotation_colors =ann_colors,
                    show_rownames = F,#labels_row = rownames(raw.prot),
                    show_colnames = F,#labels_col = meta[gnums,"name"],
                    color = rev(brewer.pal(8,"GnBu")),fontsize_col = 5,
                    cluster_cols=T,cluster_rows=T,
                    filename = heatmapFile,width = 18, height = 10)

# exclude zero rows/cols
iColsum <- colSums(mut.pairwise) > 5
iRowsum <- rowSums(mut.pairwise) > 5
mut.pairwise.nonzero <- mut.pairwise[iRowsum,iColsum]

heatmapFile=paste0(figDir,"/mutation_pair_clustering_nonzero.pdf")
gcdfOut <- pheatmap(as.matrix(t(mut.pairwise.nonzero)),
                    annotation = msk.protein.cnt.df[rownames(mut.pairwise.nonzero),c("Hugo_Symbol","source")],#annotation_colors =ann_colors,
                    show_rownames = T,labels_row = rownames(mut.pairwise.nonzero),
                    show_colnames = T,labels_col = colnames(mut.pairwise.nonzero),
                    color = rev(brewer.pal(8,"GnBu")),#fontsize_col = 5,fontsize_row = 5,
                    cluster_cols=T,cluster_rows=T,
                    filename = heatmapFile,width = 12, height = 10)


### repeat combined analysis with short list
inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/ACARE_MOAlmanac_CIVIC_COSMIC_AA_match_in_V1_ROI_summary_with_wildcard.txt"
shortlist <- read.csv(inFile,sep="\t")
print(dim(shortlist))


msk.shortlist.cnt <- msk %>%
  dplyr::inner_join(shortlist,by=c("Hugo_Symbol"="gene",
                                      "HGVSp_Short"="Mutation_AA")) %>%
  dplyr::group_by(Hugo_Symbol,protein_change) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID)) %>%#,
                   #source=paste0(unique(source),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(msk.shortlist.cnt))
print(table(msk.shortlist.cnt$source,exclude = NULL))
msk.shortlist.cancer <- msk %>%
  dplyr::inner_join(shortlist,by=c("Hugo_Symbol"="gene",
                                 "HGVSp_Short"="Mutation_AA")) 

msk.protein.cnt.df <- data.frame(msk.protein.cnt)
msk.protein.cnt.df$var.id <- paste0(msk.protein.cnt.df$Hugo_Symbol,"-",msk.protein.cnt.df$protein_change)
msk.shortlist.cancer$var.id <- paste0(msk.shortlist.cancer$Hugo_Symbol,"-",msk.shortlist.cancer$protein_change)
mut.pairwise <- data.frame()
for (i in rownames(msk.protein.cnt.df)) {
  for (j in rownames(msk.protein.cnt.df)) {
    iID <- msk.protein.cnt.df[i,"var.id"]
    jID <- msk.protein.cnt.df[j,"var.id"]
    #iGene <- msk.protein.cnt.df[i,"Hugo_Symbol"]
    #yGene <- msk.protein.cnt.df[j,"Hugo_Symbol"]
    #iAA <- msk.protein.cnt.df[i,"protein_change"]
    #jAA <- msk.protein.cnt.df[j,"protein_change"]
    
    #what fraction of patients contain both variants? 
    iMatch1 <- msk.shortlist.cancer$var.id == iID
    iMatch2 <- msk.shortlist.cancer$var.id == jID
    df.match <- msk.shortlist.cancer[iMatch1 | iMatch2,]
    subj.cnts <- df.match %>%
      dplyr::group_by(PATIENT_ID) %>%
      dplyr::summarize(n.muts=dplyr::n_distinct(var.id))
    # how many subjects had both mutations
    n.pairs <- sum(subj.cnts$n.muts==2)
    mut.pairwise[i,j] <- n.pairs
  }
}
rownames(mut.pairwise) <- msk.protein.cnt.df$var.id
colnames(mut.pairwise) <- msk.protein.cnt.df$var.id
rownames(msk.protein.cnt.df) <- msk.protein.cnt.df$var.id

# plot raw output
heatmapFile=paste0(figDir,"/mutation_shortlist_pair_clustering.pdf")
gcdfOut <- pheatmap(as.matrix(t(mut.pairwise)),annotation = msk.protein.cnt.df[,c("Hugo_Symbol","protein_change")],#annotation_colors =ann_colors,
                    show_rownames = F,#labels_row = rownames(raw.prot),
                    show_colnames = F,#labels_col = meta[gnums,"name"],
                    color = rev(brewer.pal(8,"GnBu")),fontsize_col = 5,
                    cluster_cols=T,cluster_rows=T,
                    filename = heatmapFile,width = 18, height = 10)

# exclude zero rows/cols
iColsum <- colSums(mut.pairwise) > 5
iRowsum <- rowSums(mut.pairwise) > 5
mut.pairwise.nonzero <- mut.pairwise[iRowsum,iColsum]

heatmapFile=paste0(figDir,"/mutation_shortlist_pair_clustering_nonzero.pdf")
gcdfOut <- pheatmap(as.matrix(t(mut.pairwise.nonzero)),
                    #annotation = msk.protein.cnt.df[rownames(mut.pairwise.nonzero),c("Hugo_Symbol","source")],#annotation_colors =ann_colors,
                    show_rownames = T,labels_row = rownames(mut.pairwise.nonzero),
                    show_colnames = T,labels_col = colnames(mut.pairwise.nonzero),
                    color = rev(brewer.pal(8,"GnBu")),#fontsize_col = 5,fontsize_row = 5,
                    cluster_cols=T,cluster_rows=T,
                    filename = heatmapFile,width = 12, height = 10)

xy <- t(combn(colnames(mut.pairwise.nonzero), 2))
mut.list <- data.frame(xy, cnt=mut.pairwise.nonzero[xy]) %>%
  dplyr::arrange(desc(cnt))
outF <-  paste0(figDir,"/MSK_top_shortlist_mutation_combinations.txt")
write.table(mut.list,outF,row.names=F,quote=F,sep="\t")

### write out CNV info
#iAmp <- clinical$variant == "Amplification" | clinical$variant == "AMPLIFICATION"
#iDel <- clinical$variant == "Deletion" | clinical$variant == "DELETION"
#cn.clinical <- clinical[iAmp | iDel, ]
#outFile <- paste0(figDir,"/CIViC_copy_number.txt")
#write.table(cn.clinical,outFile,sep="\t",row.names = F)
