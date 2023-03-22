#library(vcfR)
library(pryr)
library(stringr)
library(ggplot2)
library(DBI)
library(RSQLite)

baseDir <- "/data/larsonh/hartwig"
#baseDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Hartwig/data" 
figDir <- paste0(baseDir,"/output_v4")

inFile <- paste0(baseDir,"/samples.txt")
sampleList <- read.csv(inFile,sep=",", header = F)

inFile <- paste0(baseDir,"/metadata.tsv")
hMeta <- read.csv(inFile,sep="\t")

inFile <- paste0(baseDir,"/pre_biopsy_drugs.tsv")
hTreatment <- read.csv(inFile,sep="\t")

########################################
### create connection to results db ###
########################################

bDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/resultsDb"
mydb <- DBI::dbConnect(RSQLite::SQLite(), paste0(bDir,"/actionable-biomarker-db.sqlite"))

dbWriteTable(mydb, "mtcars", mtcars)

dbDisconnect(mydb)

########################################
### Load PCGR output for each sample ###
########################################

# load list of clinically actionable entries
#inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/mutation_incidence_20220913/civic_MOA_clinically_actionable_list.txt"
inFile <- "/data/larsonh/analysis/mutation_incidence_20220913/civic_MOA_clinically_actionable_list.txt"
dbRules <- read.csv(inFile,sep="\t")

### create fields for AA and genomic pos matching
dbRules$AAMatchStr <- paste0(dbRules$gene,"-",dbRules$AAChange)
dbRules$GenPosMatchStr <- paste0(dbRules$chr,"-",
                                 dbRules$pos,"-",
                                 dbRules$ref,"-",
                                 dbRules$alt)

cVariantTable <- data.frame()
ncVariantTable <- data.frame()
varsPerSubject <- data.frame()
for (sample in as.character(sampleList$V1)) {
  print(sample)
  inFile <- paste0(baseDir,"/",sample,"/purple/pcgr_out/",sample,".pcgr_acmg.grch37.pass.tsv.gz")
  
  if (file.info(inFile)$size > 15000000) { # if file is too big then skip
    print(paste0(sample, " PCGR output exceeds size limit: 10.5MB"))
  } else {
    df.pcgr <- read.csv(inFile,sep="\t",skip=1) %>%
      dplyr::mutate(AAChange=stringr::str_replace(HGVSp_short,"p.", "")) %>%
      data.frame()
    df.pcgr$sample <- sample
    #colSubset <- c("TVAF")
    print(dim(df.pcgr))
    # track variant counts per subject
    #varsPerSubject[sample,"subject"] <- sample
    #varsPerSubject[sample,"numberOfVariants"] <- dim(df.pcgr)[[1]]
    # variant summary info
    varSum <- df.pcgr %>%
      dplyr::summarise(numberOfVariants=dplyr::n(),
                       TSGEvents=sum(TUMOR_SUPPRESSOR=="True"),
                       ONCEvents=sum(TUMOR_SUPPRESSOR=="True"),
                       CodingEvents=sum(CODING_STATUS=="coding"),
                       HotspotEvents=sum(!MUTATION_HOTSPOT=="."),
                       LoFEvents=sum(!LoF=="."),
                       CanonicalEvents=sum(!CANONICAL=="."),
                       TCGADriverEvents=sum(!TCGA_DRIVER=="."))
    varSum$sample <- sample
    
    ### Look for clinically actionable entries ###
    # make match strings -- AA change
    df.pcgr$AAMatchStr <- paste0(as.character(df.pcgr$SYMBOL),"-",df.pcgr$AAChange)
    df.pcgr[df.pcgr$HGVSp_short==".","AAMatchStr"]<- "N/A"
    # make match strings -- genomic coordinate
    df.pcgr$GenPosMatchStr <- paste0("chr",as.character(df.pcgr$CHROM),"-"
                                     ,as.character(df.pcgr$POS),"-"
                                     ,as.character(df.pcgr$REF),"-"
                                     ,as.character(df.pcgr$ALT))
    
    # check for clinically actionable genome or AA matches
    iAAMatch <- df.pcgr$AAMatchStr %in% dbRules$AAMatchStr
    df.pcgr$AAClinEntryMatch <- F
    df.pcgr[iAAMatch,"AAClinEntryMatch"] <- T
    
    iGenomeMatch <- df.pcgr$GenPosMatchStr %in% dbRules$GenPosMatchStr
    df.pcgr$GenomeClinEntryMatch <- F
    df.pcgr[iGenomeMatch,"GenomeClinEntryMatch"] <- T
    
    iClinMatch <- iAAMatch | iGenomeMatch
    
    print(paste0("number of clinically actionable matches: ", as.character(sum(iClinMatch))))
    outCols <- colnames(df.pcgr)[!colnames(df.pcgr) %in% c("TNC")] # exclude TNC column present in only some subjects
    cVariantTable <- rbind(cVariantTable,df.pcgr[iClinMatch,outCols])
    varSum$numberClinicallyActionableVariants<- sum(iClinMatch)
    varsPerSubject <- rbind(varsPerSubject,varSum)
    
    ### track the MAF values of N non-actionable variants as background
    outColsRand <- c("CHROM",
                     "POS",
                     "ID",
                     "REF",
                     "ALT",
                     "TVAF",
                     "TDP",
                     "CVAF",
                     "CDP")
    nVarToGrab <- 15 # select N per sample coding + non/coding
    iCoding <- df.pcgr$CODING_STATUS=="coding"
    # non-coding
    ncDfNonCoding <- df.pcgr[!iClinMatch & !iCoding,outColsRand]
    nVarToGrabNC <- min(nVarToGrab,dim(ncDfNonCoding)[[1]])
    iRandVarsNC <- sample(dim(ncDfNonCoding)[[1]],nVarToGrabNC,replace=F)
    ncOut <- ncDfNonCoding[iRandVarsNC,]
    ncOut$coding_status <- "non-coding"
    # coding
    ncDfCoding <- df.pcgr[!iClinMatch & iCoding,outColsRand]
    nVarToGrabC <- min(nVarToGrab,dim(ncDfCoding)[[1]])
    cOut <- data.frame()
    if (nVarToGrabC > 0) {
      iRandVarsC <- sample(dim(ncDfCoding)[[1]],nVarToGrabC,replace=F)
      cOut <- ncDfCoding[iRandVarsC,]
      cOut$coding_status <- "coding"
    }
    ncVariantTable <- rbind(ncVariantTable,ncOut,cOut)
  }
  
}
outFile <- paste0(figDir,"/hartwig_clinically_actionable_pcgr_entries.txt")
write.table(cVariantTable,outFile,row.names=F,quote=F,sep="\t")
#cVariantTable <- read.csv(outFile,sep="\t") # get local cashed file

outFile <- paste0(figDir,"/hartwig_random_non_clinically_actionable_variants_per_patient.txt")
write.table(ncVariantTable,outFile,row.names=F,quote=F,sep="\t")
#ncVariantTable <- read.csv(outFile,sep="\t") # get local cashed file

outFile <- paste0(figDir,"/hartwig_variant_counts_per_subject.txt")
write.table(varsPerSubject,outFile,row.names=F,quote=F,sep="\t")
#varsPerSubject <- read.csv(outFile,sep="\t") # get local cashed file

#########################
### metadata analysis ###
#########################

iInSamples <- hMeta$sampleId %in% sampleList$V1
table(iInSamples)

# samples per patient
subject.sample.cnts <- hMeta %>%
  dplyr::group_by(hmfPatientId) %>%
  dplyr::summarise(n=dplyr::n())
table(subject.sample.cnts$n)

# primary type counts
primary.tissue.cnts <- hMeta %>%
  dplyr::group_by(primaryTumorType) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::arrange(desc(n))
table(primary.tissue.cnts$n)
outFile <- paste0(figDir,"/hartwig_metadata_top_primaryTymorTypes.txt")
write.table(primary.tissue.cnts,outFile,row.names=F,quote=F,sep="\t")

# primary location counts
primary.location.cnts <- hMeta %>%
  dplyr::group_by(primaryTumorLocation) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::arrange(desc(n))
table(primary.location.cnts$n)
outFile <- paste0(figDir,"/hartwig_metadata_top_primaryTumorLocation.txt")
write.table(primary.location.cnts,outFile,row.names=F,quote=F,sep="\t")

# primary location counts
biopsySite.cnts <- hMeta %>%
  dplyr::group_by(biopsySite) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::arrange(desc(n))
table(biopsySite.cnts$n)
outFile <- paste0(figDir,"/hartwig_metadata_biopsySite.txt")
write.table(biopsySite.cnts,outFile,row.names=F,quote=F,sep="\t")

### tumor purity
#tumorPurity

#######################################################
### join clinically actionable info to variant info ###
#######################################################

# inFile <- paste0(baseDir,"/output_v4/hartwig_clinically_actionable_pcgr_entries.txt")
# cVariantTable <- read.csv(inFile, sep="\t")

# inFile <- paste0(figDir,"/hartwig_variant_counts_per_subject.txt")
# varsPerSubject <- read.csv(inFile, sep="\t")

aa.genome.cnts <- cVariantTable %>%
  dplyr::group_by(AAClinEntryMatch,GenomeClinEntryMatch) %>%
  dplyr::summarise(n=dplyr::n())

aa.rank <- cVariantTable %>%
  dplyr::group_by(AAMatchStr) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::arrange(desc(n))
outFile <- paste0(figDir,"/hartwig_clinically_actionable_AA_match_rank.txt")
write.table(aa.rank,outFile,row.names=F,quote=F,sep="\t")

gene.rank <- cVariantTable %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(n=dplyr::n()) %>%
  dplyr::arrange(desc(n))
outFile <- paste0(figDir,"/hartwig_clinically_actionable_gene_rank.txt")
write.table(gene.rank,outFile,row.names=F,quote=F,sep="\t")


### join on metadata
dim(cVariantTable)
varTbl <- cVariantTable %>%
  dplyr::left_join(hMeta,by=c("sample"="sampleId"))
dim(varTbl)

table(varTbl$primaryTumorLocation)

### join on  clinically actionable
var.match <- varTbl %>%
  dplyr::inner_join(dbRules,by="AAMatchStr")

aa.match.cnts <- var.match %>%
  dplyr::group_by(sample,AAMatchStr) %>%
  dplyr::summarise(nMatches=dplyr::n(),
                   nDisease=dplyr::n_distinct(Disease))

table(var.match$Disease == var.match$primaryTumorLocation,exclude=NULL)

iLungVar <- var.match$primaryTumorLocation == "Lung"
iLungDB <- var.match$Disease %in% c("Lung",
                                    "Lung Adenocarcinoma",
                                    "Lung Cancer",
                                    "Lung Non-small Cell Carcinoma")
lung.match <- var.match[iLungVar & iLungDB,]
head(lung.match[,c("AAMatchStr",colnames(dbRules))])
## select representative prognostic entry
lung.match[3,200:228]
lung.match[3,"sample"]


iCRCVar <- var.match$primaryTumorLocation == "Colorectum"
iCRCDB <- var.match$Disease %in% c("Colorectal Cancer")
crc.match <- var.match[iCRCVar & iCRCDB,]

crc.match[crc.match$AAMatchStr=="PIK3CA-E545K","sample"]

## corr between variant count and # actionable
outF <-  paste0(figDir,"/number_var_number_actionable.png")
ggplot(varsPerSubject,aes(x=numberOfVariants,y=numberClinicallyActionableVariants))+
  geom_jitter(alpha=.2,height=.2)+
  theme_bw()+
  #facet_grid(metric~.,scale="free_y",)+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  xlab("raw variant count")+
  ylab("Number of actionable mutations")+
  ggtitle(paste0("Hartwig SNV and indel events"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 8,width = 5)



# ###############################################
# ### subset to variant and cancer type match ###
# ###############################################
#
#
#
#
# outF <-  paste0(figDir,"/variants_per_subject.png")
# ggplot(varsPerSubject,aes(x=numberOfVariants))+
#   geom_histogram(alpha=.4)+
#   theme_bw()+
#   #facet_grid(source~.,scale="free_y",)+
#   xlab("variants per sample")+
#   ggtitle(paste0("100 Hartwig samples\n variants per sample"))+
#   scale_fill_brewer(palette="Set1",drop=FALSE)+
#   theme(plot.title = element_text(hjust = 0.5))
#   #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   ggsave(outF,height = 5,width = 5)
#
# outF <-  paste0(figDir,"/variants_per_subject_log.png")
#   ggplot(varsPerSubject,aes(x=numberOfVariants))+
#     geom_histogram(alpha=.4)+
#     theme_bw()+
#     #facet_grid(source~.,scale="free_y",)+
#     scale_x_log10(
#       breaks = scales::trans_breaks("log10", function(x) 10^x),
#       labels = scales::trans_format("log10", scales::math_format(10^.x))
#     )+
#     xlab("variants per sample")+
#     ggtitle(paste0("100 Hartwig samples\n variants per sample"))+
#     scale_fill_brewer(palette="Set1",drop=FALSE)+
#     theme(plot.title = element_text(hjust = 0.5))
#   #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   ggsave(outF,height = 5,width = 5)
#
#
#
# #####################
# ### fusion events ###
# #####################
#
# #inFile <- "/data/larsonh/hartwig/ACTN01020005T/linx/ACTN01020005T.linx.fusion.tsv"
# #fusion.df <- read.csv(inFile,sep="\t")
#
ca.fusion.gene.list <- c("RET","ALK","ROS1","FGFR2","FGFR3","NTRK1","PDGFRA","PDGFRB") # list genes with clinically actionable fusions
ca.fusion.gene.concat <- paste(ca.fusion.gene.list, collapse = '|')
ca.fusion.gene.pairs <- c("BCR_ABL1","EML4_ALK","BCR_PDGFRA","FGFR2_TACC3","FGFR3_NSD2","TMPRSS2_ERG","FIP1L1_PDGFRA","COL1A1_PDGFRB","EWSR1_FLI1","ESRP1_RAF1","RUNX1_RUNX1T1","SLC45A3_BRAF")

fusionsPerSubject <- data.frame()
caFusions <- data.frame()
for (sample in as.character(sampleList$V1)) {
  print(sample)
  inFile <- paste0(baseDir,"/",sample,"/linx/",sample,".linx.fusion.tsv")
  fusion.df <- read.csv(inFile,sep="\t")
  if (dim(fusion.df)[[1]] > 0) {
    # per subject summary
    varSum <- fusion.df %>%
      dplyr::summarise(numberOfFusions=dplyr::n(),
                       PhasedInframeEvents=sum(phased=="INFRAME"),
                       ReportedEvents=sum(reported=="false"),
                       NonZeroChainLength=sum(chainLength>0))
    varSum$sample <- sample

    iSingleGeneMatch <- grepl(ca.fusion.gene.concat,fusion.df$name)
    fusion.df$ClinEntryMatch <- F
    fusion.df[iSingleGeneMatch,"SignleGeneClinEntryMatch"] <- T

    iPairGeneMatch <- fusion.df$name %in% ca.fusion.gene.pairs
    fusion.df$ClinEntryMatch <- F
    fusion.df[iPairGeneMatch,"GenePairClinEntryMatch"] <- T
    fusion.df$sample <- sample

    iClinMatch <- iSingleGeneMatch | iPairGeneMatch

    print(paste0("number of clinically actionable fusion matches: ", as.character(sum(iClinMatch))))
    if (sum(iClinMatch) > 0) {
      caFusions <- rbind(caFusions,fusion.df[iClinMatch,])
    }

    varSum$numberClinicallyActionableVariants<- sum(iClinMatch)
    fusionsPerSubject <- rbind(fusionsPerSubject,varSum)
  }
}

outFile <- paste0(figDir,"/hartwig_clinically_actionable_fusions_linx.txt")
write.table(caFusions,outFile,row.names=F,quote=F,sep="\t")

outFile <- paste0(figDir,"/hartwig_fusions_linx_summary.txt")
write.table(fusionsPerSubject,outFile,row.names=F,quote=F,sep="\t")

#######################
### fusion analysis ###
#######################

# inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Hartwig/data/output_v2/hartwig_fusions_linx_summary.txt"
# fusionsPerSubject <- read.csv(inFile,sep="\t")
# 
# inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Hartwig/data/output_v2/hartwig_clinically_actionable_fusions_linx.txt"
# caFusions <- read.csv(inFile,sep="\t")

# select a few fusion features for plotting
fussion.cnt.long <- fusionsPerSubject[,c("numberOfFusions","PhasedInframeEvents","NonZeroChainLength")] %>%
  tidyr::pivot_longer(c("numberOfFusions","PhasedInframeEvents","NonZeroChainLength"),
                      names_to = "metric", values_to = "value")
fussion.cnt.long$metric <- factor(fussion.cnt.long$metric,levels=c("numberOfFusions","PhasedInframeEvents","NonZeroChainLength"))

outF <-  paste0(figDir,"/fussion_metrics_per_subject_log.png")
  ggplot(fussion.cnt.long,aes(x=value),fill=metric)+
    geom_histogram(alpha=.4)+
    theme_bw()+
    facet_grid(metric~.,scale="free_y",)+
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )+
    xlab("events per sample")+
    ggtitle(paste0("Hartwig fussion events"))+
    scale_fill_brewer(palette="Set1",drop=FALSE)+
    theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 8,width = 5)

# look at CA fusions
fusion.full <- caFusions %>%
  dplyr::left_join(hMeta,by=c("sample"="sampleId"))# %>%
  #dplyr::inner_join(dbRules,by="AAMatchStr") # need to add MOA

iMatchPair <- caFusions$name %in% ca.fusion.gene.pairs
table(iMatchPair)

iMatchALKPair <- fusion.full$name %in% c("EML4_ALK")
iMatchLung <- fusion.full$primaryTumorLocation == "Lung"
exCAFusion <- fusion.full[iMatchALKPair & iMatchLung,]
table(iMatchALKPair)
head(caFusions[iMatchALKPair,])

######################################
### Hartwig patient treatment data ###
######################################

print(dim(hTreatment))
print(length(unique(hTreatment$patientIdentifier)))

per.patient.summary <- hTreatment %>%
  dplyr::group_by(patientIdentifier) %>%
  dplyr::summarise(number_of_therapies=dplyr::n_distinct(name),
                   number_of_therapy_class=dplyr::n_distinct(type),
                   number_of_mechanism=dplyr::n_distinct(mechanism))

plt.df <- per.patient.summary %>%
  tidyr::pivot_longer(!patientIdentifier,names_to = "class",values_to = "value")

outF <-  paste0(figDir,"/pre_biopsy_drug_treatment_counts.png")
ggplot(plt.df, aes(x=value)) + 
  geom_histogram(bins=33)+
  facet_grid(class ~ .)+
  xlim(c(0,13))+
  theme_bw()+
  xlab("number of drugs")+
  ylab("number of patients")+
  ggtitle("Hartwig pre-biopsy therapies per subject")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(outF,height = 6,width = 5)
