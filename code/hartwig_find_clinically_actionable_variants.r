library(vcfR)
library(pryr)
library(stringr)
library(ggplot2)

baseDir <- "/data/larsonh/hartwig"
#baseDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Hartwig/data" 
figDir <- paste0(baseDir,"/output_v2")

inFile <- paste0(baseDir,"/samples.txt")
sampleList <- read.csv(inFile,sep=",", header = F)

inFile <- "/data/larsonh/hartwig/metadata.tsv"
hMeta <- read.csv(inFile,sep="\t")

########################################
### Load PCGR output for each sample ###
########################################

# load list of clinically actionable entries
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
outFile <- paste0(figDir,"/hartwig_clinically_actionable_pcgr_entries.txt")
write.table(cVariantTable,outFile,row.names=F,quote=F,sep="\t")

outFile <- paste0(figDir,"/hartwig_random_non_clinically_actionable_variants_per_patient.txt")
write.table(cVariantTable,outFile,row.names=F,quote=F,sep="\t")

outFile <- paste0(figDir,"/hartwig_variant_counts_per_subject.txt")
write.table(varsPerSubject,outFile,row.names=F,quote=F,sep="\t")

#######################################################
### join clinically actionable info to variant info ###
#######################################################



###############################################
### subset to variant and cancer type match ###
###############################################




outF <-  paste0(figDir,"/variants_per_subject.png")
ggplot(varsPerSubject,aes(x=numberOfVariants))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(source~.,scale="free_y",)+
  xlab("variants per sample")+
  ggtitle(paste0("100 Hartwig samples\n variants per sample"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 5,width = 5)
  
outF <-  paste0(figDir,"/variants_per_subject_log.png")
  ggplot(varsPerSubject,aes(x=numberOfVariants))+    
    geom_histogram(alpha=.4)+
    theme_bw()+
    #facet_grid(source~.,scale="free_y",)+
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )+
    xlab("variants per sample")+
    ggtitle(paste0("100 Hartwig samples\n variants per sample"))+
    scale_fill_brewer(palette="Set1",drop=FALSE)+
    theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 5,width = 5)  

#########################
### metadata analysis ###
#########################

subject.sample.cnts <- hMeta %>%
  dplyr::group_by(hmfPatientId) %>%
  dplyr::summarise(n=dplyr::n())
table(subject.sample.cnts$n)


### scratch code from previous matching efforts

# ### aggregate across samples 
# msk.genome.cancer.cnt <- msk.genome.cancer %>%
#   dplyr::group_by(Chromosome,Start_Position,Tumor_Seq_Allele2,CANCER_TYPE) %>%
#   dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID),
#                    source=paste0(unique(source),collapse=";"),
#                    Drugs=paste0(unique(Drugs),collapse=";"),
#                    ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
#   dplyr::arrange(desc(n.patients.mutated))
# # join clinically actionable alterations
# df.protein <- df.pcgr %>%
#   dplyr:inner_join(dbRules,by=c("Hugo_Symbol"="gene",
#                                      "protein_change"="AAChange")) 
# 
# ### join with MOA - by protein change
# dbAlteration$MOAset <- "Y"
# msk.protein <- msk %>%
#   dplyr::left_join(dbAlteration,by=c("Hugo_Symbol"="gene",
#                                      "protein_change"="AAChange"))
# msk.protein.cancer <- msk %>%
#   dplyr::inner_join(dbRules,by=c("Hugo_Symbol"="gene",
#                                  "protein_change"="AAChange",
#                                  "CANCER_TYPE"="MskCancerType")) 
# 
# ### impose cancer type matches (genome)
# msk.genome.cancer <- msk %>%
#   dplyr::inner_join(dbRules,by=c("Chromosome.str"="chr",
#                                  "Start_Position"="pos",
#                                  "Tumor_Seq_Allele2"="alt",
#                                  "CANCER_TYPE"="MskCancerType"))
# 
