library(vcfR)
library(pryr)

baseDir <- "/data/larsonh/hartwig"
#baseDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Hartwig/data" 
figDir <- paste0(baseDir,"/output")

inFile <- paste0(baseDir,"/samples2.txt")
sampleList <- read.csv(inFile,sep=",", header = F)
# load R in cancell - module load R/4.1.2-foss-2021b
# first need to unzip vcf files

outVariantTable <- data.frame()
for (sample in as.character(sampleList$V1)) {
  print(sample)
  #inFile <- "/data/larsonh/hartwig/CPCT02060137T/purple/CPCT02060137T.purple.somatic.test.vcf"
  inFile <- paste0(baseDir,"/",sample,"/purple/",sample,".purple.somatic.vcf.gz")
  
  #local path
  #inFile <- paste0(baseDir,"/",sample,".purple.somatic.vcf.gz")
  print(inFile)
  vcf <- read.vcfR( inFile, verbose = FALSE )
  df <- vcfR::vcfR2tidy(vcf)
  
  #df.meta <- df$meta
  #dim(df.meta)
  #head(data.frame(df.meta))
  
  df.fix <- df$fix
  print(dim(df.fix))
  
  df.gt <- df$gt
  #dim(df.gt)
  
  # parse Ref and alt count entries 
  gtSize <- dim(df.gt)[[1]]
  g1 <- df.gt[1:(gtSize/2),]
  g2 <- df.gt[((gtSize/2)+1):(gtSize),]
 
  af1 <- g1[,c("gt_AF","gt_GT_alleles","gt_DP","gt_RAD")]
  colnames(af1) <- c("Ref_AF","gt_REF_alleles","gt_REF_DP","gt_REF_RAD")
  af1$gt_RAD_1 <- as.numeric(sapply(strsplit(g1$gt_RAD, "\\,"), "[[", 1))
  af1$gt_RAD_2 <- as.numeric(sapply(strsplit(g1$gt_RAD, "\\,"), "[[", 2))
  #af1$Ref_AF_alternative <- af1$gt_RAD_1/af1$gt_REF_DP
  af1$Ref_AF_alternative2 <- af1$gt_RAD_2/af1$gt_REF_DP
  af2 <- g2[,c("gt_AF","gt_GT_alleles","gt_DP","gt_RAD")]
  colnames(af2) <- c("tumor_AF","gt_tumor_alleles","gt_tumor_DP","gt_tumor_RAD")
  af2$gt_RAD_1 <- as.numeric(sapply(strsplit(g2$gt_RAD, "\\,"), "[[", 1))
  af2$gt_RAD_2 <- as.numeric(sapply(strsplit(g2$gt_RAD, "\\,"), "[[", 2))
  #af2$Alt_AF_alternative <- af2$gt_RAD_1/af2$gt_ALT_DP
  af2$tumor_AF_alternative2 <- af2$gt_RAD_2/af2$gt_tumor_DP
  
  dd <- cbind(df.fix, af1, af2)
  dd$sample <- sample
  outCols <- c("sample","ChromKey","CHROM","POS","REF","ALT" ,"QUAL", "FILTER","BIALLELIC",
               "Ref_AF","gt_REF_alleles","Ref_AF_alternative2","gt_REF_DP",
               "tumor_AF","gt_tumor_alleles","tumor_AF_alternative2","gt_tumor_DP")
  dd.filter <- dd[dd$FILTER=="PASS",] # dd$gt_ALT_DP > 10
  print(dim(dd.filter))
  outVariantTable <- rbind(outVariantTable,dd.filter[,outCols])
  print(pryr::object_size(outVariantTable))
  
  #vcfR::write.vcf
}

outFile <- paste0(baseDir,"/pancan_variant_table_subset_filtered_v1.txt")
write.table(outVariantTable,outFile,sep="\t",row.names = F)
#outVariantTable <- read.csv(outFile,sep="\t") # load data if previously generated


### summary data
variant.cnts.per.subject <- outVariantTable %>%
  group_by(sample) %>%
  summarise(number.of.variants=dplyr::n())

outF <-  paste0(figDir,"/variants_per_subject.png")
ggplot(variant.cnts.per.subject,aes(x=number.of.variants))+
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
  ggplot(variant.cnts.per.subject,aes(x=number.of.variants))+
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

  
table(outVariantTable$HOTSPOT)
table(outVariantTable$BIALLELIC)

outF <-  paste0(figDir,"/variants_alt_af.png")
ggplot(outVariantTable,aes(x=Alt_AF))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(source~.,scale="free_y",)+
  xlab("Variant AF")+
  ggtitle(paste0("100 Hartwig samples\n all ALT variants"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
#theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggsave(outF,height = 5,width = 5)


outF <-  paste0(figDir,"/variants_Ref_af_log.png")
ggplot(outVariantTable,aes(x=Ref_AF))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(source~.,scale="free_y",)+
  xlab("Variant AF")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  ggtitle(paste0("100 Hartwig samples\n all REF variants"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
#theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggsave(outF,height = 5,width = 5)


### number of subjects seen with a given variant
variant.cnts <- outVariantTable %>%
  group_by(CHROM,POS,REF,ALT) %>%
  summarise(number.of.samples=dplyr::n()) %>%
  dplyr::arrange(desc(number.of.samples))


# ### scratch
# inFile <- "/data/larsonh/hartwig/CPCT02060137T/purple/CPCT02060137T.purple.somatic.test.vcf"
# vcf <- read.vcfR( inFile, verbose = FALSE )
# df <- vcfR::vcfR2tidy(vcf)
# 
# # sample has 4,499 variants
# df.fix <- df$fix
# head(df.fix$PURPLE_AF)
# head(df.fix$PON_MAX / df.fix$PON_COUNT)
# 
# # PURPLE_AF - "PURPLE_AF: Purity adjusted allelic frequency of variant"  
# 
# df.gt <- df$gt
# head(data.frame(df.gt))
# 
# df.meta <- df$meta
# head(data.frame(df.meta))
# 
# 
