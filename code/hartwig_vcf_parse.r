library(vcfR)

baseDir <- "/data/larsonh/hartwig"
#baseDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Hartwig/data" 
figDir <- paste0(baseDir,"/output")

inFile <- paste0(baseDir,"/samples2.txt")
sampleList <- read.csv(inFile,sep=",", header = F)
# load R in cancell - module load R/4.1.2-foss-2021b
# first need to unzip vcf files

outVariantTable <- data.frame()
for (sample in sampleList$V1) {
  print(sample)
  #inFile <- "/data/larsonh/hartwig/CPCT02060137T/purple/CPCT02060137T.purple.somatic.test.vcf"
  inFile <- paste0(baseDir,"/",sample,"/purple/",sample,".purple.somatic.vcf.gz")
  
  #local path
  #inFile <- paste0(baseDir,"/",sample,".purple.somatic.vcf.gz")
  print(inFile)
  vcf <- read.vcfR( inFile, verbose = FALSE )
  df <- vcfR::vcfR2tidy(vcf)
  
  df.fix <- df$fix
  dim(df.fix)
  
  df.gt <- df$gt
  dim(df.gt)
  
  ### check for matching entries
  #head(df.gt[,c("POS","gt_GT_alleles")])
  #head(df.fix[,c("POS","REF","ALT")])
  
  # parse redundant entries
  #all(df.fix[,"POS"] == df.gt[1:3095,"POS"])
  #g1 <- df.gt[1:3095,]
  #head(data.frame(g1))
  #g2 <- df.gt[3096:6190,]
  #head(data.frame(g2))
  #all(g1 == g2)
  
  # parse Ref and alt count entries 
  gtSize <- dim(df.gt)[[1]]
  g1 <- df.gt[1:(gtSize/2),]
  #head(data.frame(g1))
  g2 <- df.gt[((gtSize/2)+1):(gtSize),]
  #head(data.frame(g2))
 
  # scratch 
  #head(g2[,c("POS","gt_GT_alleles")])
  #head(df.fix[,c("POS","REF","ALT")])
  #head(g1[,c("gt_RDP")])
  #head(g2[,c("gt_RDP")])
  #head(g1[,c("gt_RAD")])
  #head(g2[,c("gt_RAD")])
  
  af1 <- g1[,c("gt_AF","gt_GT_alleles")]
  colnames(af1) <- c("Ref_AF","gt_REF_alleles")
  af2 <- g2[,c("gt_AF","gt_GT_alleles")]
  colnames(af2) <- c("Alt_AF","gt_alt_alleles")
  
  dd <- cbind(df.fix, af1, af2)
  dd$sample <- sample
  outCols <- c("sample","ChromKey","CHROM","POS","REF","ALT" ,"QUAL", "FILTER","BIALLELIC", "HOTSPOT","IMPACT",
               "Ref_AF","gt_REF_alleles","Alt_AF","gt_alt_alleles")
  
  outVariantTable <- rbind(outVariantTable,dd[,outCols])
  
  #df.meta <- df$meta
  #dim(df.meta)
  #head(data.frame(df.meta))

}

outFile <- paste0(baseDir,"/pancan_variant_table_subset.txt")
write.table(outVariantTable,outFile,sep="\t",row.names = F)

outVariantTable <- read.csv(outFile,sep="\t") # load data if previously generated


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
