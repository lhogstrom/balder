library(vcfR)
library(pryr)
library(stringr)

baseDir <- "/data/larsonh/hartwig"
#baseDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Hartwig/data" 
figDir <- paste0(baseDir,"/output")

inFile <- paste0(baseDir,"/samples3.txt")
sampleList <- read.csv(inFile,sep=",", header = F)

inFile <- "/data/larsonh/hartwig/metadata.tsv"
hMeta <- read.csv(inFile,sep="\t")

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
  af2$tumor_AF_alternative2 <- af2$gt_RAD_2/af2$gt_tumor_DP
  
  dd <- cbind(df.fix, af1, af2)
  dd$sample <- sample
  outCols <- c("sample","ChromKey","CHROM","POS","REF","ALT" ,"QUAL", "FILTER","BIALLELIC",
               "Ref_AF","gt_REF_alleles","Ref_AF_alternative2","gt_REF_DP",
               "tumor_AF","gt_tumor_alleles","tumor_AF_alternative2","gt_tumor_DP")
  iPass <- dd$FILTER=="PASS"
  dd.filter <- dd[iPass,] # dd$gt_ALT_DP > 10
  #outVariantTable <- rbind(outVariantTable,dd.filter[,outCols])
  #print(pryr::object_size(outVariantTable))

  # test out vcf for pcgr  
  outvcf <- vcf
  infoOut <- paste0("TVAF=",as.character(dd.filter$tumor_AF),
                    ";TDP=",as.character(dd.filter$gt_tumor_DP),
                    ";CVAF=",as.character(dd.filter$Ref_AF),
                    ";CDP=",as.character(dd.filter$gt_REF_DP))
  
  # filter to only PASS entries
  outvcf@fix <- outvcf@fix[iPass,]
  outvcf@fix[,8] <- infoOut
  nVariants <- dim(dd.filter)[[1]]
  outvcf@gt <- matrix("",nVariants)
  
  ### remove FORMAT headers
  
  outvcf@meta <- c(outvcf@meta[1:28],
                  "##INFO=<ID=TVAF,Number=.,Type=Float,Description=\"Allelic fraction of alternative allele in tumor\">",
                  "##INFO=<ID=TDP,Number=.,Type=Integer,Description=\"Read depth across variant site in tumor\">",
                  "##INFO=<ID=CVAF,Number=.,Type=Float,Description=\"Allelic fraction of alternative allele in control\">",
                  "##INFO=<ID=CDP,Number=.,Type=Integer,Description=\"Read depth across variant site in control\">",
                   outvcf@meta[61:90])
  outFile <- paste0(baseDir,"/",sample,"/purple/",sample,".pcgr_format.vcf.gz")
  vcfR::write.vcf(outvcf,file=outFile)
  #
  pcgrDir <- "/data/sigven/pcgr"
  outDir <- paste0(baseDir,"/",sample,"/purple/testpcgr")
  str.pcgr <- paste0("pcgr --input_vcf ",outFile,
                     " --output_dir ", outDir,
                     " --pcgr_dir ", pcgrDir,
                     " --genome_assembly grch37",
                     " --sample_id ", sample,
                     " --tumor_dp_tag TDP --tumor_af_tag TVAF --basic")
  print(str.pcgr)
  
}

#outFile <- paste0(baseDir,"/pancan_variant_table_subset_filtered_v1.txt")
#write.table(outVariantTable,outFile,sep="\t",row.names = F)
#outVariantTable <- read.csv(outFile,sep="\t") # load data if previously generated





#   
# table(outVariantTable$HOTSPOT)
# table(outVariantTable$BIALLELIC)
# 
# outF <-  paste0(figDir,"/variants_alt_af.png")
# ggplot(outVariantTable,aes(x=Alt_AF))+
#   geom_histogram(alpha=.4)+
#   theme_bw()+
#   #facet_grid(source~.,scale="free_y",)+
#   xlab("Variant AF")+
#   ggtitle(paste0("100 Hartwig samples\n all ALT variants"))+
#   scale_fill_brewer(palette="Set1",drop=FALSE)+
#   theme(plot.title = element_text(hjust = 0.5))
# #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
# ggsave(outF,height = 5,width = 5)
# 
# 
# outF <-  paste0(figDir,"/variants_Ref_af_log.png")
# ggplot(outVariantTable,aes(x=Ref_AF))+
#   geom_histogram(alpha=.4)+
#   theme_bw()+
#   #facet_grid(source~.,scale="free_y",)+
#   xlab("Variant AF")+
#   scale_x_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   )+
#   ggtitle(paste0("100 Hartwig samples\n all REF variants"))+
#   scale_fill_brewer(palette="Set1",drop=FALSE)+
#   theme(plot.title = element_text(hjust = 0.5))
# #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
# ggsave(outF,height = 5,width = 5)
# 
# 
# ### number of subjects seen with a given variant
# variant.cnts <- outVariantTable %>%
#   group_by(CHROM,POS,REF,ALT) %>%
#   summarise(number.of.samples=dplyr::n()) %>%
#   dplyr::arrange(desc(number.of.samples))
# 
# #########################
# ### metadata analysis ###
# #########################
# 
# subject.sample.cnts <- hMeta %>%
#   dplyr::group_by(hmfPatientId) %>%
#   dplyr::summarise(n=dplyr::n())
# table(subject.sample.cnts$n)

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
