# Databricks notebook source
#library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#remotes::install_github('sigven/geneOncoX')
#library(geneOncoX)

# COMMAND ----------
figDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/Pancan_biomarker_incidence_20221004"

inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_cna.txt"
cna.full <- read.csv(inFile,sep="\t")

inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/ICGC_TCGA_WGS_2020/CNA_Genes.txt"
cyto <- read.csv(inFile,sep="\t")

inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_mutations.txt"
sv <- read.csv(inFile,sep="\t",skip = 2)

inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/ICGC_TCGA_WGS_2020/pancan_pcawg_2020/data_clinical_sample.txt"
pancanSampInfo <- read.csv(inFile,sep="\t",skip = 4)

inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/MOAlmanac_43018_2021_243_MOESM2_ESM.txt"
moa <- read.csv(inFile,sep="\t")

# CIViC
inFile <- "/Users/larsonhogstrom/Documents/variant_annotation/CIViC-01-Dec-2021-ClinicalEvidenceSummaries.tsv"
clinical <- read.csv(inFile,sep="\t",quote = "")

# combined MOA+CIVIC file
inFile <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/mutation_incidence_20220913/civic_MOA_clinically_actionable_list.txt"
dbRules <- read.csv(inFile,sep="\t")

# load geneoOncoX basic
inFile <- "/Users/larsonhogstrom/Documents/code/oncology_biomarkers/geneOncoX_basic_records.txt"
geneOncox <- read.csv(inFile,sep="\t",quote = "",stringsAsFactors = FALSE)

# get cyto labels for each gene
cyto.label <- cyto %>%
  dplyr::group_by(Gene) %>%
  dplyr::filter(dplyr::row_number()==1)
cyto.label$chr_arm <- sapply(strsplit(as.character(cyto.label$Cytoband), "\\."), "[[", 1)
  
# pivot cna longer
cna <- cna.full %>% 
  tidyr::pivot_longer(!Hugo_Symbol, names_to = "sample", values_to = "value") %>%
  dplyr::left_join(cyto.label[,c("Gene","Cytoband","chr_arm")],by=c("Hugo_Symbol"="Gene")) %>%
  dplyr::left_join(pancanSampInfo,by=c("sample"="SAMPLE_ID"))
print(dim(cna))


## how many cyto labesl are there in CNA?
table(cna$Hugo_Symbol %in% cyto$Gene)
print(length(unique(cyto$Cytoband)))
print(length(unique(cna$chr_arm)))

## what is the distribution of CNA values at the gene-level?
# outF <-  paste0(figDir,"/CNA_value_gene_dist.png")
# ggplot(cna,aes(x=value))+
#   geom_histogram(alpha=.4)+
#   theme_bw()+
#   #facet_grid(source~.,scale="free_y",)+
#   xlab("CNA value")+
#   ggtitle(paste0("Observed gene-level CNA values in PANCAN"))+
#   scale_fill_brewer(palette="Set1",drop=FALSE)+
#   theme(plot.title = element_text(hjust = 0.5))+
#   #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   ggsave(outF,height = 5,width = 5)



### write out CNV info
iAmp <- clinical$variant == "Amplification" | clinical$variant == "AMPLIFICATION"
iDel <- clinical$variant == "Deletion" | clinical$variant == "DELETION"
clinical$direction <- "Amplification"
clinical[iDel,"direction"] <- "Deletion"
clinical$gene <- as.character(clinical$gene)
cn.clinical <- clinical[iAmp | iDel, ]
outFile <- paste0(figDir,"/CIViC_copy_number.txt")
write.table(cn.clinical,outFile,sep="\t",row.names = F)

table(as.character(cn.clinical$gene),as.character(cn.clinical$variant))
table(as.character(cn.clinical$clinical_significance))
table(as.character(cn.clinical$disease))


### what are the MOA copy number entries
moa.cna <- moa[moa$feature_type=="Copy Number",]
moa.cna$gene <- as.character(moa.cna$gene)

outFile <- paste0(figDir,"/MOA_copy_number.txt")
write.table(moa.cna,outFile,sep="\t",row.names = F)

table(moa.cna$direction)
table(as.character(moa.cna$gene))
table(as.character(moa.cna$gene),moa.cna$direction)
table(as.character(moa.cna$predictive_implication))
print(table(moa.cna$gene %in% cyto.label$Gene))
# combine CIVIC and MOA by gene and direction
cmoa <- rbind(cn.clinical[,c("gene","direction")],moa.cna[,c("gene","direction")])
table(as.character(cmoa$gene),cmoa$direction)

# pick a single direction for each gene
#moa.cna.rep <- moa.cna %>%
#  group_by(gene) %>%
#  filter(row_number()==1)
#table(moa.cna$gene,moa.cna$direction)

cmoa.cna.rep <- cmoa %>%
  group_by(gene) %>%
  filter(row_number()==1)
table(cmoa$gene,cmoa$direction)

cna.actionable <- cna %>% 
  filter(Hugo_Symbol %in% cmoa.cna.rep$gene) %>%
  left_join(cmoa.cna.rep[,c("gene","direction")],by=c("Hugo_Symbol"="gene"))
print(dim(cna.actionable))

outF <-  paste0(figDir,"/CNA_actionable_by_gene_v1.png")
ggplot(cna.actionable,aes(x=value))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  facet_grid(Hugo_Symbol~direction,scale="free_y",)+
  xlab("mean CNA value")+
  ggtitle(paste0("Observed crhom-arm CNA values in PANCAN"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 35,width = 5)

# get perc of subjects falling into each bin across genes
nSamp <- length(unique(cna$sample))
full.cna.bin.cnt <- cna %>%
  #dplyr::filter(!is.na(perc_subj)) %>%
  dplyr::group_by(value) %>% 
  dplyr::summarise(bin_cnts=n(),
                   n_genes=n_distinct(Hugo_Symbol),
                   av_perc_subj_across_genes=100*(bin_cnts/(nSamp*n_genes)))

### make cna value bins pas perc of subjects
cna.act.perc <- cna.actionable %>%
  group_by(Hugo_Symbol,value) %>%
  summarise(bin_cnts=n(),
            perc_subj=100*(bin_cnts/nSamp)) %>%
  dplyr::left_join(full.cna.bin.cnt[,c("value","av_perc_subj_across_genes")],by="value") %>%
  left_join(cmoa.cna.rep[,c("gene","direction")],by=c("Hugo_Symbol"="gene"))
cna.act.perc$relative_change <- cna.act.perc$perc_subj - cna.act.perc$av_perc_subj_across_genes

cna.act.perc.cType <- cna.actionable %>%
  group_by(CANCER_TYPE,Hugo_Symbol,value) %>%
  summarise(bin_cnts=n(),
            perc_subj=100*(bin_cnts/nSamp),
            nSubjects=dplyr::n_distinct(PATIENT_ID)) %>%
  dplyr::left_join(full.cna.bin.cnt[,c("value","av_perc_subj_across_genes")],by="value") %>%
  left_join(cmoa.cna.rep[,c("gene","direction")],by=c("Hugo_Symbol"="gene"))
cna.act.perc.cType$relative_change <- cna.act.perc.cType$perc_subj - cna.act.perc.cType$av_perc_subj_across_genes



# check they add up to 100%
#tmp.perc <- cna.act.perc %>%
#  filter(!is.na(perc_subj)) %>%
#  group_by(Hugo_Symbol) %>%
#  summarise(total=sum(perc_subj))
#head(tmp.perc)

outF <-  paste0(figDir,"/CNA_actionable_by_gene_v2.png")
ggplot(cna.actionable,aes(x=value))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  facet_grid(Hugo_Symbol~direction,scale="free_y",)+
  xlab("mean CNA value")+
  ggtitle(paste0("Observed crhom-arm CNA values in PANCAN"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 35,width = 5)

outF <-  paste0(figDir,"/CNA_differential_actionable_by_gene.png")
ggplot(cna.act.perc,aes(x=value,y=relative_change))+
  geom_point(alpha=.4)+
  theme_bw()+
  facet_grid(Hugo_Symbol~direction,scale="free_y",)+
  geom_hline(yintercept = 0,linetype="dotted")+
  xlab("mean CNA value")+
  ggtitle(paste0(""))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 35,width = 5)

### what are the biggest changes in relative change
largest.relative.change <- cna.act.perc %>%
  dplyr::arrange(desc(relative_change))
largestChangeGenes <- largest.relative.change$Hugo_Symbol[1:10]
outFile <- paste0(figDir,"/largest_relative_CNA_change.txt")
write.table(largest.relative.change,outFile,sep="\t",row.names = F)


outF <-  paste0(figDir,"/CNA_differential_actionable_by_gene_most_change.png")
ggplot(cna.act.perc[cna.act.perc$Hugo_Symbol %in% largestChangeGenes,],aes(x=value,y=relative_change))+
  geom_point(alpha=.4)+
  theme_bw()+
  facet_grid(Hugo_Symbol~direction,scale="free_y",)+
  geom_hline(yintercept = 0,linetype="dotted")+
  xlab("mean CNA value")+
  ggtitle(paste0(""))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 15,width = 5)


### make cna value bins for all genes
nSamp <- length(unique(cna$sample))
cna.perc <- cna %>%
  group_by(Hugo_Symbol,value) %>%
  summarise(bin_cnts=n(),
            perc_subj=100*(bin_cnts/nSamp)) %>%
  dplyr::left_join(full.cna.bin.cnt[,c("value","av_perc_subj_across_genes")],by="value") %>%
  left_join(cmoa.cna.rep[,c("gene","direction")],by=c("Hugo_Symbol"="gene")) %>%
  dplyr::left_join(geneOncox,by=c("Hugo_Symbol"="symbol"))
cna.perc$relative_change <- cna.perc$perc_subj - cna.perc$av_perc_subj_across_genes
outF <-  paste0(figDir,"/CNA_differential_actionable_by_all_genes.png")
ggplot(cna.perc,aes(x=relative_change))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  facet_grid(direction~value,scale="free_y",)+
  geom_vline(xintercept = 0,linetype="dotted")+
  xlab("perc differential")+
  ylab("number of genes")+
  ggtitle(paste0("Percent of subject enrichment for CNA events \n according to clinical actionability"))+
  xlim(-25,25)+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 4,width = 10)

summary.cna.perc <- cna.perc[,c("Hugo_Symbol","value","relative_change")] %>%
  tidyr::pivot_wider(names_from = value,values_from=relative_change)
colnames(summary.cna.perc) <- c("Hugo_Symbol","neg2","neg1","zero","pos1","pos2","NA")
summary.cna.perc <- summary.cna.perc %>%
  dplyr::arrange(desc(pos2)) %>%
  left_join(cmoa.cna.rep[,c("gene","direction")],by=c("Hugo_Symbol"="gene")) %>%
  dplyr::left_join(geneOncox,by=c("Hugo_Symbol"="symbol")) %>%
  data.frame()

summary.cna.perc$ncg_onc_tsg <- "null"
summary.cna.perc[summary.cna.perc$ncg_tsg==T & !is.na(summary.cna.perc$cgc_tsg),"ncg_onc_tsg"] <- "TSG"
summary.cna.perc[summary.cna.perc$ncg_oncogene==T & !is.na(summary.cna.perc$cgc_oncogene),"ncg_onc_tsg"] <- "ONC"

summary.cna.perc$cgc_onc_tsg <- "null"
summary.cna.perc[summary.cna.perc$cgc_tsg==T & !is.na(summary.cna.perc$cgc_tsg),"cgc_onc_tsg"] <- "TSG"
summary.cna.perc[summary.cna.perc$cgc_oncogene==T & !is.na(summary.cna.perc$cgc_oncogene),"cgc_onc_tsg"] <- "ONC"

my_colour = list(
  direction = c(Amplification = "#5977ff", Deletion = "#f74747"),
  cgc_oncogene = c(T = "#82ed82", F = "#9e82ed")
)

cna.perc.subset <- summary.cna.perc[!is.na(summary.cna.perc$direction),]
heatmapFile=paste0(figDir,"/cna_clustering3.pdf")
gcdfOut <- pheatmap(as.matrix((cna.perc.subset[,c("neg2","neg1","zero","pos1","pos2")])),
                    #annotation_row = cna.perc.subset[,c("direction","cgc_onc_tsg","gene_indispensability_pred","sanchezvega2018_signaling_pathway")], 
                    annotation_row = cna.perc.subset[,c("direction","cgc_onc_tsg")], 
                    #nnotation_colors = my_colour,
                    show_rownames = T,
                    labels_row = data.frame(cna.perc.subset)[,c("Hugo_Symbol")],
                    show_colnames = T,#labels_col = meta[gnums,"name"],
                    color = rev(brewer.pal(8,"GnBu")),fontsize_col = 5,
                    cluster_cols=F,cluster_rows=F,
                    filename = heatmapFile,width = 9, height = 9)

outF <-  paste0(figDir,"/CNA_differential_gene_annotation.png")
#ggplot(summary.cna.perc[summary.cna.perc$cgc_onc_tsg=="null",],aes(x=neg1,y=pos2,color=cgc_onc_tsg))+
ggplot()+
  geom_point(data=summary.cna.perc[summary.cna.perc$cgc_onc_tsg=="null",],aes(x=neg1,y=pos2,color=cgc_onc_tsg),alpha=.2)+
  geom_point(data=summary.cna.perc[!summary.cna.perc$cgc_onc_tsg=="null",],aes(x=neg1,y=pos2,color=cgc_onc_tsg),alpha=.8)
  #theme_bw()+
  #facet_grid(direction~value,scale="free_y",)+
  #geom_vline(xintercept = 0,linetype="dotted")+
  #xlab("perc differential")+
  #ylab("number of genes")+
  #ggtitle(paste0("Percent of subject enrichment for CNA events \n according to clinical actionability"))+
  #xlim(-25,25)+
  #scale_fill_brewer(palette="Set1",drop=FALSE)+
  #theme(plot.title = element_text(hjust = 0.5))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 4,width = 10)

outF <-  paste0(figDir,"/CNA_differential_gene_gene_indispensability_score.png")
ggplot()+
  geom_point(data=summary.cna.perc,aes(x=neg1,y=pos2,color=gene_indispensability_score),alpha=.2)
  ggsave(outF,height = 4,width = 10)

outF <-  paste0(figDir,"/CNA_differential_prob_exac_lof_intolerant.png")
ggplot()+
  geom_point(data=summary.cna.perc,aes(x=neg1,y=pos2,color=prob_exac_lof_intolerant),alpha=.2)
  ggsave(outF,height = 4,width = 10)

#### look at perc enrichment grouped by cancer types ####

# get per cance type gene enrichements
cna.perc.cType <- cna %>%
  dplyr::filter(Hugo_Symbol %in% cmoa.cna.rep$gene) %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::mutate(nSampC=dplyr::n_distinct(sample))  %>%
  dplyr::ungroup() %>%
  group_by(CANCER_TYPE,Hugo_Symbol,value) %>%
  summarise(bin_cnts=n(),
            perc_subj=100*(bin_cnts/max(nSampC))) %>%
  dplyr::left_join(full.cna.bin.cnt[,c("value","av_perc_subj_across_genes")],by="value") %>% # should we look at average bins by cancer type or across all cancers? 
  left_join(cmoa.cna.rep[,c("gene","direction")],by=c("Hugo_Symbol"="gene")) #%>%
  #dplyr::left_join(geneOncox,by=c("Hugo_Symbol"="symbol"))
cna.perc.cType$relative_change <- cna.perc.cType$perc_subj - cna.perc.cType$av_perc_subj_across_genes

summary.cna.perc.cType <- cna.perc.cType[,c("CANCER_TYPE","Hugo_Symbol","value","relative_change")] %>%
  tidyr::pivot_wider(names_from = value,values_from=relative_change)
colnames(summary.cna.perc.cType) <- c("CANCER_TYPE","Hugo_Symbol","neg2","neg1","zero","pos1","pos2","NA")
summary.cna.perc.cType <- summary.cna.perc.cType %>%
  dplyr::arrange(CANCER_TYPE,desc(pos2)) %>%
  left_join(cmoa.cna.rep[,c("gene","direction")],by=c("Hugo_Symbol"="gene")) %>%
  #dplyr::left_join(geneOncox,by=c("Hugo_Symbol"="symbol")) %>%
  data.frame()

### make a heatmap that shows for each cancer type what are the top 3 most deleted or amplified genes? 

cna.perc.cType.subset <- summary.cna.perc.cType[!is.na(summary.cna.perc.cType$direction),] %>%
  dplyr::arrange(CANCER_TYPE,desc(neg2)) %>%
  dplyr::group_by(CANCER_TYPE) %>%
  dplyr::filter(row_number() %in% c(1,2,3)) %>%
  data.frame()
rownames(cna.perc.cType.subset) <- paste0(cna.perc.cType.subset$CANCER_TYPE, " , ",cna.perc.cType.subset$Hugo_Symbol)
heatmapFile=paste0(figDir,"/cna_cancer_type_clustering.pdf")
gcdfOut <- pheatmap(as.matrix((cna.perc.cType.subset[,c("neg2","neg1","zero","pos1","pos2")])),
                    #annotation_row = cna.perc.subset[,c("direction","cgc_onc_tsg","gene_indispensability_pred","sanchezvega2018_signaling_pathway")], 
                    #annotation_row = cna.perc.cType.subset[,c("direction")], 
                    #nnotation_colors = my_colour,
                    show_rownames = T,
                    #labels_row = data.frame(rownames(cna.perc.cType.subset)),
                    show_colnames = T,#labels_col = meta[gnums,"name"],
                    color = rev(brewer.pal(8,"GnBu")),fontsize_col = 5,
                    cluster_cols=F,cluster_rows=F,
                    filename = heatmapFile,width = 11, height = 14)


# loop over genes
# moa.df <- data.frame(moa.cna.rep)
# for (gene in unique(moa.cna$gene)) {
#   print(gene)
#   cna.gene <- cna.perc[cna.perc$Hugo_Symbol==gene,]
#   direct <- as.character(moa.df[moa.df$gene==gene,"direction"])
#   iDirection <- cna.gene$direction == direct
#   iNa <- is.na(cna.gene$direction)
#   df.subset <- cna.gene[iDirection & iNa,]
#   table(df.subset$direction,exclude=NULL)
# }

# loop over bins and directions
for (bin in c(-2,-1,0,1,2)) {
  for (direction in c("Amplification","Deletion")) {
    print(paste0(bin," and ", direction))
    iNullBin <- is.na(cna.perc$value)
    iBin <- cna.perc$value==bin
    iDir <- cna.perc$direction==direction
    iNA <- is.na(cna.perc$direction)
    cna.match <- data.frame(cna.perc[(!iNullBin & iBin) & (iNA | iDir),])
    cna.match$direction2 <- as.character(cna.match$direction)
    cna.match[is.na(cna.match$direction2),"direction2"] <- "background"
    table(cna.match$direction2)
    cna.match$direction2 <- factor(cna.match$direction2)
    print(table(cna.match$direction2))
    pval.greater <- wilcox.test(relative_change~direction2,
                        data=cna.match,
                        alternative = "greater",
                        paired=FALSE)$p.value
    pval.less <- wilcox.test(relative_change~direction2,
                                data=cna.match,
                                alternative = "less",
                                paired=FALSE)$p.value
    pval.twosided <- wilcox.test(relative_change~direction2,
                                data=cna.match,
                                alternative = "two.sided",
                                paired=FALSE)$p.value
    #pval <- wilcox.test(cna.match$av_perc_subj_across_genes~cna.match$direction2)$p.value
    #print(pval.greater)
    #print(pval.less)
    print(pval.twosided)
  }
}

#####################
### short variants ###
#####################

# exclude events that aren't short variants
#exList <- c("expression","EXPRESSION","TRANSLOCATION","omoter","OMOTER","Copy Number Variation")
exList <- c("xpression|EXPRESSION|TRANSLOCATION|omoter|OMOTER|Copy Number Variation|mplification|ATION|REARRANGEMENT|DELETION|eletion")

iNotSV <- grepl(exList,dbRules$AAChange)
table(iNotSV)
genome.change.cnts <- dbRules[!iNotSV,] %>%
  dplyr::group_by(pos,ref,alt) %>%
  dplyr::summarise(n=n(),
                   nSource=dplyr::n_distinct(source),
                   nCancers=dplyr::n_distinct(Disease)) %>%
  dplyr::arrange(desc(n))

aaChange.cnts <- dbRules[!iNotSV,] %>%
  dplyr::group_by(gene, AAChange) %>%
  dplyr::summarise(n=n(),
                   nSource=dplyr::n_distinct(source),
                   nCancers=dplyr::n_distinct(Disease)) %>%
  dplyr::arrange(desc(n))

# select representative genome entry for general matches
genome.rep <- dbRules[!iNotSV,] %>%
  dplyr::group_by(chr,pos,alt) %>%
  dplyr::filter(dplyr::row_number()==1)

print(dim(sv))
sv.full <- sv %>%
  dplyr::left_join(pancanSampInfo,by=c("Tumor_Sample_Barcode"="SAMPLE_ID")) %>%
  dplyr::mutate(HGVSp_Short_mod=gsub("p.","",HGVSp_Short),
                Chromosome.str=paste0("chr",as.character(Chromosome)))
print(dim(sv.full))

#sv.full$MAF <- 100*(sv.full$t_alt_count / (sv.full$t_ref_count+sv.full$t_alt_count))


### subset pancan SV to actionable mutations by AA (in at least once cancer)
sv.aa.actionable <- sv.full %>% 
  filter(Hugo_Symbol %in% unique(as.character(aaChange.cnts$gene))) %>% # aa.rep
  inner_join(aaChange.cnts[,c("gene","AAChange")],by=c("Hugo_Symbol"="gene","HGVSp_Short_mod"="AAChange"))
print(dim(sv.aa.actionable))
aaPancanActionableCnts <- sv.aa.actionable %>%
  dplyr::group_by(Hugo_Symbol,HGVSp_Short_mod) %>%
  dplyr::summarise(nSubjects=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(nSubjects))
print(dim(aaPancanActionableCnts))
outFile <- paste0(figDir,"/top_actionable_pancan_matuations_by_AA_change.txt")
write.table(aaPancanActionableCnts,outFile,sep="\t",row.names = F)

### subset pancan SV to actionable mutations by genome change (in at least once cancer)
sv.genome.actionable <- sv.full %>% 
  filter(Hugo_Symbol %in% unique(as.character(genome.rep$gene))) %>%
  dplyr::inner_join(genome.rep,by=c("Chromosome.str"="chr",
                                 "Start_Position"="pos",
                                 "Tumor_Seq_Allele2"="alt"))
                                 #"CANCER_TYPE"="MskCancerType"))
print(dim(sv.genome.actionable))
genomePancanActionableCnts <- sv.genome.actionable %>%
  dplyr::group_by(Hugo_Symbol,Chromosome,Start_Position,Tumor_Seq_Allele2,HGVSp_Short_mod) %>%
  dplyr::summarise(nSubjects=dplyr::n_distinct(Tumor_Sample_Barcode)) %>%
  dplyr::arrange(desc(nSubjects))
print(dim(genomePancanActionableCnts))
outFile <- paste0(figDir,"/top_actionable_pancan_matuations_by_genomic_change.txt")
write.table(genomePancanActionableCnts,outFile,sep="\t",row.names = F)

### what is the MAF of matches compared to all other mutations? 
aaConcat <- paste0(sv.aa.actionable$Hugo_Symbol,"-",sv.aa.actionable$HGVSp_Short_mod)
length(unique(aaConcat))
table(aaConcat)

### how do cancer types line up for PANCAN vs. MOA/CIVIC
pancanCancerTypes <- as.character(unique(pancanSampInfo$CANCER_TYPE))
civicCancerTypes <- as.character(unique(clinical$disease))
moaCancerTypes <- as.character(unique(moa$disease))
mskTypes <- as.character(unique(dbRules$MskCancerType))

### which MSK types are not in pancan labels? 
mskTypes[!mskTypes %in% pancanCancerTypes]
moaCancerTypes[!moaCancerTypes %in% pancanCancerTypes]
civicCancerTypes[!civicCancerTypes %in% pancanCancerTypes]

########################################################
### compare co-occurence of mutations and cna evnets ###
########################################################

cna.actionable.gene <- cna %>%
  dplyr::filter(Hugo_Symbol %in% cmoa.cna.rep$gene) %>%
  dplyr::filter(value == 2 | value == -2)
length(unique(cna.actionable.gene$Hugo_Symbol))
length(unique(cna.actionable.gene$ICGC_SAMPLE_ID))
table(as.character(cna.actionable.gene$Hugo_Symbol))
table(cna.actionable.gene$value)

### join with SV
co.occur.tbl <- cna.actionable.gene %>%
  full_join(sv.aa.actionable,by="PATIENT_ID")
  
### create list of genes that occur most frequently together
gene.co.occr.cnts <- co.occur.tbl %>%
  filter(!is.na(HGVSp_Short_mod)) %>%
  dplyr::group_by(Hugo_Symbol.x,value,Hugo_Symbol.y,HGVSp_Short_mod) %>%
  summarize(n=n()) %>%
  dplyr::arrange(desc(n))
outFile <- paste0(figDir,"/top_co_occurance_results.txt")
write.table(gene.co.occr.cnts,outFile,sep="\t",row.names = F)


## MAF of cancer type specific actionabile variants 

### MAF results
maf.df.pancan <- data.frame(MAF=sv.full$DNA_VAF)
maf.df.pancan$source <- "PANCAN all variants"
maf.df.pancan$sample_type <- sv.full$SAMPLE_TYPE
maf.df.pancan$Tumor_Sample_Barcode <- sv.full$Tumor_Sample_Barcode

# maf of clinically actionable based on protein match
maf.df.protein<- data.frame(MAF=sv.aa.actionable$DNA_VAF)
maf.df.protein$source <- "clinically actionable"
maf.df.protein$sample_type <- sv.aa.actionable[,"SAMPLE_TYPE"]
maf.df.protein$Tumor_Sample_Barcode <- sv.aa.actionable$Tumor_Sample_Barcode

maf.df <- rbind(maf.df.pancan,maf.df.protein)

# plot clinically actionable distribution
outF <-  paste0(figDir,"/MAF_distribution_clinically_actionable.png")
ggplot(maf.df,aes(x=MAF,fill=source))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  facet_grid(source~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n in PANCAN"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 5,width = 5)

outF <-  paste0(figDir,"/MAF_distribution_clinically_actionable_density.png")
ggplot(maf.df,aes(x=MAF,fill=source))+
  geom_density(alpha=.4)+
  theme_bw()+
  #facet_grid(source~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n in PANCAN"))+
  scale_fill_brewer(palette="Set2",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 5,width = 5)

#Density plot of primary vs. met 
outF <-  paste0(figDir,"/MAF_distribution_clinically_actionable_density.png")
iSubset <- !maf.df$sample_type=="Cell line" & !is.na(maf.df$sample_type)
ggplot(maf.df[iSubset,],aes(x=MAF,fill=source))+
  geom_density(alpha=.4)+
  theme_bw()+
  facet_grid(sample_type~.,scale="free_y",)+
  xlab("MAF (%)")+
  ggtitle(paste0("MAF of clinically actioanable variants \n in PANCAN"))+
  scale_fill_brewer(palette="Set2",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(outF,height = 5,width = 5)


### more descriptive info

### how do distributions change when grouped by cyto band
cna.band <- cna %>%
  dplyr::group_by(Cytoband,sample) %>%
  dplyr::summarize(mean_value=mean(value))
print(head(cna.band))
print(dim(cna.band))

## what is the distribution of CNA values at the gene-level?
outF <-  paste0(figDir,"/CNA_value_cytoband_dist.png")
ggplot(cna.band,aes(x=mean_value))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(source~.,scale="free_y",)+
  xlab("mean CNA value")+
  ggtitle(paste0("Observed cytoband-level CNA values in PANCAN"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
#theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggsave(outF,height = 5,width = 5)

### check non-null values
iSampleNull <- is.na(cna$sample)
iBandNull <- is.na(cna$Cytoband)
table(iBandNull)
# are the entries the same across bands for the same sample 
table(cna[cna$Cytoband=="10p11.1" & cna$sample=="SP10084" & !iBandNull,"value"])
table(cna[cna$Cytoband=="6q27" & cna$sample=="SP10084" & !iBandNull,"value"])
table(cna[cna$Cytoband=="7p14.2" & cna$sample=="SP10084" & !iBandNull,"value"])
table(cna[cna$chr_arm=="14q13" & cna$sample=="SP10084" & !iBandNull,"value"])

### how do distributions change when grouped by chr arm
cna.arm <- cna %>%
  dplyr::group_by(chr_arm,sample) %>%
  dplyr::summarize(mean_value=mean(value))
print(head(cna.arm))
print(dim(cna.arm))
## what is the distribution of CNA values at the gene-level?
outF <-  paste0(figDir,"/CNA_value_chr_arm_dist.png")
ggplot(cna.arm,aes(x=mean_value))+
  geom_histogram(alpha=.4)+
  theme_bw()+
  #facet_grid(source~.,scale="free_y",)+
  xlab("mean CNA value")+
  ggtitle(paste0("Observed crhom-arm CNA values in PANCAN"))+
  scale_fill_brewer(palette="Set1",drop=FALSE)+
  theme(plot.title = element_text(hjust = 0.5))
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggsave(outF,height = 5,width = 5)

outRFile <- paste0(figDir,"/pancan_work_space.RData")
save.image(file = outRFile)


