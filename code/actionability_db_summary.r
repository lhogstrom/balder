library(dplyr)
library(tidyr)
library(ggplot2)
library(DBI)
library(RSQLite)

# load tmp data structures from 'actionability_db_curration.r'
figDir <- "../../output/actionability_db_curration_20231220"
inF <- paste0(figDir,"/cancerTypeMatching20231212.RData")
load(inF)

### MOA summary
print(table(moa$gene))
print(table(moa$predictive_implication))
print(table(moa$source_type))
print(table(moa$feature_type))
print(table(moa$variant_annotation))
print(table(moa$cosmic_signature_number,exclude=NULL))

signatures.moa <- moa[moa$feature_type=="Mutational signature",]

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


# count of alteration sources from CIViC or MOA for AA changes vs. other alterations
table(dbProtein[,"source"]) # except 
table(dbOtherAlt[,"source"])


