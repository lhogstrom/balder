library(dplyr)

inFile <- "/data/common/hartwig_pcgr/24_Skin/CPCT02100033T.pcgr.grch37.maf"
df.pcgr <- read.csv(inFile,sep="\t",skip=1) %>%
  dplyr::mutate(AAChange=stringr::str_replace(HGVSp_Short,"p.", "")) %>%
  data.frame()


# consider filtering on:
# gene name in oncokb list
# n_alt_count and n_ref_count

