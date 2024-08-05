library(dplyr)

tissueList <- "/data/larsonh/hartwig/output/tissue_list.txt"
tList <- read.csv(tissueList,sep="\t",header=F)
okbGeneList <- "/data/larsonh/hartwig/output/cancerGeneList.tsv"
oncokbGenes <- read.csv(okbGeneList,sep="\t")

baseDir <- "/data/common/hartwig_pcgr"

mafTable <- data.frame()
varCntTbl <- data.frame()
for (dir in tList$V1) {
	print(dir)
	dirPath <- paste0(baseDir,"/",dir)
	maf_files <- list.files(path = dirPath, pattern = "\\pcgr.grch37.maf$", full.names = TRUE) 
	for (f in maf_files) {
		# get sample name
		file_name <- basename(f)
		parts <- strsplit(file_name, "\\.")[[1]]
		sample_name <- parts[1]
		print(paste0(sample_name, ", ", dir))
		# read maf table
		df.pcgr <- read.csv(f,sep="\t",skip=1) %>%
		  dplyr::mutate(AAChange=stringr::str_replace(HGVSp_Short,"p.", ""),
				filename=f,
				sample=sample_name) %>%
 		  data.frame()
	  	# filter out non-coding and variants not in oncokb list
	  	df.subset <- df.pcgr[df.pcgr$Hugo_Symbol %in% oncokbGenes$Hugo.Symbol,] %>%
			dplyr::filter(HGVSp != "")
	  	print(dim(df.pcgr))
		print(".")
		print(dim(df.subset))
		print("----")
		varCntTbl <- rbind(varCntTbl,data.frame(rawVariantCnt=c(dim(df.pcgr)[[1]]),
				filteredVariantCnt=c(dim(df.subset)[[1]]),
				sample=sample_name))
		mafTable <- rbind(mafTable,df.subset)
	}
}
outF <- "/data/larsonh/hartwig/output/hartwig_sample_pcgr_maf_oncokb_gene_subset.maf" 
write.table(mafTable,outF,row.names=F,quote=F,sep="\t")
outF <- "/data/larsonh/hartwig/output/hartwig_sample_pcgr_maf_oncokb_gene_subset_cnt_tbl.csv" 
write.table(varCntTbl,outF,row.names=F,quote=F,sep=",")

