#!/software/software/easybuild/software/R/4.2.1-foss-2022a/bin/Rscript

library(dplyr)

# define inputs
tissueList <- "/data/larsonh/hartwig/output/tissue_list.txt"
tList <- read.csv(tissueList,sep="\t",header=F)
okbGeneList <- "/data/larsonh/hartwig/output/cancerGeneList.tsv"
oncokbGenes <- read.csv(okbGeneList,sep="\t")

baseDir <- "/data/common/hartwig_pcgr"

# out files
outMAF <- "/data/larsonh/hartwig/output/hartwig_sample_pcgr_maf_oncokb_gene_subset.maf" 
outOncokb <- "/data/larsonh/hartwig/output/hartwig_sample_pcgr_maf_oncokb_full.maf" 
outTbl <- "/data/larsonh/hartwig/output/hartwig_sample_pcgr_maf_oncokb_gene_subset_cnt_tbl.csv" 

# allowable VEP consequences
conseqList <- c("splice_acceptor_variant",
		"splice_donor_variant",
		"stop_gained",
		"frameshift_variant",
		"missense_variant",
		"frameshift_variant",
		"inframe_insertion",
		"inframe_deletion",
		"stop_gained",
		"stop_lost")

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
	  	df.oncokb <- df.pcgr[df.pcgr$Hugo_Symbol %in% oncokbGenes$Hugo.Symbol,] 
	  	df.subset <- df.oncokb %>%
			mutate(contains_vep_str = sapply(Consequence, function(x) any(sapply(conseqList, function(y) grepl(y, x))))) %>%
			dplyr::filter(HGVSp != "" | contains_vep_str) %>%
			dplyr::select(-contains_vep_str)
	  	print(dim(df.pcgr))
		print(".")
		print(dim(df.subset))
		print("----")
		write.table(df.subset, file = outMAF, append = TRUE, sep = "\t", quote =F, row.names = FALSE, col.names = !file.exists(outMAF))
		#write.table(df.oncokb, file = outOncokb, append = TRUE, sep = "\t", quote=F, row.names = FALSE, col.names = !file.exists(outOncokb))
		varSumRow <- data.frame(rawVariantCnt=c(dim(df.pcgr)[[1]]),
				oncokbVariantCnt=c(dim(df.oncokb)[[1]]),
				filteredVariantCnt=c(dim(df.subset)[[1]]),
				sample=sample_name)
		write.table(varSumRow, file = outTbl, append = TRUE, sep = "\t", quote=F, row.names = FALSE, col.names = !file.exists(outTbl))

	}
}

