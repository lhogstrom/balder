#!/bin/bash

baseDir=/data/larsonh/hartwig
inFile=$baseDir/samples2.txt

while read p; do
	echo $p
	inVCFgz=$baseDir/${p}/purple/${p}.pcgr_format.vcf.gz
	inVCF=$baseDir/${p}/purple/${p}.pcgr_format.vcf
	gunzip $inVCFgz
	bgzip $inVCF
	tabix $inVCFgz
	pcgr --input_vcf /data/larsonh/hartwig/${p}/purple/${p}.pcgr_format.vcf.gz --output_dir /data/larsonh/hartwig/${p}/purple/testpcgr --pcgr_dir /data/sigven/pcgr --genome_assembly grch37 --sample_id $p --tumor_dp_tag TDP --tumor_af_tag TVAF --control_dp_tag CDP --control_af_tag CVAF --basic
done <$inFile

