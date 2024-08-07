#!/software/software/easybuild/software/R/4.2.1-foss-2022a/bin/Rscript

#This script performs the following actions:
#- reads in Hartwig sample metadata - maps primaryTumorLocation to the PCGR primary site
#- reads in a list of all somatic (SNV/InDel) VCF filenames provided by Hartwig (Purple) - and from this pulls out the tumor sample name
#- For each tumor sample:
#  - reads in copy number data and purity/ploidy (Purple) information from Hartwig
#  - converts the somatic VCF provided by Hartwig to TSV with vcf2tsvpy (easier to work with TSV files in R - data.frames)
#    - ignore LOW_CONFIDENCE calls
#    - creates the MAF and DP for control and tumor sample (i.e. TVAF, TDP, CVAF, CDP)
#  - prints a simplified “PCGR-ready” VCF per sample
#  - makes a PCGR command, running with Singularity/Apptainer (many paths (path to refdata, VEP etc) hard-coded)

pcgr_tt_codes <- read.table(file="pcgr_tt_codes.tsv", sep="\t", header = T)

metadata <- "metadata.tsv"

sample_data <- read.table(metadata, 
                          header = T, 
                          sep = "\t") |>
  dplyr::select(sampleId, primaryTumorLocation) |>
  dplyr::rename(tumor_sample_id = sampleId, site = primaryTumorLocation) |>
  dplyr::mutate(site = dplyr::case_when(
    site == "Ovary" | site == "Fallopian tube" ~ "Ovary/Fallopian Tube",
    site == "Colorectum" | site == "Anus" ~ "Colon/Rectum",
    site == "Stomach" | site == "Small intestine" | 
      site == "Gastroesophageal" | site == "Appendix" |
      site == "Esophagus" |
      site == "Gastrointestinal tract" ~ "Esophagus/Stomach",
    site == "Lymphoid tissue" ~ "Lymphoid",
    site == "Adrenal gland" ~ "Adrenal Gland",
    site == "Bile duct" | site == "Gallbladder" ~ "Biliary Tract",
    site == "Bone marrow" ~ "Myeloid",
    site == "Lung" | site == "Trachea" | site == "Bronchus" ~ "Lung",
    site == "Breast" ~ "Breast",
    site == "Nervous system" ~ "CNS/Brain",
    site == "Prostate" ~ "Prostate",
    site == "Pancreas" ~ "Pancreas",
    site == "Mesothelium" ~ "Pleura",
    site == "Vagina" | site == "Vulva" ~ "Vulva/Vagina",
    site == "Kidney" ~ "Kidney",
    site == "Eye" ~ "Eye",
    site == "Uterus" ~ "Uterus",
    site == "Skin" ~ "Skin",
    site == "Testis" ~ "Testis",
    site == "Liver" | site == "Hepatobiliary system" ~ "Liver",
    site == "Head and neck" ~ "Head and Neck",
    site == "Thyroid gland" ~ "Thyroid",
    site == "Thymus" ~ "Thymus",
    site == "Bone/Soft tissue" ~ "Soft Tissue",
    site == "Urothelial tract" ~ "Bladder/Urinary Tract",
    TRUE ~ "Any"
  )) |>
  dplyr::left_join(
    pcgr_tt_codes, by = "site"
  )

vcf_file_list <- "vcf_files_hartwig.tsv"
vcf_files <- read.table(vcf_file_list, sep="\t", header = T)
i <- 1

while(i <= NROW(vcf_files)){
  vep_buffer_size <- 500
  vcf_purple_fname <- vcf_files[i,"fname"]
  tumor_sample_id <- stringr::str_split(vcf_purple_fname,"/")[[1]][5]
  cna_purple_fname <- file.path("/data", "common", "hartwig",
                         tumor_sample_id, "purple",
                         paste0(tumor_sample_id,".purple.cnv.somatic.tsv"))
  pp_purple_fname <- file.path("/data", "common", "hartwig",
                        tumor_sample_id, "purple",
                        paste0(tumor_sample_id,".purple.purity.tsv"))
  vcf2tsv_fname <- file.path(here::here(), 'vcf2tsv', paste0(tumor_sample_id,".vcf2tsv.tsv"))
  vcf_fname <- file.path(here::here(), 'pcgr_input', paste0(tumor_sample_id,".vcf"))
  cna_fname <- file.path(here::here(), 'pcgr_input', paste0(tumor_sample_id,".cna.tsv"))
  purity <- NA
  ploidy <- NA
  
  tsite_code <- 0
  if(tumor_sample_id %in% sample_data$tumor_sample_id){
    tsite_code <- sample_data[sample_data$tumor_sample_id == tumor_sample_id, "pcgr_site_code"]
  }
 
  tsite_folder <- "0_Any"
  if(tumor_sample_id %in% sample_data$tumor_sample_id){
    tsite_folder <- sample_data[sample_data$tumor_sample_id == tumor_sample_id, "pcgr_site_folder"]
  } 

  if(!(dir.exists(file.path(here::here(),"pcgr_output",tsite_folder)))){
     dir.create(file.path(here::here(),"pcgr_output",tsite_folder))
  }

  if(file.exists(cna_purple_fname) &
     !file.exists(cna_fname)){

    cna_data <- read.table(
      cna_purple_fname, 
      comment.char = "#", header = T) |>
      dplyr::mutate(
        nMinor = round(minorAlleleCopyNumber),
        nMajor = round(majorAlleleCopyNumber),
        Chromosome = chromosome,
        Start = start,
        End = end) |>
      dplyr::select(
        Chromosome, Start, End, 
        nMajor, nMinor
      ) |>
      dplyr::distinct()
    
    options(scipen=999)
    
    write.table(cna_data, 
                cna_fname,
                quote = F, 
                row.names = F, 
                col.names = T, 
                sep = "\t")
  }
  
  if(file.exists(pp_purple_fname)){
    purity_ploidy <- read.table(
      pp_purple_fname, 
      comment.char = "#", header = T) |>
      dplyr::mutate(
        purity = round(purity, 4),
        ploidy = round(ploidy, 4)
      ) 
    
    purity <- purity_ploidy$purity
    ploidy <- purity_ploidy$ploidy
  }
  
  
  ## convert to TSV with vcf2tsvpy
  command = paste0("vcf2tsvpy --input_vcf ",
                   vcf_purple_fname," --out_tsv ", 
                   vcf2tsv_fname," --compress")
  
  if(!file.exists(paste0(vcf2tsv_fname,".gz")) & 
     file.exists(vcf_purple_fname)){
    system(command)
  }
  
  if(file.exists(paste0(vcf2tsv_fname,".gz")) &
     !file.exists(paste0(vcf_fname,".gz"))){
    
    ## read in allelic support for tumor sample
    tumor_data <- read.table(paste0(vcf2tsv_fname,".gz"), 
                             comment.char = "#", header = T) |> 
      dplyr::filter(VCF_SAMPLE_ID == tumor_sample_id) |> 
      dplyr::filter(TIER != "LOW_CONFIDENCE")
    
    ## read in allelic support for control sample
    control_data <- read.table(paste0(vcf2tsv_fname,".gz"), 
                               comment.char = "#", header = T) |> 
      dplyr::filter(VCF_SAMPLE_ID != tumor_sample_id) |> 
      dplyr::filter(TIER != "LOW_CONFIDENCE")
    
    ## join tumor and control data
    if(NROW(tumor_data) > 0 & 
       NROW(control_data) > 0 &
       NROW(tumor_data) == NROW(control_data)){
      control_data <- control_data |>
        dplyr::mutate(CVAF = AF, CDP = DP) |>
        dplyr::select(CHROM, POS, REF, ALT, CVAF, 
                      CDP)
      
      tumor_data <- tumor_data |>
        dplyr::mutate(TVAF = AF, TDP = DP) |>
        dplyr::select(CHROM, POS, REF, ALT, TVAF, 
                      TDP, VCF_SAMPLE_ID, TIER,
                      PURPLE_AF, PURPLE_CN) |>
        dplyr::rename(PURPLE_CALL_TIER = TIER)
      
      tumor_data <- tumor_data |>
        dplyr::mutate(PURPLE_CALL_TIER = dplyr::case_when(
          PURPLE_CALL_TIER == "HIGH_CONFIDENCE" ~ "HC",
          TRUE ~ as.character(PURPLE_CALL_TIER))
        ) |>
        dplyr::left_join(control_data, 
                         by = c("CHROM", "POS", "REF", "ALT")) |>
        dplyr::mutate(
          TVAF = paste0("TVAF=",TVAF),
          TDP = paste0("TDP=",TDP),
          CVAF = paste0("CVAF=",CVAF),
          CDP = paste0("CDP=",CDP),
          PURPLE_AF = paste0("PURPLE_AF=",PURPLE_AF),
          PURPLE_CN = paste0("PURPLE_CN=",PURPLE_CN),
          PURPLE_CALL_TIER = paste0("PURPLE_CALL_TIER=",PURPLE_CALL_TIER)
        ) |>
        dplyr::mutate(INFO = paste(
          TVAF, TDP, CVAF, CDP, PURPLE_AF, 
          PURPLE_CN, PURPLE_CALL_TIER, sep=";")) |>
        dplyr::mutate(QUAL = '.', FILTER = 'PASS', ID = '.') |>
        dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) |>
        dplyr::distinct()
      
      options(scipen=999)
      
      if(NROW(tumor_data) <= 10000 & NROW(tumor_data) > 5000){
	 vep_buffer_size <- 1000
      }
      if(NROW(tumor_data) <= 20000 & NROW(tumor_data) > 10000){
         vep_buffer_size <- 2000
      }
      if(NROW(tumor_data) <= 30000 & NROW(tumor_data) > 20000){
         vep_buffer_size <- 4000
      }
      if(NROW(tumor_data) <= 100000 & NROW(tumor_data) > 30000){
         vep_buffer_size <- 5000
      }
      if(NROW(tumor_data) > 100000){
         vep_buffer_size <- 10000
      }
      system(paste0('cat hartwig_vcf_header.txt > ',vcf_fname))
      write.table(tumor_data, 
                  vcf_fname, 
                  quote = F, 
                  row.names = F, 
                  col.names = F, 
                  append = T, 
                  sep = "\t")
      system(paste0('bgzip ', vcf_fname))
      system(paste0('tabix -p vcf ', vcf_fname,'.gz'))
    }
  }
  
  apptainer_command <- paste0(
    "apptainer exec --writable-tmpfs --no-home ",
    "--env \"XDG_CACHE_HOME=/data/sigven/quarto_cache_home\" ",
    "-B /data/sigven/hartwig_pcgr_annotation/pcgr_input:/mnt/input ",
    "-B /data/sigven/hartwig_pcgr_annotation/pcgr_output:/mnt/output ",
    "-B /data/sigven/pcgr:/mnt/bundle ",
    "-B /data/sigven/data__vep/.vep:/mnt/.vep ",
    "/data/sigven/hartwig_pcgr_annotation/pcgr_2.0.3.singularity.sif ",
    "pcgr --input_vcf /mnt/input/",tumor_sample_id,".vcf.gz ",
    "--input_cna /mnt/input/",tumor_sample_id,".cna.tsv ",
    "--refdata_dir /mnt/bundle ",
    "--vep_dir /mnt/.vep ",
    "--output_dir /mnt/output/",
    tsite_folder," ",
    "--genome_assembly grch37 ",
    "--tumor_site ",tsite_code," ",
    "--sample_id ",tumor_sample_id," ",
    "--retained_info_tag PURPLE_CALL_TIER ",
    "--tumor_af_tag TVAF --tumor_dp_tag TDP ",
    "--control_af_tag CVAF --control_dp_tag CDP ",
    "--estimate_msi --estimate_tmb --estimate_signatures ",
    "--vep_buffer_size ",vep_buffer_size," ",
    "--assay WGS ",
    "--force_overwrite ",
    "--vcf2maf ")
  
  if(!is.na(purity) & !is.na(ploidy)){
    apptainer_command <- paste0(
      apptainer_command,
      "--tumor_ploidy ",ploidy," ",
      "--tumor_purity ",purity," ")
  }
  
  apptainer_command <- paste0(
    apptainer_command,
    "> /data/sigven/hartwig_pcgr_annotation/pcgr_output/",tsite_folder,"/",tumor_sample_id,".pcgr.log 2>&1")
  
  write(apptainer_command, 
              file = "pcgr_commands.txt", sep = "\n", append = TRUE)
  i <- i + 1
}
  
  
#split --additional-suffix .sh -d -l 200 pcgr_commands.txt run_hartwig_
