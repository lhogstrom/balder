
#' Combine patient observed snv/indels with sample info, oncotree codes, and 
#' OncoKB results
#'
#' @param poVariants patient observed snv/indel table
#' @param sampleInfoCompiled compiled sample information with cancer type
#' @param ot_code_full table of oncotree codes
#' @param oncokbRes oncokb results matching patient observed variants
#'
#' @return compiled variant and sample information
combine_patient_table_and_oncokb <-
  function(poVariants,
           sampleInfoCompiled,
           ot_code_full,
           oncokbRes) {
    
    # oncotree code level names
    otCodeCols <- c("ot_code", "ot_name", "Highest_Non_Null_Level", "oncotree_level",
                    "level_1_disease", "ot_code_level_1",
                    "level_2_disease", "ot_code_level_2",
                    "level_3_disease", "ot_code_level_3",
                    "level_4_disease", "ot_code_level_4",
                    "level_5_disease", "ot_code_level_5",
                    "level_6_disease", "ot_code_level_6")
    
    ### join sample info to variant data
    svCompiledPrelim <- poVariants %>%
      dplyr::inner_join(sampleInfoCompiled[,!colnames(sampleInfoCompiled)=="SourceStudy"],
                        by=c("Tumor_Sample_Barcode"="SAMPLE_ID")) %>% # exclude samples that do not have any sample information
      dplyr::left_join(ot_code_full[,otCodeCols],by=c("ONCOTREE_CODE"="ot_code")) %>% # join ot_codes to annotation
      dplyr::filter(!SAMPLE_TYPE == "Cell line")
    print(dim(svCompiledPrelim))
    
    # confirm oncokb entries match order of annotated variants
    # table(oncokbRes$ONCOKB_Tumor_Sample_Barcode == svCompiledPrelim$Tumor_Sample_Barcode,exclude=NULL)
    # table(oncokbRes$ONCOKB_HGVSp == svCompiledPrelim$HGVSp,exclude=NULL)
    # table(oncokbRes$ONCOKB_ONCOTREE_CODE == svCompiledPrelim$ONCOTREE_CODE,exclude=NULL)
    
    # assume oncokb row order matches patient observed variants
    svCompiled <- cbind(svCompiledPrelim,oncokbRes[,6:ncol(oncokbRes)]) # exclude columns in oncokb output that are already in data product
    
    return(svCompiled)
    
  }

#' Function that joins snv/indel data with clinical annotations according 
#' to AA change and genomic coordinates 
#'
#' @param svCompiled Compiled snv/indel data with patient's cancer type and oncotree code
#' @param dbRules List of compiled clinical annotation rules
#'
#' @return exhaustive list of annotated variants based on AA coordinates and genomic matching
#' @export
#'
#' @examples
annotation_matching_genomic_coordinate_and_AA_change <-
  function(svCompiled,
           dbRules) {
    aa.match.exhaustive <- svCompiled %>%
      dplyr::inner_join(dbRules,by=c("Hugo_Symbol"="gene",
                                     "AAChangeObserved"="AAChange"),
                        relationship = "many-to-many") %>%
      dplyr::mutate(matchAndVarID=paste0(balderVariantID,"-",balderRuleID)) %>%
      dplyr::mutate(AAChange=AAChangeObserved,
                    gene=Hugo_Symbol) # create dummy variable for combining
    
    genome.match.exhaustive <- svCompiled %>%
      dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome_annotation",
                                     "Start_Position"="pos",
                                     "Reference_Allele"="ref",
                                     "Tumor_Seq_Allele2"="alt"),
                        relationship = "many-to-many") %>%
      dplyr::mutate(matchAndVarID=paste0(balderVariantID,"-",balderRuleID)) %>%
      dplyr::mutate(chromosome_annotation=Chromosome,
                    pos=Start_Position,
                    ref=Reference_Allele,
                    alt=Tumor_Seq_Allele2) # create dummy variable for combining
    
    ### AA mismatch events between observed and annotation names? 
    # how often is there a mismatch b/t AA change calls for genomic-based targets? 
    #table(genome.match.exhaustive$AAChangeObserved == genome.match.exhaustive$AAChange)
    
    aa.genome.exahustive <- rbind(aa.match.exhaustive,genome.match.exhaustive) %>% 
      dplyr::group_by(balderVariantID,balderRuleID) %>%
      dplyr::filter(dplyr::row_number()==1) %>%
      dplyr::mutate(cancerTypeMatchPrimary=oncotree_code_annotation == ONCOTREE_CODE & !is.na(ONCOTREE_CODE), # ot code from variants (CAPS) & ot code from annotation (not caps)
                    cancerTypeMatchLevel1=oncotree_code_annotation == ot_code_level_1 & !is.na(ot_code_level_1),
                    cancerTypeMatchLevel2=oncotree_code_annotation == ot_code_level_2 & !is.na(ot_code_level_2),
                    cancerTypeMatchLevel3=oncotree_code_annotation == ot_code_level_3 & !is.na(ot_code_level_3),
                    cancerTypeMatchLevel4=oncotree_code_annotation == ot_code_level_4 & !is.na(ot_code_level_4),
                    cancerTypeMatchLevel5=oncotree_code_annotation == ot_code_level_5 & !is.na(ot_code_level_5),
                    cancerTypeMatchLevel6=oncotree_code_annotation == ot_code_level_6 & !is.na(ot_code_level_6),
                    cancerTypeMatch = cancerTypeMatchPrimary | cancerTypeMatchLevel1 | cancerTypeMatchLevel2 | cancerTypeMatchLevel3 | cancerTypeMatchLevel4 | cancerTypeMatchLevel5 | cancerTypeMatchLevel6)
    
    aa.genome.exahustive$annotationMatchGenomicCoord <- aa.genome.exahustive$matchAndVarID %in% genome.match.exhaustive$matchAndVarID
    aa.genome.exahustive$annotationMatchAAChange <- aa.genome.exahustive$matchAndVarID %in% aa.match.exhaustive$matchAndVarID
    
    return(aa.genome.exahustive)
    
  }
