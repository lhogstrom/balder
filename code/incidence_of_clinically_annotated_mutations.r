library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)

### connect to result DB and get variants
bDir <- "../../data/processed/balderResultsDb"
dbName <- paste0(bDir,"/balder-harmonized-biomarker-data-v20240311.sqlite")
harmonizedDb <- DBI::dbConnect(RSQLite::SQLite(), dbName)

# load clinical annotations
dbRules <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM MoaCiVICRuleEntries')
dbAlterations <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM actionableSNVsByAAChange')
dbGenome <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM actionableSNVsByGenomicCoordinate')

# load patient observed variants
poVariants<- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM patientObservedVariantTable') %>%
  dplyr::mutate(AAChange=stringr::str_replace(HGVSp_Short,"p.", ""))
print(dim(poVariants))

# load sample info for observed variants
sampleInfoCompiled <- RSQLite::dbGetQuery(harmonizedDb, 'SELECT * FROM harmonizedSampleInfo')
print(dim(sampleInfoCompiled))
RSQLite::dbDisconnect(mydb)

### join sample info to variant data
svCompiled <- poVariants %>%
  dplyr::left_join(sampleInfoCompiled[,!colnames(sampleInfoCompiled)=="SourceStudy"],
                   by=c("Tumor_Sample_Barcode"="SAMPLE_ID"))


# Computing actionable variant counts for different match approaches
# 
# Different actionablity matching approaches were computed:
#   
# Unrestricted where no cancer type matching was imposed
# 
# Restricted with raw cancer type string matching

###
### join with observed variant table with clinically annotated - by genome position
genome.matched.unrestricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome",
                                 "Start_Position"="pos",
                                 "Reference_Allele"="ref",
                                 "Tumor_Seq_Allele2"="alt")) %>%
  dplyr::group_by(Chromosome,Start_Position,Tumor_Seq_Allele2) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(genome.matched.unrestricted))
print(table(genome.matched.unrestricted$source,exclude = NULL))

### Unristricted match events - AA Change
aa.matched.unrestricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Hugo_Symbol"="gene",
                                      "AAChange"="AAChange")) %>%
  dplyr::group_by(Hugo_Symbol,AAChange) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(aa.matched.unrestricted))
print(table(aa.matched.unrestricted$source,exclude = NULL))

### impose cancer type matches (genome) - oncotree
genome.matched.ot.restricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Chromosome"="chromosome",
                                 "Start_Position"="pos",
                                 "Reference_Allele"="ref",
                                 "Tumor_Seq_Allele2"="alt",
                                 "ONCOTREE_CODE"="oncotree_code")) %>%
  dplyr::group_by(Chromosome,Start_Position,Tumor_Seq_Allele2,ONCOTREE_CODE) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   disease=paste0(unique(Disease),collapse=";"),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(genome.matched.unrestricted))
print(table(genome.matched.unrestricted$source,exclude = NULL))

### impose cancer type matches (protein) - oncotree
aa.matched.ot.restricted <- svCompiled %>%
  dplyr::inner_join(dbRules,by=c("Hugo_Symbol"="gene",
                                 "AAChange"="AAChange",
                                 "ONCOTREE_CODE"="oncotree_code")) %>%
  dplyr::group_by(Hugo_Symbol,AAChange,ONCOTREE_CODE) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(Tumor_Sample_Barcode),
                   disease=paste0(unique(Disease),collapse=";"),
                   source=paste0(unique(source),collapse=";"),
                   Drugs=paste0(unique(Drugs),collapse=";"),
                   ReferenceOrTrialID=paste0(unique(ReferenceOrTrialID),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(aa.matched.unrestricted))
print(table(aa.matched.unrestricted$source,exclude = NULL))

### Reporting table

n.patients <- length(unique(svCompiled$Tumor_Sample_Barcode))
n.variants <- length(unique(dbAlterations$AAChange))

c1 <- c("Restricted","Restricted","Unestricted","Unestricted")
c2 <- c("genome","AA change", "genome", "AA change")
c3 <- c(dim(genome.matched.ot.restricted)[[1]], 
        dim(aa.matched.ot.restricted)[[1]], 
        dim(genome.matched.unrestricted)[[1]], 
        dim(aa.matched.unrestricted)[[1]])
#c4 <- paste0(round(100*(c3/n.variants),1),"%")
c5 <- c(sum(genome.matched.ot.restricted$n.patients.mutated), 
        sum(aa.matched.ot.restricted$n.patients.mutated), 
        sum(genome.matched.unrestricted$n.patients.mutated), 
        sum(aa.matched.unrestricted$n.patients.mutated))
c6 <- paste0(round(100*(c5/n.patients),1),"%") 
# Create a data frame
df <- data.frame(c1=c1, c2=c2, c3=c3, c5=c5, c6=c6)
colnames(df) <- c("Cancer type match restriction", "Cancer type matching approach", "number of variants","number of patients","percent of patients")
# Print the data frame
#print(df)
knitr::kable(df, )

