library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)

### connect to result DB and get variants
bDir <- "../../data/processed/balderResultsDb"
mydb <- DBI::dbConnect(RSQLite::SQLite(), paste0(bDir,"/actionable-biomarker-db.sqlite"))

# load clinical annotations
dbRules <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM MoaCiVICRuleEntries')
dbAlterations <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM actionableSNVsByAAChange')
dbGenome <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM actionableSNVsByGenomicCoordinate')

# load patient observed variants
poVariants<- RSQLite::dbGetQuery(mydb, 'SELECT * FROM patientObservedVariantTable')
print(dim(svCompiled))

# load sample info for observed variants
sampleInfoCompiled <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM harmonizedSampleInfo')
print(dim(svCompiled))
RSQLite::dbDisconnect(mydb)

# columns needed: Hugo_Symbol, AA change

### join sample info to variant data
svCompiled <- poVariants %>%
  dplyr::left_join(sampleInfoCompiled[,!colnames(sampleInfoCompiled)=="SourceStudy"],
                   by=c("sample"="SAMPLE_ID"))

###
### join with observed variant table with clinically annotated - by genome position
genome.matched <- svCompiled %>%
  dplyr::left_join(dbGenome,by=c("chrom"="chr",
                            "pos"="pos",
                            "ref"="ref",
                            "alt"="alt"))
table(is.na(genome.matched$source),exclude=NULL)
# 
# ###
# ### join with MOA - by protein change
# dbAlteration$MOAset <- "Y"
# msk.protein <- msk %>%
#   dplyr::left_join(dbAlteration,by=c("Hugo_Symbol"="gene",
#                                      "protein_change"="AAChange"))
# print(table(msk.protein$MOAset,msk.protein$source,exclude = NULL))
# print(table(msk.protein$feature_type,exclude = NULL))
# 
# msk.protein.cnt <- msk %>%
#   dplyr::inner_join(dbAlteration,by=c("Hugo_Symbol"="gene",
#                                       "protein_change"="AAChange")) %>%
#   dplyr::group_by(Hugo_Symbol,protein_change) %>%
#   dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID),
#                    source=paste0(unique(source),collapse=";")) %>%
#   dplyr::arrange(desc(n.patients.mutated))
# print(dim(msk.protein.cnt))
# print(table(msk.protein.cnt$source,exclude = NULL))