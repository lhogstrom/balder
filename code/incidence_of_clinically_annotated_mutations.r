library(dplyr)
library(tidyr)
library(DBI)
library(RSQLite)

### connect to result DB and get variants
bDir <- "../../data/processed/balderResultsDb"
mydb <- DBI::dbConnect(RSQLite::SQLite(), paste0(bDir,"/actionable-biomarker-db.sqlite"))
svCompiled <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM patientObservedVariantTable')
print(dim(svCompiled))

sampleInfoCompiled <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM harmonizedSampleInfo')
print(dim(svCompiled))
RSQLite::dbDisconnect(mydb)

# columns needed: Hugo_Symbol, AA change, cancer type

### join sample info to variant data
svCompiled <- svCompiled %>%
  dplyr::left_join(sampleInfoCompiled,by=c("sample"="SAMPLE_ID",
                                     "protein_change"="AAChange"))




###
### join with MOA - by protein change
dbAlteration$MOAset <- "Y"
msk.protein <- msk %>%
  dplyr::left_join(dbAlteration,by=c("Hugo_Symbol"="gene",
                                     "protein_change"="AAChange"))
print(table(msk.protein$MOAset,msk.protein$source,exclude = NULL))
print(table(msk.protein$feature_type,exclude = NULL))

msk.protein.cnt <- msk %>%
  dplyr::inner_join(dbAlteration,by=c("Hugo_Symbol"="gene",
                                      "protein_change"="AAChange")) %>%
  dplyr::group_by(Hugo_Symbol,protein_change) %>%
  dplyr::summarize(n.patients.mutated=dplyr::n_distinct(PATIENT_ID),
                   source=paste0(unique(source),collapse=";")) %>%
  dplyr::arrange(desc(n.patients.mutated))
print(dim(msk.protein.cnt))
print(table(msk.protein.cnt$source,exclude = NULL))