library(DBI)
library(RSQLite)

## Goal ##
# The goal of this script is to summarize the contents of the results db

bDir <- "/Users/larsonhogstrom/Documents/oncology_biomarkers/resultsDb"
mydb <- DBI::dbConnect(RSQLite::SQLite(), paste0(bDir,"/actionable-biomarker-db.sqlite"))

# show metadata table
subsetHmeta <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM hartwigMetadata LIMIT 5')
print(head(subsetHmeta))

# view a head and dim of each table
tableList <- RSQLite::dbListTables(mydb)
for (tbl in tableList) {
  print(tbl)
  sQuery <- paste0('SELECT * FROM ',tbl)
  tblDf <- RSQLite::dbGetQuery(mydb, sQuery)
  print(dim(tblDf))
  nColMax <- min(15,dim(tblDf)[[2]])
  print(head(tblDf[1:5,colnames(tblDf)[1:nColMax]]))
}

RSQLite::dbDisconnect(mydb)
