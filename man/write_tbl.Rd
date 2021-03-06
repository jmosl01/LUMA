% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database-functions.R
\name{write_tbl}
\alias{write_tbl}
\title{Writes Peak.list to database}
\usage{
write_tbl(mydf, peak.db, myname)
}
\arguments{
\item{mydf}{tbl or dataframe to write}

\item{peak.db}{Formal class SQLiteConnection}

\item{myname}{character what should the table be called}
}
\value{
a tbl object in the remote source
}
\description{
Writes Peak.list to an RSQLite database. Peak.list can be a tibble or data frame.
}
\examples{
library(LUMA)
if(require(RSQLite, quietly = TRUE)) {
file <- system.file("extdata","Sample_Data.csv", package =  "LUMA") # is case sensitive on Linux
sample_data <- read.table(file, header = TRUE, sep = ",")
mzdatafiles <- sample_data$CT.ID
file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode =
"Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
peak_db <- connect_peakdb(file.base = file.base, mem = TRUE)
dbIsValid(peak_db) #Database is valid
dbListTables(peak_db) #But no tables yet
mydf <- Peaklist_Pos$From_CAMERA
write_tbl(mydf = mydf, myname = "From_CAMERA", peak.db = peak_db)
dbListTables(peak_db) #Now has tables
dbDisconnect(peak_db)
}
}
