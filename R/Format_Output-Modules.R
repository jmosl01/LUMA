#' @title Format Peaklists for input into SIMCA
#'
#' @export
#' @description Formats Peaklists processed by LUMA and exports a csv file for direct import into the multivariate statistical software SIMCA.
#' See format_simca for more details.
#' @param from.table from which table should LUMA read the Peaklist
#' @param to.csv to what filename (excluding .csv extension) should LUMA save the formatted Peaklist
#' @param peak.db what database contains the Peaklists to be combined.
#' Default is 'Peaklist_db'
#' @param db.dir directory containing the database.
#' Default is 'db'
FormatForSIMCA <- function(from.table,to.csv,peak.db,db.dir) {

  #Set default values
  if(missing(peak.db))
    peak.db = "Peaklist_db"
  if(missing(db.dir))
    db.dir = "db"

  cat("Formatting for SIMCA.")

  #Update Peaklist database connection globally
  peak_db <<- peak_db <- connect_peakdb(file.base = peak.db,db.dir = db.dir)

  SIMCA.data <- format_simca(Peak.list = NULL,
                             Sample.df = data.frame(Sex = Sexes,
                                                    Class = Classes,
                                                    n = no.Samples,
                                                    Endogenous = Endogenous),
                             Sample.data = cbind.data.frame(CT.ID,
                                                            Plate.Number,
                                                            Plate.Position,
                                                            Sample.phenodata),
                             tbl.id = from.table,
                             peak.db = peak_db)

  write.table(SIMCA.data, file = paste(to.csv,".csv",sep = ""), sep = ",", row.names = FALSE)


}
