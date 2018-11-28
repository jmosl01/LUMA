#' @title Normalize Peaklists from different ion modes
#'
#' @export
#' @description Normalizes Peaklists from positive and negative ionization modes separately.
#' See replace_zeros_normalize for more details.
#' @param from.table from which table should LUMA read the Peaklist
#' @param to.table to which table should LUMA save the modified Peaklist
#' @param peak.db what database contains the Peaklists to be combined.
#' Default is 'Peaklist_db'
#' @param db.dir directory containing the database.
#' Default is 'db'
NormalizePeaklists <- function(from.table,to.table,peak.db,db.dir) {

  #Set default values
  if(missing(peak.db))
    peak.db = "Peaklist_db"
  if(missing(db.dir))
    db.dir = "db"

  cat("Replacing Zeros and Normalizing Peaklists.")

  #Update Peaklist database connection globally
  peak_db <<- peak_db <- connect_peakdb(file.base = peak.db,db.dir = db.dir)

  Norm.Peaklist <- replace_zeros_normalize(Peak.list = NULL,
                                           Sample.df = data.frame(Sex = Sexes,
                                                                  Class = Classes,
                                                                  n = no.Samples,
                                                                  Endogenous = Endogenous),
                                           tbl.id = from.table,
                                           peak.db = peak_db,
                                           asdf = TRUE)

  write_tbl(Norm.Peaklist,
            peak.db = peak_db,
            myname = to.table)


}
