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

#' @title Formatting of metabolite data for MetaboAnalystR.
#'
#' @export
#' @description This function initializes objects that will hold the metabolite data, formats peak intensity data into one of the formats acceptable by MetaboAnalystR, and sets the metabolite data object.
#' @param Peak.list data frame containing combined ion mode peaklist with ion mode duplicates removed.
#' @param Sample.df data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Sample.data ata frame with phenotype data as columns and a row for each study sample.  First column must be a unique sample identifier with the header 'CT-ID'.  Phenotype columns may vary, but must include two columns called 'Plate Number' and 'Plate Position' for determining run order.
#' @param data.type NULL
#' @param anal.type NULL
#' @param paired NULL
#' @return mSetObj
FormatForMetaboAnalystR <- function(Peak.list, Sample.df, Sample.data, data.type = "pktable", anal.type = "stat", paired = FALSE)
{
  ##-----------------------------------------------------------------------------------------
  ## Initialize the metabolite object.
  ##-----------------------------------------------------------------------------------------
  dataSet <- list();
  dataSet$type <- data.type;
  dataSet$design.type <- "regular";    # one factor to two factor
  dataSet$cls.type <- "disc";
  dataSet$format <- "rowu";
  dataSet$paired <- paired;
  analSet <- list();
  analSet$type <- anal.type;

  mSetObj <- list();
  mSetObj$dataSet <- dataSet;
  mSetObj$analSet <- analSet;
  mSetObj$imgSet <- list();
  mSetObj$msgSet <- list();                                 # store various message during data processing
  mSetObj$msgSet$msg.vec <- vector(mode = "character");     # store error messages
  mSetObj$cmdSet <- vector(mode = "character");             # store R command

  # Set the class type for the mSetObj.
  class(mSetObj) <- data.type

  mSetObj <- format_MetabolomicData(mSetObj, Peak.list, Sample.df, Sample.data)

  return(mSetObj)
}
