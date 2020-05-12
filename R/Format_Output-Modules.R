#' @title Format Peaklists for input into SIMCA
#'
#' @export
#' @description Formats Peaklists processed by LUMA and exports a csv file for
#'   direct import into the multivariate statistical software SIMCA. See
#'   format_simca for more details.
#' @param from.table from which table should LUMA read the Peaklist
#' @param to.csv to what filename (excluding .csv extension) should LUMA save
#'   the formatted Peaklist
#' @param peak.db what database contains the Peaklists to be combined. Default
#'   is 'Peaklist_db'
#' @param db.dir directory containing the database. Default is 'db'
#' @examples
#' \dontrun{
#' library(LUMA)
#' db.dir <- system.file('extdata/', package = "LUMA")
#' InitWorkflow(db.dir = db.dir)
#' AnnotatePeaklist(from.table = "From CAMERA", to.table = "Annotated")
#' FormatForSIMCA(from.table = "Annotated", to.csv = "Peaklist_for_SIMCA",
#' peak.db = peak_db, db.dir = db.dir)
#' }
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
#' @description This function initializes objects that will hold the metabolite
#'   data, formats peak intensity data into one of the formats acceptable by
#'   MetaboAnalystR, and sets the metabolite data object.
#' @param from.table from which table should LUMA read the Peaklist
#' @param to.csv to what filename (excluding .csv extension) should LUMA save
#'   the formatted Peaklist
#' @param data.type What type of data will be generated. See usage and
#'   format_MetabolomicData for options.
#' @param anal.type character Indicates the analysis module the data will be
#'   used for. See usage and documentation for MetaboAnalystR::InitDataObjects()
#'   for options.
#' @param paired logical Indicate if the data is paired or not. Default is FALSE
#' @param peak.db what database contains the Peaklists to be combined. Default
#'   is 'Peaklist_db'
#' @param db.dir directory containing the database. Default is 'db'
#' @return mSetObj
#' @examples
#' \dontrun{
#' library(LUMA)
#' db.dir <- system.file('extdata/', package = "LUMA")
#' InitWorkflow(db.dir = db.dir)
#' AnnotatePeaklist(from.table = "From CAMERA", to.table = "Annotated")
#' FormatForMetaboAnalystR(from.table = "Annotated", to.csv = "Peaklist_for_MetaboAnalyist", peak.db = peak_db, db.dir = db.dir)
#' }
FormatForMetaboAnalystR <- function(from.table, to.csv, data.type = "pktable",
                                    anal.type = "stat", paired = FALSE, peak.db, db.dir)
{

  #Set default values
  if(missing(paired))
    paired = FALSE
  if(missing(peak.db))
    peak.db = "Peaklist_db"
  if(missing(db.dir))
    db.dir = "db"

  cat("Formatting for MetaboAnalystR.")

  #Update Peaklist database connection globally
  peak_db <<- peak_db <- connect_peakdb(file.base = peak.db,db.dir = db.dir)



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

  mSetObj <- format_MetabolomicData(mSetObj, Peak.list = NULL,

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


  conc <- mSetObj$dataSet$orig
  MT.data <- cbind(Sample = rownames(conc),Class = mSetObj$dataSet$orig.cls,conc)

  write.table(MT.data, file = paste(to.csv,".csv",sep = ""), sep = ",", row.names = FALSE)

  #Generate Metadata files for mSetObj
  peak.list <- read_tbl(from.table, peak_db)

  output_MetaData(mSetObj, peak.list,
                  Sample.df = data.frame(Sex = Sexes,
                                         Class = Classes,
                                         n = no.Samples,
                                         Endogenous = Endogenous),

                  Sample.data = cbind.data.frame(CT.ID,
                                                 Plate.Number,
                                                 Plate.Position,
                                                 Sample.phenodata))


  return(mSetObj)
}
