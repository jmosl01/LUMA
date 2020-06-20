#' @title Culls the Peak.list by Void volume
#'
#' @export
#' @description Culls the Peak.list by removing all components coming from the
#'   void volume. See \code{remove_void_volume} for more details. For working
#'   examples, see \code{InitWorkflow, AnnotatePeaklist, CombineFeatures,
#'   CombinePeaklists, CullBackground, CullMF, CullCV, CullVoidVolume,
#'   FormatForMetaboAnalystR, FormatForSIMCA, NormalizePeaklists, ParseCAMERA,
#'   SimplyPeaklists, FinalWorkflow}.
#' @param from.table from which table should LUMA pull the Peak.list
#' @param to.table to which should LUMA save the modified Peak.list
#' @param method which method to apply to trim by retention time.  See
#'   \code{trim_rt} for details.
#' @return NULL
#' @examples
#' \dontrun{
#' library(LUMA)
#' db.dir <- system.file('extdata', package = "LUMA")
#' InitWorkflow(db.dir = db.dir)
#' CullVoidVolume(from.table = "From_CAMERA", to.table = "Trimmed_by_RT")
#' }
CullVoidVolume <- function(from.table,to.table,method) {
  ##Culls Peaklist by RT > void volume
  Peak.list <- read_tbl(mytable = from.table,
                        peak.db = peak_db,
                        asdf = TRUE)

  attributes(Peak.list)
  Peak.list.trimmed <- remove_void_volume(Peak.list = Peak.list,
                                          search.par = data.frame(ppm = ppm.cutoff,
                                                                  rt = rt.cutoff,
                                                                  Voidrt = Voidrt,
                                                                  Corr.stat.pos = Corr.stat.pos,
                                                                  Corr.stat.neg = Corr.stat.neg,
                                                                  CV = cv.cutoff,
                                                                  Minfrac = mf.cutoff,
                                                                  Endogenous = Endogenous.thresh,
                                                                  Solvent = Solvent.ratio,
                                                                  gen.plots = gen.plots,
                                                                  keep.singletons = keep.singletons),
                                          method = method)
  write_tbl(mydf = Peak.list.trimmed,
            peak.db = peak_db,
            myname = to.table)
}

#' @title Culls the Peak.list by Percent CV
#'
#' @export
#' @description Culls the Peak.list by removing all components with a
#'   coefficient of variation (CV) greater than the user-specified cutoff across
#'   all Pooled QC samples. See \code{trim_cv} for more details. For working
#'   examples, see \code{InitWorkflow, AnnotatePeaklist, CombineFeatures,
#'   CombinePeaklists, CullBackground, CullMF, CullCV, CullVoidVolume,
#'   FormatForMetaboAnalystR, FormatForSIMCA, NormalizePeaklists, ParseCAMERA,
#'   SimplyPeaklists, FinalWorkflow}.
#' @param from.table from which table should LUMA pull the Peak.list
#' @param to.table to which should LUMA save the modified Peak.list
#' @param QC.id character vector specifying identifier in filename designating a
#'   Pooled QC sample.  Only the first value will be used.  Default is
#'   \code{"Pooled_QC_"}
#' @return NULL
#' @examples
#' \dontrun{
#' library(LUMA)
#' db.dir <- system.file('extdata', package = "LUMA")
#' InitWorkflow(db.dir = db.dir)
#' CullCV(from.table = "From_CAMERA", to.table = "Trimmed_by_CV")
#' }
CullCV <- function(from.table,to.table,QC.id) {

  #Set default values
  if(missing(QC.id))
    QC.id <- "Pooled_QC_"

  #Culls Peaklist by CV
  Peak.list <- read_tbl(mytable = from.table,
                        peak.db = peak_db,
                        asdf = TRUE)

  Peak.list.trimmed <- trim_cv(Peak.list = Peak.list,
                               QC.id = QC.id,
                               search.par = data.frame(ppm = ppm.cutoff,
                                                       rt = rt.cutoff,
                                                       Voidrt = Voidrt,
                                                       Corr.stat.pos = Corr.stat.pos,
                                                       Corr.stat.neg = Corr.stat.neg,
                                                       CV = cv.cutoff,
                                                       Minfrac = mf.cutoff,
                                                       Endogenous = Endogenous.thresh,
                                                       Solvent = Solvent.ratio,
                                                       gen.plots = gen.plots,
                                                       keep.singletons = keep.singletons))


  if(!all(is.na(Peak.list.trimmed))) { #Only writes table to database if Peak.list.trimmed contains actual values

    write_tbl(mydf = Peak.list.trimmed,
              peak.db = peak_db,
              myname = to.table)

  }
}

#' @title Culls the Peak.list by Minimum Fraction
#'
#' @export
#' @description Culls the Peak.list by removing all components with minimum
#'   fraction less than the user-specified cutoff. See \code{trim_minfrac} for
#'   more details. For working examples, see \code{InitWorkflow,
#'   AnnotatePeaklist, CombineFeatures, CombinePeaklists, CullBackground,
#'   CullMF, CullCV, CullVoidVolume, FormatForMetaboAnalystR, FormatForSIMCA,
#'   NormalizePeaklists, ParseCAMERA, SimplyPeaklists, FinalWorkflow}.
#' @param from.table from which table should LUMA pull the Peak.list
#' @param to.table to which should LUMA save the modified Peak.list
#' @param method which method to apply to trim by minimum fraction values.  See
#'   \code{trim_minfrac} for details.
#' @return NULL
#' @examples
#' \dontrun{
#' library(LUMA)
#' db.dir <- system.file('extdata', package = "LUMA")
#' InitWorkflow(db.dir = db.dir)
#' CullMF(from.table = "From_CAMERA_with_MinFrac", to.table = "Trimmed_by_MinFrac")
#' }
CullMF <- function(from.table,to.table,method) {
  #Trims Peaklist by MinFrac
  Peak.list <- read_tbl(mytable = from.table,
                        peak.db = peak_db,
                        asdf = TRUE)

  class(method) <- method

  Peak.list.trimmed <- trim_minfrac(Peak.list = Peak.list,
                                    search.par = data.frame(ppm = ppm.cutoff,
                                                            rt = rt.cutoff,
                                                            Voidrt = Voidrt,
                                                            Corr.stat.pos = Corr.stat.pos,
                                                            Corr.stat.neg = Corr.stat.neg,
                                                            CV = cv.cutoff,
                                                            Minfrac = mf.cutoff,
                                                            Endogenous = Endogenous.thresh,
                                                            Solvent = Solvent.ratio,
                                                            gen.plots = gen.plots,
                                                            keep.singletons = keep.singletons),
                                    object = method)


  if(!is.null(Peak.list.trimmed)) {

    write_tbl(mydf = Peak.list.trimmed,
              peak.db = peak_db,
              myname = to.table)

  }
}

#' @title Culls the Peak.list by Background
#'
#' @export
#' @description Culls the Peak.list by removing background components. See
#'   \code{remove_background_peaks} for more details. For working examples, see
#'   \code{InitWorkflow, AnnotatePeaklist, CombineFeatures, CombinePeaklists,
#'   CullBackground, CullMF, CullCV, CullVoidVolume, FormatForMetaboAnalystR,
#'   FormatForSIMCA, NormalizePeaklists, ParseCAMERA, SimplyPeaklists,
#'   FinalWorkflow}.
#' @param from.tables character vector of table names to draw from databases.
#'   First value should be table name from peak database, second should be table
#'   name from solvent database.
#' @param to.tables character vector of table names to write to new database.
#'   Should be twice as long as the number of processed ion modes.
#' @param method which method to use to remove background components.
#' @param db.list vector of database names containing results from processing
#'   modules. Can be left blank. See connect_lumadb for details.
#' @param db.dir directory containing the databases. Default is 'db'
#' @param new.db what should the new database be called. Default is
#'   'Peaklist_db'
#' @param lib.db what should the library containing the background components be
#'   called. Default is 'Background Components Library'
#' @return NULL
#' @examples
#' \dontrun{
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)){
#'
#'   db.dir <- system.file('extdata', package = "lcmsfishdata")
#'   InitWorkflow(db.dir = db.dir)
#'   CullBackground(from.tables =
#'   c("Trimmed_by_MinFrac","Combined_Isotopes_and_Adducts"),
#'                to.tables = c("Peaklist_Pos_Solvent_Peaks_Removed",
#'                              "Peaklist_Neg_Solvent_Peaks_Removed",
#'                              "Peaklist_Pos_Solvent Peaks Only",
#'                              "Peaklist_Neg_Solvent Peaks Only"),
#'                method = "monoMass")
#'
#'   }
#' }
CullBackground <- function(from.tables,to.tables,method,db.list,db.dir,new.db,lib.db) {

  #Initialize all other global variables
  peak_db <- NULL

  #Set default values
  if(missing(db.list))
    db.list = c("Peaklist_Pos","Peaklist_Neg","Blanks_Pos","Blanks_Neg")
  if(missing(db.dir))
    db.dir = "db"
  if(missing(new.db))
    new.db = "Peaklist_db"
  if(missing(lib.db))
    lib.db = "Background Components Library"

  Peak.list.trimmed <- remove_background_peaks(Peak.list = NULL,
                                       Sample.df = data.frame(Sex = Sexes,
                                                              Class = Classes,
                                                              n = no.Samples,
                                                              Endogenous = Endogenous),
                                       search.par = data.frame(ppm = ppm.cutoff,
                                                               rt = rt.cutoff,
                                                               Voidrt = Voidrt,
                                                               Corr.stat.pos = Corr.stat.pos,
                                                               Corr.stat.neg = Corr.stat.neg,
                                                               CV = cv.cutoff,
                                                               Minfrac = mf.cutoff,
                                                               Endogenous = Endogenous.thresh,
                                                               Solvent = Solvent.ratio,
                                                               gen.plots = gen.plots,
                                                               keep.singletons = keep.singletons),
                                       tbl.id = from.tables,
                                       ion.mode = ion.mode,
                                       method = method,
                                       db.list = db.list,
                                       lib.db = lib.db,
                                       db.dir = db.dir,
                                       new.db = new.db)

  #Update Peaklist database connection globally
  peak_db <<- peak_db <- connect_peakdb(file.base = new.db,db.dir = db.dir)

  write_tbl(mydf = Peak.list.trimmed$Positive$Peaklist,
            peak.db = peak_db,
            myname = to.tables[1])

  write_tbl(mydf = Peak.list.trimmed$Negative$Peaklist,
            peak.db = peak_db,
            myname = to.tables[2])

  write_tbl(mydf = Peak.list.trimmed$Positive$Solvent_Peaks,
            peak.db = peak_db,
            myname = to.tables[3])

  write_tbl(mydf = Peak.list.trimmed$Negative$Solvent_Peaks,
            peak.db = peak_db,
            myname = to.tables[4])
}
