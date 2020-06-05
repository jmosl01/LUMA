#' @title Combine Peaklists from different ion modes
#'
#' @export
#' @description Combines Peaklists from positive and negative ionization modes while removing ion mode duplicates.
#' See combine_ion_modes and search_IonDup for more details.
#' @param from.tables character vector of table names to draw from databases to be combined.
#' If gen.plots = TRUE, should be of length 2. Otherwise, only the first element will be used.
#' @param to.table to which table should LUMA save the modified Peak.list
#' @param method which method to use for removing ion mode duplicates.
#' See search_IonDup for available methods.
#' @param peak.db what database contains the Peaklists to be combined.
#' Default is 'Peaklist_db'
#' @param db.dir directory containing the database.
#' Default is 'db'
#' @param gen.plots logical indicating whether LUMA needs to generate plots to inspect ion mode duplicates.
#' Default is to check whether this variable was assigned by a call to InitWorkflow
#' @param CAMERA.pos which CAMERA object to use to plot metabolite groups for positive mode.
#' Default is to read from saved R objects by call to InitWorkflow
#' @param CAMERA.neg which CAMERA object to use to plot metabolite groups for negative mode.
#' Default is to read from saved R objects by call to InitWorkflow
#' @examples
#' \dontrun{
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'   db.dir <- system.file("extdata/", package = "lcmsfishdata")
#'   InitWorkflow(db.dir = db.dir)
#'   AnnotatePeaklist(from.table = "From_CAMERA", to.table = "Annotated")
#'   ParseCAMERA(from.table = "Annotated", to.table = "output_parsed", CAMERA.obj
#'   = anposGa)
#'   CombineFeatures(from.table = "output_parsed", to.table =
#'   "Combined_Isotopes_and_Adducts")
#'   CombinePeaklists(from.tables = c("Peaklist_Pos_Solvent Peaks Removed",
#'   "Peaklist_Neg_Solvent Peaks Removed"),
#'   to.table = "Peaklist_Combined",
#'   method = "monoMass", peak.db = peak_db, db.dir = db.dir)
#'   FinalWorkflow(peak_db = peak_db, lib_db = lib_db)
#'   }
#' }
CombinePeaklists <- function(from.tables,to.table,method,peak.db,db.dir,gen.plots,CAMERA.pos,CAMERA.neg) {

  #Set default values
  if(missing(peak.db))
    peak.db = "Peaklist_db"
  if(missing(db.dir))
    db.dir = "db"
  if(missing(CAMERA.pos))
    CAMERA.pos <- NULL
  if(missing(CAMERA.neg))
    CAMERA.neg <- NULL
  if(missing(gen.plots)) {
    if(exists("gen.plots",where = .GlobalEnv)) {
      mygen.plots <- get("gen.plots", envir = .GlobalEnv)
    } else {
      stop("You must set gen.plots manually. See InitWorkflow for more details!\n\n")
    }
  } else {
    mygen.plots = gen.plots
  }
  if(exists("BLANK",where = .GlobalEnv)) {
    #Update BLANK to FALSE globally, because it doesn't make sense to plot duplicates for blanks!
    BLANK <<- BLANK <- FALSE
  }


  #Update Peaklist database connection globally
  peak_db <<- peak_db <- connect_peakdb(file.base = peak.db,db.dir = db.dir)

  #Ion modes sanity check
  if(length(ion.mode) != 2) {
    stop(paste("ion.mode must be a vector of length 2. \nIt is currently of length ",length(ion.mode),sep = ""))
  }

  if(mygen.plots) {

    ##Check for existing positive mode CAMERA objects. If not specified, check for saved CAMERA objects.
    ##If none exist, return error.
    if(is.null(CAMERA.pos)) {
      PreProcesslist <- .set_PreProcessFileNames(ion.mode[1],BLANK)
      CAMERA.pos <- PreProcesslist[[2]]

      #CAMERA sanity check
      if(file.exists(CAMERA.pos)){
        cat("Reading in positive mode CAMERA files.\n\n")
        load(file = CAMERA.pos)
      } else {
        stop("No saved positive mode CAMERA files exist. \nBe sure to call InitWorkflow module for positive mode first!\n\n")
      }
    } else {
      if(!exists(CAMERA.pos)) {
        stop("Did not find the specified CAMERA object in the global environment!\n\n")
      }
    }
    CAMERA.pos <- .CAMERASanityCheck(CAMERA.pos)

    ##Check for existing negative mode CAMERA objects. If not specified, check for saved CAMERA objects.
    ##If none exist, return error.
    if(is.null(CAMERA.neg)) {
      PreProcesslist <- .set_PreProcessFileNames(ion.mode[2],BLANK)
      CAMERA.neg <- PreProcesslist[[2]]

      #CAMERA sanity check
      if(file.exists(CAMERA.neg)){
        cat("Reading in negative mode CAMERA files.\n\n")
        load(file = CAMERA.neg)
      } else {
        stop("No saved negative mode CAMERA files exist. \nBe sure to call InitWorkflow module for negative mode first!\n\n")
      }
    } else {
      if(!exists(CAMERA.neg)) {
        stop("Did not find the specified CAMERA object in the global environment!\n\n")
      }
    }
    CAMERA.neg <- .CAMERASanityCheck(CAMERA.neg)

    Peak.list.combined <- combine_ion_modes(Peak.list = NULL,
                                            search.par = data.frame(ppm = ppm.cutoff,
                                                                    rt = rt.cutoff,
                                                                    Voidrt = Voidrt,
                                                                    Corr.stat.pos = Corr.stat.pos,
                                                                    Corr.stat.neg = Corr.stat.neg,
                                                                    CV = cv.cutoff,
                                                                    Minfrac = mf.cutoff,
                                                                    Endogenous = Endogenous.thresh,
                                                                    Solvent = Solvent.ratio,
                                                                    gen.plots = mygen.plots,
                                                                    keep.singletons = keep.singletons),
                                            tbl.id = from.tables,
                                            method = method,
                                            peak.db = peak_db)

    write_tbl(Peak.list.combined,
              peak.db = peak_db,
              myname = to.table)

    Key.list <- plot_ionduplicate(anposGa = anposGa,
                                  annegGa = annegGa,
                                  mytable = to.table,
                                  peak.db = peak_db,
                                  gen.plots = mygen.plots,
                                  maxEIC = 2,
                                  maxQC = 1)

    write.table(Key.list[[1]], file = "EIC_index_pos.txt", sep = "\t", row.names = TRUE, quote = FALSE)
    write.table(Key.list[[2]], file = "EIC_index_neg.txt", sep = "\t", row.names = TRUE, quote = FALSE)

  } else {

    #Combines the two ion mode tables with compounds into one for EIC plotting
    Key.pos <- read.table(file = "EIC_index_pos.txt", sep = "\t", header = TRUE)
    Key.neg <- read.table(file = "EIC_index_neg.txt", sep = "\t", header = TRUE)
    Key.list <- list(Key.pos,Key.neg)

    Peak.list.combined <- remove_ion_dup(Peak.list = NULL,
                                         Key.list = Key.list,
                                         tbl.id = from.tables[1],
                                         peak.db = peak_db)

    write_tbl(Peak.list.combined,
              peak.db = peak_db,
              myname = to.table)


  }


}

#' @title Combine Peaklists simply
#'
#' @export
#' @description Combines Peaklists from positive and negative ionization modes
#'   without removing ion mode duplicates.
#' @param from.tables character vector of table names to draw from databases to
#'   be combined simply.
#' @param to.table to which table should LUMA save the modified Peaklist
#' @param peak.db what database contains the Peaklists to be combined simply.
#'   Default is \code{"Peaklist_db"}
#' @param db.dir directory containing the database. Default is 'db'
#' @return the file \code{"Peaklist_Combined.csv"} is written to the working
#'   directory
#' @examples
#' \dontrun{
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'   db.dir <- system.file("extdata/", package = "lcmsfishdata")
#'   InitWorkflow(db.dir = db.dir)
#'   SimplyPeaklists(from.tables = c("Peaklist_Pos_Solvent Peaks Removed",
#'   "Peaklist_Neg_Solvent Peaks Removed"),
#'   to.table = "Peaklist_Combined",
#'   peak.db = peak_db, db.dir = db.dir)
#'   FinalWorkflow(peak_db = peak_db, lib_db = lib_db)
#'   }
#' }
SimplyPeaklists <- function(from.tables,to.table,peak.db,db.dir){

  #Set default values
  if(missing(peak.db))
    peak.db = "Peaklist_db"
  if(missing(db.dir))
    db.dir = "db"

  #Update Peaklist database connection globally
  peak_db <<- peak_db <- connect_peakdb(file.base = peak.db,db.dir = db.dir)

  Peak.list.combined <- pre_combine_ion_modes(Peak.list = NULL,
                                              tbl.id = c(from.tables[1],
                                                         from.tables[2]),
                                              peak.db = peak_db)

  write_tbl(Peak.list.combined,
            peak.db = peak_db,
            myname = to.table)


  write.table(Peak.list.combined, file = "Peaklist_Combined.csv", sep = ",", row.names = FALSE)

}
