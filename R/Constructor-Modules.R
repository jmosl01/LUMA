#' @title Constructor function to initiate LUMA Workflow
#'
#' @export
#' @description All LUMA workflows must start with this function. Creates the first Peaklist and sets up storing and passing Peaklists and ancillary data between modules.
#' @param ion.id character vector specifying identifier in mz data filenames designating positive or negative ionization or both.
#' Default is c('Pos','Neg')
#' @param blanks.dir character name of subdirectory containing process blank data files.
#' Default is 'Blanks'
#' @param db.dir character name of subdirectory to store databases
#' Default is 'db'
#' Positive identifier must come first. Default is c('Pos','Neg')
#' @param adduct.files character vector specifying file names of csv files containing adduct rules, whose file structure is inherited from CAMERA.
#' Positive mode adducts must come first. Default is c("primary_adducts_pos.csv","primary_adducts_neg.csv")
#' @param use.CAMERA logical indicating whether to use existing CAMERA object in global environment.
#' Default is to look for CAMERA objects saved by previous calls to this function and run CAMERA if missing.
#' @param use.XCMS logical indicating whether to use existing XCMS object in global environment.
#' Default is to look for XCMS objects saved by previous calls to this function and run XCMS if missing.
#' @param CAMERA.obj which CAMERA object to use to initialize LUMA workflow.
#' Only relevant if use.CAMERA is TRUE
#' @param XCMS.obj which XCMS object to use to initialize LUMA workflow.
#' Only relevant if use.XCMS is TRUE
#' @param graph.method graphing method to use for CAMERA.
#' Default is 'lpc'. See CAMERA documentation for details.
#' @param ion.modes which ion modes will be processed in the workflow. Must be 'Positive', 'Negative', or both.
#' Default is both; i.e. c('Positive','Negative').
#' @return global variables are returned
#' @importFrom utils read.csv
InitWorkflow <- function(ion.id,blanks.dir,db.dir,adduct.files,use.CAMERA,use.XCMS,CAMERA.obj,XCMS.obj,graph.method) {

  #Initialize hidden global variables for tracking
  .LUMAmsg <- ""

  #Initialize all other global variables
  BLANK <- NULL
  opt.dir <- NULL
  mzdatapath <- NULL
  ion.mode <- NULL
  ppm.cutoff <- NULL
  rt.cutoff <- NULL
  Voidrt <- NULL
  Corr.stat.pos <- NULL
  Corr.stat.neg <- NULL
  cv.cutoff <- NULL
  mf.cutoff <- NULL
  Endogenous.thresh <- NULL
  Solvent.ratio <- NULL
  gen.plots <- NULL
  keep.singletons <- NULL
  DataFiles <- NULL
  Sexes <- NULL
  Classes <- NULL
  no.Samples <- NULL
  Endogenous <- NULL
  rules <- NULL
  peak_db <- NULL
  ion.modes <- NULL



  #Set default values for constructor function arguments
  if(missing(ion.id))
    ion.id <- c("Pos","Neg")
  if(missing(adduct.files))
    adduct.files <- c("primary_adducts_pos.csv","primary_adducts_neg.csv")
  if(missing(blanks.dir))
    blanks.dir <- "Blanks"
  if(missing(graph_method))
    graph.method <- "lpc"
  if(missing(db.dir))
    db.dir <- "db"
  if(missing(use.CAMERA))
    use.CAMERA = FALSE
  if(missing(use.XCMS))
    use.XCMS = FALSE
  if(missing(ion.modes))
    ion.modes = c("Positive","Negative")

  #Set Script Info globally
  if(file.exists("Script Info.txt")) {
    .LUMAmsg <<- "Initiating LUMA Workflow!\n\n"
    cat(.LUMAmsg)
    mydir <- read.table(file = "Script Info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    opt.dir <<- mydir[1,1]
    mzdatapath <<- mydir[1,2]
    ion.mode <<- mydir[1,3]
    BLANK <<- as.logical(mydir[1,4])

  } else {
    if(is.null(opt.dir) || is.null(mzdatapath) || is.null(ion.mode) || is.null(BLANK)) {
      .set_LUMAerror("Please place \"Script Info.txt\" into your working directory. \nAlternatively, you should set the opt.dir, mzdatapath, ion.mode and BLANK arguments.\n\n")
    }
  }

  #Set mzdatafiles, adduct rules, and ion.modes globally
  DataFiles <<- .get_DataFiles(mzdatapath,ion.mode,BLANK,ion.id)
  rules <<- .get_rules(ion.mode,adduct.files)
  ion.modes <<- ion.modes

  #Set search parameters globally
  if(file.exists("Search Parameters.txt")) {
  .set_LUMAmsg("Setting search parameters globally.\n\n")
    search.par <- read.table(file = "Search Parameters.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
    ppm.cutoff <<- search.par[,1]
    rt.cutoff <<- search.par[,2]
    Voidrt <<- search.par[,3]
    Corr.stat.pos <<- search.par[,4]
    Corr.stat.neg <<- search.par[,5]
    cv.cutoff <<- search.par[,6]
    mf.cutoff <<- search.par[,7]
    Endogenous.thresh <<- search.par[,8]
    Solvent.ratio <<- search.par[,9]
    gen.plots <<- search.par[,10]
    keep.singletons <<- search.par[,11]

  } else .set_LUMAmsg("Search parameters not set globally for LUMA. \nYou must set search parameters manually for all modules in this workflow.\n\n")

  #Set sample class info globally
  if(file.exists("Sample Class.txt"))  {
  .set_LUMAmsg("Setting sample class info globally.\n\n")
    Sample.df <- read.table(file = "Sample Class.txt", sep = "\t", header = TRUE)
    Sexes <<- Sample.df[,1]
    Classes <<- Sample.df[,2]
    no.Samples <<- Sample.df[,3]
    Endogenous <<- Sample.df[,4]
      } else {
        if(is.null(Sexes) || is.null(Classes) || is.null(no.Samples) || is.null(Endogenous)) {
          .set_LUMAerror("Please place \"Sample Class.txt\" into your working directory. \nAlternatively, you should set the Sex, Class, no.Samples and Endogenous arguments.\n\n")
      }
    }

  #Initialize SQLite database connection
  peak_db <<- connect_peakdb(DataFiles,db.dir)

  ##Check for existing XCMS and CAMERA objects. If not specified, check for saved XCMS and CAMERA objects.
  ##If none exist, runs XCMS and CAMERA.
  PreProcesslist <- .set_PreProcessFileNames(ion.mode,BLANK)
  XCMS.file <- PreProcesslist[[1]]
  CAMERA.file <- PreProcesslist[[2]]

  #XCMS sanity check
  if(use.XCMS) {
    if(missing(XCMS.obj)) {
      .set_LUMAerror("You must set XCMS.obj if use.XCMS is true.\n\n")
    } else {
      XCMS.obj <- .xcmsSanityCheck(XCMS.obj)
    }
  }

  #CAMERA sanity check
  if(use.CAMERA) {
    if(missing(CAMERA.obj)) {
      .LUMAmsg <<- "You must set CAMERA.obj if use.CAMERA is true. \nSee the LUMA vignette for details.\n\n"
      cat(.LUMAmsg)
    } else {
      CAMERA.obj <- .CAMERASanityCheck(CAMERA.obj)
    }
  }

  #Pre-process DataFiles
  .PreProcess_Files(XCMS.file,CAMERA.file)
}

#' @title Constructor function to calculate minimum fractions across samples
#'
#' @export
#' @description Modifies the Peak.list by calculating the minimum number of samples within a treatment class to contain a given feature.
#' See calc_minfrac for more details
#' @param from.table from which table should LUMA pull the Peak.list
#' @param to.table to which should LUMA save the modified Peak.list
AddMinFrac <- function(from.table,to.table) {
  .set_LUMAmsg("Adding Mininum Fraction values to Peaklist.\n\n")

  Peak.list <- read_tbl(mytable = from.table,
                        peak.db = peak_db)

  #Calculate Minfrac for sample classes
  new.Peak.list <- calc_minfrac(Sample.df = data.frame(Sex = Sex,
                                                       Class = Class,
                                                       n = no.Samples,
                                                       Endogenous = Endogenous),
                                xset4 = xset4,
                                BLANK = BLANK,
                                Peak.list = Peak.list)
  write_tbl(mydf = new.Peak.list,
            peak.db = peak_db,
            myname = to.table)
}

#' @title Constructor function to finalize LUMA workflow
#'
#' @export
#' @description All LUMA workflows must end with this function.  Records the names of existing data tables and writes out the LUMA log file for traceability.
#' @param my.log name of log file to write.
#' Default is 'LUMA_log.txt'
#' @param peak_db existing peak database connection
#' Default is to look for existing database connections.
#' @param lib_db existing library database connection
#' Default is to look for existing database connections.
#' @return disconnected SQLite connections and written log file
#' @importFrom DBI dbListTables dbDisconnect
FinalWorkflow <- function(my.log,peak_db,lib_db) {

  #Set default values
  if(missing(my.log))
    my.log <- "LUMA_log.txt"

  if(missing(peak_db)) {
    if(is.null(peak_db))
      .set_LUMAmsg("No peak database connection exists. Therefore did not close database connection.\n\n")
  } else {
    peak.tbls <- dbListTables(peak_db)
    .set_LUMAmsg("Peak database contains the following tables:")
    msg <- paste(peak.tbls, collapse = '\n')
    .set_LUMAmsg(msg)
    .set_LUMAmsg("\n\nClosing the Peak database connection.")
    dbDisconnect(peak_db)
  }

  if(missing(lib_db)) {
    if(is.null(lib_db))
      .set_LUMAmsg("No library database connection exists. Therefore did not close database connection.\n\n")
  } else {
    lib.tbls <- dbListTables(lib_db)
    .set_LUMAmsg("Library database contains the following tables:")
    msg <- paste(lib.tbls, collapse = '\n')
    .set_LUMAmsg(msg)
    .set_LUMAmsg("\n\nClosing the Library database connection.")
  }


sink(file = my.log)
writeLines(.LUMAmsg, sep = "\n")
sink()

}
