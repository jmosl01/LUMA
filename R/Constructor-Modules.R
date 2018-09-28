#' @title Constructor function to initiate LUMA Workflow
#'
#' @export
#' @description Creates the first Peak.list and sets up storing and passing Peak.lists and ancillary data between modules.
#' @param ion.id character vector specifying identifier in mz data filenames designating positive or negative ionization or both.
#' Default is c('Pos','Neg')
#' @param blanks.dir character name of subdirectory containing process blank data files.
#' Default is 'Blanks'
#' @param db.dir character name of subdirectory to store databases
#' Default is 'db'
#' Positive identifier must come first. Default is c('Pos','Neg')
#' @param adduct.files character vector specifying file names of csv files containing adduct rules, whose file structure is inherited from CAMERA.
#' Positive mode adducts must come first. Default is c("primary_adducts_pos.csv","primary_adducts_neg.csv")
#' @param use.CAMERA logical indicating whether to run CAMERA.
#' Default is to look for saved CAMERA objects and run CAMERA if missing.
#' @param use.XCMS logical indicating whether to run XCMS.
#' Default is to look for saved XCMS objects and run XCMS if missing.
#' @param graph.method graphing method to use for CAMERA.
#' Default is 'lpc'. See CAMERA documentation for details.
#' @return global variables are returned
#' @importFrom utils read.csv
start_LUMAWorkflow = function(ion.id,blanks.dir,db.dir,adduct.files,use.CAMERA,use.XCMS,graph.method) {

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
  file.base <- NULL
  Sexes <- NULL
  Classes <- NULL
  n <- NULL
  Endogenous <- NULL
  rules <- NULL
  peak_db <- NULL



  #Set default values for constructor function arguments
  if(missing(ion.id))
    ion.id <- c("Pos","Neg")
  if(missing(files))
    files <- c("primary_adducts_pos.csv","primary_adducts_neg.csv")
  if(missing(blanks.dir))
    blanks.dir <- "Blanks"
  if(missing(graph_method))
    graph.method <- "lpc"
  if(missing(db.dir))
    db.dir <- "db"

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
      .LUMAmsg <<- "Please place \"Script Info.txt\" into your working directory. \nAlternatively, you should set the opt.dir, mzdatapath, ion.mode and BLANK arguments. \nSee the LUMA vignette for details.\n\n"
      cat(.LUMAmsg)
      stop()
    }
  }

  #Set mzdatafiles and adduct rules globally
  file.base <<- .get_DataFiles(mzdatapath,ion.mode,BLANK,ion.id)
  rules <<- .get_rules(ion.mode,files)

  #Set search parameters globally
  if(file.exists("Search Parameters.txt")) {
    msg <- "Setting search parameters globally.\n\n"
    .LUMAmsg <<- c(.LUMAmsg,msg)
    cat(msg)
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

  } else {
    msg <- "Search parameters not set globally for LUMA. \nYou must set search parameters for all modules in this workflow.\n\n"
    .LUMAmsg <<- c(.LUMAmsg,msg)
    cat(msg)
  }

  #Set sample class info globally
  if(file.exists("Sample Class.txt"))  {
    msg <- "Setting sample class info globally.\n\n"
    .LUMAmsg <<- c(.LUMAmsg,msg)
    cat(msg)
    Sample.df <- read.table(file = "Sample Class.txt", sep = "\t", header = TRUE)
    Sexes <<- Sample.df[,1]
    Classes <<- Sample.df[,2]
    n <<- Sample.df[,3]
    Endogenous <<- Sample.df[,4]
      } else {
        if(is.null(Sexes) || is.null(Classes) || is.null(n) || is.null(Endogenous)) {
          .LUMAmsg <<- "Please place \"Sample Class.txt\" into your working directory. \nAlternatively, you should set the Sex, Class, n and Endogenous arguments. \nSee the LUMA vignette for details.\n\n"
          cat(.LUMAmsg)
          stop()
      }
    }

  #Initialize SQLite database connection
  peak_db <<- connect_peakdb(file.base,db.dir)

}

