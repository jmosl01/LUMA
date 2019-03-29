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
#' @param mytable character name of the first Peak.list table
#' Default is 'From CAMERA'
#' @param calc.minfrac logical should LUMA calculate the minimum fraction values for the initial Peak.list
#' Default is TRUE
#' @return global variables and Peaklist in database are returned
#' @importFrom utils read.csv
InitWorkflow <- function(ion.id,blanks.dir,db.dir,adduct.files,use.CAMERA,use.XCMS,CAMERA.obj,XCMS.obj,
                         graph.method,ion.modes,mytable,calc.minfrac) {

  #Initialize all global variables
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
  Name <- NULL
  Formula <- NULL
  Molecular.Weight <- NULL
  RT..Min. <- NULL
  XCMS.par <- NULL
  CAMERA.par <- NULL
  CT.ID <- NULL
  Plate.Number <- NULL
  Plate.Position <- NULL
  Sample.phenodata <- NULL
  Library.phenodata <- NULL

  #Set default values for constructor function arguments
  if(missing(ion.id))
    ion.id <- c("Pos","Neg")
  if(missing(adduct.files))
    adduct.files <- c("primary_adducts_pos.csv","primary_adducts_neg.csv")
  if(missing(blanks.dir))
    blanks.dir <- "Blanks"
  if(missing(graph.method))
    graph.method <- "lpc"
  if(missing(db.dir))
    db.dir <- "db"
  if(missing(use.CAMERA))
    use.CAMERA <- FALSE
  if(missing(use.XCMS))
    use.XCMS <- FALSE
  if(missing(ion.modes))
    ion.modes <- c("Positive","Negative")
  if(missing(mytable))
    mytable <- "From CAMERA"
  if(missing(calc.minfrac))
    calc.minfrac <- TRUE

  #Set Script Info globally
  if(file.exists("Script Info.txt")) {
    cat("Initiating LUMA Workflow!\n\n")
    mydir <- read.table(file = "Script Info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    opt.dir <<- opt.dir <- mydir[1,1]
    mzdatapath <<- mzdatapath <- mydir[1,2]
    ion.mode <<- ion.mode <- mydir[1,3]
    BLANK <<- BLANK <- as.logical(mydir[1,4])

  } else {
    if(is.null(opt.dir) || is.null(mzdatapath) || is.null(ion.mode) || is.null(BLANK)) {
      stop("Please place \"Script Info.txt\" into your working directory. \nAlternatively, you should set the opt.dir, mzdatapath, ion.mode and BLANK arguments.\n\n")
    }
  }

  #Set metadata globally
  DataFiles <<- DataFiles <- .get_DataFiles(mzdatapath,ion.mode,BLANK,ion.id,blanks.dir)
  rules <<- rules <- .get_rules(ion.mode,adduct.files)
  ion.modes <<- ion.modes
  ion.id <<- ion.id

  #Set search parameters globally
  if(file.exists("Search Parameters.txt")) {
  cat("Setting search parameters globally.\n\n")
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

  } else cat("Search parameters not set globally for LUMA. \nYou must set search parameters manually for all modules in this workflow.\n\n")

  #Set sample class info globally
  if(file.exists("Sample Class.txt"))  {
    cat("Setting sample class info globally.\n\n")
    Sample.df <- read.table(file = "Sample Class.txt", sep = "\t", header = TRUE,
                            colClasses = c("character","character","numeric","logical"))
    Sexes <<- Sexes <- Sample.df[,"Sex"]
    Classes <<- Classes <- Sample.df[,"Class"]
    no.Samples <<- no.Samples <- Sample.df[,"n"]
    Endogenous <<- Endogenous <- Sample.df[,"Endogenous"]
  } else {
    if(is.null(Sexes) || is.null(Classes) || is.null(no.Samples) || is.null(Endogenous)) {
      stop("Please place \"Sample Class.txt\" into your working directory. \nAlternatively, you should set the Sex, Class, no.Samples and Endogenous arguments.\n\n")
    }
  }

  #Set sample phenotype data globally
  if(file.exists("Sample Data.csv"))  {
    cat("Setting sample phenotype data globally.\n\n")
    Sample.data <- read.table(file = "Sample Data.csv", sep = "," , header = TRUE, stringsAsFactors = FALSE)
    Sample.data <- Sample.data[order(Sample.data[,"CT.ID"]),]
    CT.ID <<- CT.ID <- Sample.data[,"CT.ID"]
    Plate.Number <<- Plate.Number <- Sample.data[,"Plate.Number"]
    Plate.Position <<- Plate.Position <- Sample.data[,"Plate.Position"]
    Sample.phenodata <<- Sample.phenodata <- Sample.data[,-which(colnames(Sample.data) %in% c("CT.ID","Plate.Number","Plate.Position"))]
  } else {
    if(is.null(CT.ID) || is.null(Plate.Number) || is.null(Plate.Position)) {
      stop("Please place \"Sample Data.csv\" into your working directory. \nAlternatively, you should set the CT-ID, Plate.Number, and Plate.Position arguments.\n\n")
    }
  }

  #Set Annotated Library info globally
  if(file.exists("Annotated Library.csv"))  {
    cat("Setting Annotated Library info globally.\n\n")
    Annotated.Library <- read.csv(file = "Annotated Library.csv", sep = ",", fill = TRUE, header = TRUE)
    Name <<- Name <- Annotated.Library[,"Name"]
    Formula <<- Formula <- Annotated.Library[,"Formula"]
    Molecular.Weight <<- Molecular.Weight <- Annotated.Library[,"Molecular.Weight"]
    RT..Min. <<- RT..Min. <- Annotated.Library[,"RT..Min."]
    Library.phenodata <<- Library.phenodata <- Annotated.Library[,-which(colnames(Annotated.Library) %in% c("Name","Formula","Molecular.Weight","RT..Min."))]
  } else {
    if(is.null(Name) || is.null(Formula) || is.null(Molecular.Weight) || is.null(RT..Min.)) {
      stop("Please place \"Sample Class.txt\" into your working directory. \nAlternatively, you should set the Name, Formula, Molecular.Weight and RT..Min. arguments.\n\n")
    }
  }


  #Initialize SQLite database connections globally
  file.base <- gen_filebase(DataFiles,BLANK,ion.id,ion.mode)
  peak_db <<- peak_db <- connect_peakdb(file.base,db.dir)

  ##Check for existing XCMS and CAMERA objects. If not specified, check for saved XCMS and CAMERA objects.
  ##If none exist, runs XCMS and CAMERA.
  PreProcesslist <- .set_PreProcessFileNames(ion.mode,BLANK)
  XCMS.file <- PreProcesslist[[1]]
  CAMERA.file <- PreProcesslist[[2]]

  ##Set XCMS parameters globally
  if(ion.mode == "Positive"){
    XCMS.par <<- XCMS.par <- read.table(file = "Best XCMS parameters_positive.csv", sep = "," , header = TRUE)
  } else {
    if(ion.mode == "Negative"){
      XCMS.par <<- XCMS.par <- read.table(file = "Best XCMS parameters_negative.csv", sep = "," , header = TRUE)
    }
  }

  #XCMS sanity check
  if(use.XCMS) {
    if(missing(XCMS.obj)) {
      stop("You must set XCMS.obj if use.XCMS is true.\n\n")
    } else {
      XCMS.obj <- .xcmsSanityCheck(XCMS.obj)
    }
  }

  ## Set CAMERA parameters globally
  if(ion.mode == "Positive"){
    CAMERA.par <<- CAMERA.par <- read.table(file = paste(opt.dir,"/Best CAMERA parameters_positive.csv", sep = ""), sep = "," , header = TRUE)
  } else {
    if(ion.mode == "Negative"){
      CAMERA.par <<- CAMERA.par <- read.table(file = paste(opt.dir,"/Best CAMERA parameters_negative.csv",sep = ""), sep = "," , header = TRUE)
    }
  }

  #CAMERA sanity check
  if(use.CAMERA) {
    if(missing(CAMERA.obj)) {
      cat("You must set CAMERA.obj if use.CAMERA is true. \nSee the LUMA vignette for details.\n\n")
    } else {
      CAMERA.obj <- .CAMERASanityCheck(CAMERA.obj,CAMERA.file)
    }
  }

  #Pre-process DataFiles
  xset4 <- .PreProcess_Files(XCMS.file,CAMERA.file,mytable)

  if(calc.minfrac) {
    ##Add minimum fraction to Peak.list
    cat("Adding Mininum Fraction values to Peaklist.\n\n")

    Peak.list <- read_tbl(mytable = mytable,
                          peak.db = peak_db)

    #Calculate Minfrac for sample classes
    new.Peak.list <- calc_minfrac(Sample.df = data.frame(Sex = Sexes,
                                                         Class = Classes,
                                                         n = no.Samples,
                                                         Endogenous = Endogenous),
                                  xset4 = xset4,
                                  BLANK = BLANK,
                                  Peak.list = Peak.list)

    write_tbl(mydf = new.Peak.list,
              peak.db = peak_db,
              myname = paste(mytable,"with Minfrac", sep = "_"))

  }
}

#' @title Constructor function to finalize LUMA workflow
#'
#' @export
#' @description All LUMA workflows must end with this function.  Records the names of existing data tables and writes out the LUMA log file for traceability.
#' @param peak_db existing peak database connection
#' @param lib_db existing library database connection
#' @return NULL
#' @importFrom DBI dbListTables dbDisconnect
FinalWorkflow <- function(peak_db,lib_db) {

  #Set default values
  if(missing(peak_db)) {
    peak_db <- NULL
  }
  if(missing(lib_db)) {
    lib_db <- NULL
  }


  #Database connection sanity check
  if(is.null(peak_db)) {
      cat("No peak database connection provided. Therefore did not close database connection.\n\n")
  } else {
    test <- dbConnect(peak_db)
    if(class(test)[1] != "SQLiteConnection") {
      stop(paste("peak_db is of class ",class(test)[1],", but needs to be of class \"SQLiteConnection\".",sep = ""))
    } else {
      peak.tbls <- dbListTables(test)
      cat("Peak database contains the following tables:\n")
      msg <- paste(peak.tbls, collapse = '\n')
      cat(msg)
      cat("\n\nClosing the Peak database connection.\n\n")
      dbDisconnect(test)
      dbDisconnect(peak_db)
      mylist <- deparse(substitute(peak_db))
    }
  }

  #Library connection sanity check
  if(is.null(lib_db)) {
      cat("No library database connection provided. Therefore did not close database connection.\n\n")
  } else {
    test <- dbConnect(lib_db)
    if(class(test)[1] != "SQLiteConnection") {
      stop(paste("peak_db is of class ",class(test)[1],", but needs to be of class \"SQLiteConnection\".",sep = ""))
    } else {
      lib.tbls <- dbListTables(test)
      cat("Library database contains the following tables:\n")
      msg <- paste(lib.tbls, collapse = '\n')
      cat(msg)
      cat("\n\nClosing the Library database connection.\n\n")
      dbDisconnect(test)
      dbDisconnect(lib_db)
      mylist <- c(mylist,deparse(substitute(lib_db)))
    }
  }

  #Clean up the database connections in the global environment
  rm(list=mylist,envir = .GlobalEnv)
}
