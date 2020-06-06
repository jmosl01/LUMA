#' @title Constructor module to initiate LUMA Workflow
#'
#' @export
#' @description All LUMA workflows must start with this module Creates the
#'   first Peaklist and sets up storing and passing Peaklists and ancillary data
#'   between modules. For working examples, see \code{InitWorkflow, AnnotatePeaklist,
#'   CombineFeatures, CombinePeaklists, CullBackground, CullMF, CullCV,
#'   CullVoidVolume, FormatForMetaboAnalystR, FormatForSIMCA,
#'   NormalizePeaklists, ParseCAMERA, SimplyPeaklists, FinalWorkflow}.
#' @param ion.id character vector specifying identifier in mzdata filenames
#'   designating positive or negative ionization or both. Must be a
#'   (case-insensitive) abbreviation of the ionization mode name. Positive
#'   identifier must come first. Default is \code{c("Pos","Neg")}.
#' @param db.dir character name of subdirectory to store databases. Default is
#'   \code{"db"}
#' @param use.CAMERA logical indicating whether to use existing CAMERA object in
#'   global environment. Default is to look for CAMERA objects saved by previous
#'   calls to this function and run CAMERA if missing.
#' @param use.XCMS logical indicating whether to use existing XCMS object in
#'   global environment. Default is to look for XCMS objects saved by previous
#'   calls to this function and run XCMS if missing.
#' @param CAMERA.obj which CAMERA object to use to initialize LUMA workflow.
#'   Only relevant if \code{use.CAMERA == TRUE}.
#' @param XCMS.obj which XCMS object to use to initialize LUMA workflow. Only
#'   relevant if \code{use.XCMS == TRUE}.
#' @param graph.method graphing method to use for CAMERA. Default is
#'   \code{"lpc"}. See CAMERA documentation for details.
#' @param QC.id character vector specifying identifier in filename designating a
#'   Pooled QC sample.  Only the first value will be used.  Default is
#'   \code{"Pooled_QC_"}
#' @param ion.mode which ion mode(s) will be processed for this data. Must be
#'   one or both of \code{c("Positive","Negative")}. Default is both.
#' @param mytable character name of the first Peaklist table in the database.
#'   Default is \code{"From_CAMERA"}.
#' @param calc.minfrac logical should LUMA calculate the minimum fraction values
#'   for the initial Peaklist. Default is \code{TRUE}.
#' @param multiple logical should multiple fields be allowed in dialog boxes.
#'   Default is \code{FALSE}.
#' @return global variables and Peaklist in database are returned
#' @importFrom utils read.csv
#' @examples
#' \dontrun{
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'   db.dir <- system.file("extdata/", package = "lcmsfishdata")
#'   InitWorkflow(db.dir = db.dir)
#'   AnnotatePeaklist(from.table = "From CAMERA", to.table = "Annotated")
#'   FinalWorkflow(peak_db = peak_db)
#'   }
#' }
InitWorkflow <- function(ion.id,db.dir,use.CAMERA,use.XCMS,CAMERA.obj,XCMS.obj,
                         graph.method,QC.id,ion.mode,mytable,calc.minfrac,multiple) {

  #Initialize all global variables
  BLANK <- NULL
  opt.dir <- NULL
  IonMode <- NULL
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
  mzdatafiles <- NULL
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
  isdb <- TRUE

  #Set default values for constructor function arguments
  if(missing(ion.id))
    ion.id <- c("Pos","Neg")
  if(missing(graph.method))
    graph.method <- "lpc"
  if(missing(QC.id))
    QC.id <- "Pooled_QC_"
  if(missing(db.dir)) {
    db.dir <- "db"
    isdb <- FALSE
  }
  if(missing(use.CAMERA))
    use.CAMERA <- FALSE
  if(missing(use.XCMS))
    use.XCMS <- FALSE
  if(missing(ion.mode))
    ion.mode <- c("Positive","Negative")
  if(missing(mytable))
    mytable <- "From_CAMERA"
  if(missing(calc.minfrac))
    calc.minfrac <- TRUE
  if(missing(multiple))
    multiple <- FALSE

  #Set Script Info globally

  #Initiate Dialog Boxes
  script_dlg <- ScriptInfo_dlg(multiple = multiple, isdb, db.dir)
  BLANK = script_dlg$BLANK
  IonMode = script_dlg$IonMode

  DataFiles <- .get_DataFiles(mzdatapath = script_dlg$DataDir,
                                            BLANK = BLANK,
                                            IonMode = IonMode,
                                            ion.id = ion.id,
                                            ion.mode = ion.mode)

  input_dlg <- InputFiles_dlg(WorkingDir = script_dlg$WorkingDir, multiple = multiple)

  rules <- .get_rules(adduct.file = input_dlg$Adducts)


  #Set metadata parameters globally
  DataFiles <<- DataFiles
  rules <<- rules
  ion.mode <<- ion.mode
  ion.id <<- ion.id
  QC.id <<- QC.id
  BLANK <<- BLANK
  IonMode <<- IonMode

  #Set search parameters globally
  if(file.exists(input_dlg$SearchPar)) {
  cat("Setting search parameters globally.\n\n")
    search.par <- read.table(file = input_dlg$SearchPar, sep = "\t", stringsAsFactors = FALSE, header = TRUE)
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
  if(file.exists(input_dlg$SampleClass))  {
    cat("Setting sample class info globally.\n\n")
    Sample.df <- read.table(file = input_dlg$SampleClass, sep = "\t", header = TRUE,
                            colClasses = c("character","character","numeric","logical"))
    Sexes <<- Sexes <- Sample.df[,"Sex"]
    Classes <<- Classes <- Sample.df[,"Class"]
    no.Samples <<- no.Samples <- Sample.df[,"n"]
    Endogenous <<- Endogenous <- Sample.df[,"Endogenous"]
  } else {
    if(is.null(Sexes) || is.null(Classes) || is.null(no.Samples) || is.null(Endogenous)) {
      stop("Please place \"Sample_Class.txt\" into your working directory. \nAlternatively, you should set the Sex, Class, no.Samples and Endogenous arguments.\n\n")
    }
  }

  #Set sample phenotype data globally
  if(file.exists(input_dlg$SampleData))  {
    cat("Setting sample phenotype data globally.\n\n")
    Sample.data <- read.table(file = input_dlg$SampleData, sep = "," , header = TRUE, stringsAsFactors = FALSE)
    Sample.data <- Sample.data[order(Sample.data[,"CT.ID"]),]
    CT.ID <<- CT.ID <- Sample.data[,"CT.ID"]
    Plate.Number <<- Plate.Number <- Sample.data[,"Plate.Number"]
    Plate.Position <<- Plate.Position <- Sample.data[,"Plate.Position"]
    Sample.phenodata <<- Sample.phenodata <- Sample.data[,-which(colnames(Sample.data) %in% c("CT.ID","Plate.Number","Plate.Position"))]
  } else {
    if(is.null(CT.ID) || is.null(Plate.Number) || is.null(Plate.Position)) {
      stop("Please place \"Sample_Data.csv\" into your working directory. \nAlternatively, you should set the CT-ID, Plate.Number, and Plate.Position arguments.\n\n")
    }
  }

  #Set Annotated Library info globally
  if(file.exists(input_dlg$AnnotatedLibrary))  {
    cat("Setting Annotated Library info globally.\n\n")
    Annotated.Library <- read.csv(file = input_dlg$AnnotatedLibrary, sep = ",", fill = TRUE, header = TRUE)
    Name <<- Name <- Annotated.Library[,"Name"]
    Formula <<- Formula <- Annotated.Library[,"Formula"]
    Molecular.Weight <<- Molecular.Weight <- Annotated.Library[,"Molecular.Weight"]
    RT..Min. <<- RT..Min. <- Annotated.Library[,"RT..Min."]
    Library.phenodata <<- Library.phenodata <- Annotated.Library[,-which(colnames(Annotated.Library) %in% c("Name","Formula","Molecular.Weight","RT..Min."))]
  } else {
    if(is.null(Name) || is.null(Formula) || is.null(Molecular.Weight) || is.null(RT..Min.)) {
      stop("Please place \"Sample_Class.txt\" into your working directory. \nAlternatively, you should set the Name, Formula, Molecular.Weight and RT..Min. arguments.\n\n")
    }
  }


  #Initialize SQLite database connections globally
  file.base <- gen_filebase(DataFiles,BLANK,ion.id,IonMode)
  peak_db <<- peak_db <- connect_peakdb(file.base,db.dir)

  ##Check for existing XCMS and CAMERA objects. If not specified, check for saved XCMS and CAMERA objects.
  ##If none exist, runs XCMS and CAMERA.
  XCMS.file <- input_dlg$XCMSObj
  CAMERA.file <- input_dlg$CAMERAObj

  ##Set XCMS parameters globally
  temp_xcms <- grep(".csv",XCMS.file)

  if(length(temp_xcms) == 0){

    XCMS.par <<- XCMS.par <- input_dlg$XCMS.par

  } else {
    if(length(temp_xcms) == 1){
      if(file.exists(XCMS.file)) {
        XCMS.par <<- XCMS.par <- read.table(XCMS.file, sep = ",", header = TRUE)
      }

    } else {
      if(length(temp_xcms) > 1) {
        stop("Error: Does your XCMS parameters file have too many file extensions?")
      }
    }
  }

  #XCMS sanity check
  if(use.XCMS) {
    if(missing(XCMS.obj)) {
      stop("You must set XCMS.obj if use.XCMS is true. \nSee the LUMA vignette for details.\n\n")
    } else {
      XCMS.obj <- .xcmsSanityCheck(XCMS.obj)
    }
  }

  #CAMERA sanity check
  if(use.CAMERA) {
    if(missing(CAMERA.obj)) {
      stop("Error: You must set CAMERA.obj if use.CAMERA is true. \nSee the LUMA vignette for details.\n\n")
    } else {
      CAMERA.obj <- .CAMERASanityCheck(CAMERA.obj,CAMERA.file)
    }
  }

  #Pre-process DataFiles
  xset4 <- .PreProcess_Files(XCMS.file = XCMS.file,
                             CAMERA.file = CAMERA.file,
                             mytable = mytable,
                             file.base = file.base,
                             IonMode = IonMode)

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

#' @title Constructor module to finalize LUMA workflow
#'
#' @export
#' @description All LUMA workflows must end with this module. Records the names
#'   of existing data tables and writes out the LUMA log file for traceability.
#'   For working examples, see \code{InitWorkflow, AnnotatePeaklist,
#'   CombineFeatures, CombinePeaklists, CullBackground, CullMF, CullCV,
#'   CullVoidVolume, FormatForMetaboAnalystR, FormatForSIMCA,
#'   NormalizePeaklists, ParseCAMERA, SimplyPeaklists, FinalWorkflow}.
#' @param peak_db existing peak database connection
#' @param lib_db existing library database connection
#' @return NULL
#' @importFrom DBI dbListTables dbDisconnect
FinalWorkflow <- function(peak_db,lib_db) {

  #Initialize all global variables
  peak.tbls <- NULL
  lib.tbls <- NULL


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
    if(slot(lib_db, "dbname") != ":memory:") {

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

      } else {
      lib.tbls <- dbListTables(lib_db)
      cat("Library database contains the following tables:\n")
      msg <- paste(lib.tbls, collapse = '\n')
      cat(msg)
      cat("\n\nClosing the Library database connection.\n\n")
      dbDisconnect(lib_db)
      mylist <- c(mylist,deparse(substitute(lib_db)))
        }
    }

  #Set LUMA log variables globally
  peak.tbls <<- peak.tbls
  lib.tbls <<- lib.tbls

  #Clean up the database connections in the global environment
  rm(list=mylist,envir = .GlobalEnv)
}
