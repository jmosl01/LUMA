#' @title Removes background components from Peak.list
#'
#' @export
#' @description Searches features or compounds in Peak.list that are present in
#'   process blanks and removes them. Also flags compounds as exogenous, defined
#'   as mean abundance in endogenous samples below user-defined threshold.
#' @param Peak.list a named list of data frames (two per ionization mode)
#'   containing intensity matrices across all study samples and Pooled QCs and
#'   process blanks. Names should be c('pos','neg','blanks_pos','blanks_neg').
#'   Alternatively may use existing database connections by setting to NULL
#' @param Sample.df a data frame with class info as columns, containing a
#'   separate row entry for each unique sex/class combination. Must contain the
#'   columns 'Sex','Class','n','Endogenous'.
#' @param search.par a single-row data frame with 11 variables containing
#'   user-defined search parameters. Must contain the columns
#'   \code{"ppm","rt","Voidrt","Corr.stat.pos","Corr.stat.neg","CV","Minfrac","Endogenous","Solvent","gen.plots","keep.singletons"}.
#' @param method which method to apply to search for background components.  See
#'   find_Background for details.
#' @param lib.db character name of database to contain Solvent Library
#' @param tbl.id character vector of table names to draw from databases. First
#'   value should be table name from peak database, second should be table name
#'   from solvent database. Default is NULL
#' @param db.list list chracter names of databases containing results from
#'   processing positive mode (1,3) and negative mode (2,4) data for samples
#'   (1,2) and blanks (3,4) Default is NULL
#' @param db.dir character directory containing the databases Default is 'db'
#' @param new.db character what should the new database be called Default is
#'   'Peaklist_db'
#' @param mem logical should database be in-memory. Default is FALSE
#' @param ... Arguments to pass parameters to find_Background
#' @return nested list a list for each ionization mode, each containing a list
#'   of two dataframes: the first contains the intensity matrix for the peaklist
#'   with solvent peaks removed, the second contains the intensity matrix for
#'   the solvent peaks
#' @importFrom plyr llply
#' @importFrom dplyr filter
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#' file <- system.file('extdata','Sample_Class.txt', package = "LUMA")
#' Sample.df <- read.table(file, header = TRUE, sep = "\t")
#' file2 <- system.file('extdata','Search_parameters.txt', package = "LUMA")
#' search.par <- read.table(file2, header = TRUE, sep = "\t")
#' \donttest{
#'   #From m/z features
#'   Peak.list <- list(pos = lcmsfishdata::Peaklist_Pos$From_CAMERA, neg =
#'   lcmsfishdata::Peaklist_Neg$From_CAMERA, blanks_pos =
#'   lcmsfishdata::Blanks_Pos$From_CAMERA, blanks_neg =
#'   lcmsfishdata::Blanks_Neg$From_CAMERA)
#'   test <- remove_background_peaks(Peak.list = Peak.list, Sample.df =
#'   Sample.df, search.par = search.par, method = "mz", mem = TRUE)
#'   lapply(test, head) #Peaklists with removed background components are returned
#' }
#'
#' #From combined features
#' Peak.list <- list(pos = lcmsfishdata::Peaklist_Pos$Trimmed_by_MinFrac, neg =
#' lcmsfishdata::Peaklist_Neg$Trimmed_by_MinFrac, blanks_pos =
#' lcmsfishdata::Blanks_Pos$Combined_Isotopes_and_Adducts, blanks_neg =
#' lcmsfishdata::Blanks_Neg$Combined_Isotopes_and_Adducts)
#' test <- remove_background_peaks(Peak.list = Peak.list, Sample.df = Sample.df,
#' search.par = search.par, method = "monoMass", mem = TRUE)
#' lapply(test, head) #Peaklists with removed background components are returned
#'  }
remove_background_peaks = function(Peak.list = NULL, Sample.df, search.par, method, lib.db, tbl.id, db.list, db.dir, new.db, mem, ...) {

  # Set default values
  if (missing(tbl.id))
    tbl.id = NULL
  if (missing(db.list))
    db.list = NULL
  if (missing(db.dir))
    db.dir = "db"
  if (missing(mem))
    mem = FALSE
  if (missing(method))
    stop("Need to specify a method!", call. = FALSE)
  if (is.null(Peak.list) && is.null(tbl.id))
    stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)

  if (is.null(Peak.list)) {
    mydbs <- connect_lumadb(db.list,db.dir,new.db, mem)
    nm = names(mydbs)[-1]
    nm1 = nm[-grep("blanks", nm)]  # Return all sample databases
    nm2 = nm[grep("blanks", nm)]  #Return all blank databases
    list1 <- llply(nm1, .fun = function(x) read_tbl(tbl.id[1], mydbs[[x]]))
    list2 <- llply(nm2, .fun = function(x) read_tbl(tbl.id[2], mydbs[[x]]))
    names(list1) <- gsub("_db", "\\1", nm1)
    names(list2) <- gsub("_db", "\\1", nm2)
    Peak.list <- list1
    Solv.list <- list2
  } else {
    nm = names(Peak.list)
    nm1 = nm[-grep("blanks", nm)]  # Return all sample databases
    nm2 = nm[grep("blanks", nm)]  #Return all blank databases
    Solv.list <- Peak.list[nm2]
    Peak.list <- Peak.list[nm1]
  }
  lib_db <- connect_libdb(lib.db = lib.db, db.dir = db.dir, mem = mem)
  peak_db <- connect_lumadb(db.list = db.list, db.dir = db.dir, new.db = new.db, mem = mem)

  class(method) <- method
  masterlist <- find_Background(method, Peak.list, Solv.list, Sample.df, search.par, lib_db, ...)
  return(masterlist)
}



#' @title Removes void volume from Peak.list
#'
#' @export
#' @description Removes features or compounds found in the void volume of the
#'   chromatographic run from the Peak.list
#' @param Peak.list data frame. Must have Correlation.stat column.  Should
#'   contain output columns from XCMS and CAMERA, and additional columns from
#'   IHL.search, Calc.MinFrac, Calc.corr.stat and EIC.plotter functions.
#' @param search.par a single-row data frame with 11 variables containing
#'   user-defined search parameters. Must contain the columns
#'   \code{"ppm","rt","Voidrt","Corr.stat.pos","Corr.stat.neg","CV","Minfrac","Endogenous","Solvent","gen.plots","keep.singletons"}.
#' @param method which method to apply to trim by retention time.  See trim_rt
#'   for details
#' @param ... Arguments to pass to trim_rt
#' @return data frame Peak.list.trimmed original Peak.list without all
#'   metabolite groups containing at least one feature in the void volume
#' @md
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)){
#' file <- system.file('extdata','Search_Parameters.txt', package = "LUMA")
#' search.par <- read.table(file, sep = "\t", header = TRUE)
#' test <- remove_void_volume(Peak.list = lcmsfishdata::Peaklist_Pos$From_CAMERA, search.par =
#' search.par, method = "mz")
#'
#' #Removes 275 features within the void volume
#' nrow(lcmsfishdata::Peaklist_Pos$From_CAMERA) - nrow(test)
#' }
remove_void_volume = function(Peak.list, search.par, method,...) {
    void.rt <- as.numeric(search.par[1, "Voidrt"])
    class(method) <- method
    Peak.list.trimmed <- trim_rt(method,Peak.list,void.rt,...)
    return(Peak.list.trimmed)
}

#' @title Trims by CV
#'
#' @export
#' @description Removes metabolites with coefficient of variation greater than
#'   the user specified threshold, calculated from the QC samples
#' @param Peak.list data frame. Must have QC sample columns that contain the
#'   string 'Pooled_QC_'.  Should contain output columns from XCMS and CAMERA,
#'   and additional columns from IHL.search, Calc.MinFrac, CAMERA.parser,
#'   Calc.corr.stat and Combine.phenodata base functions.
#' @param search.par a single-row data frame with 11 variables containing
#'   user-defined search parameters. Must contain the columns
#'   \code{"ppm","rt","Voidrt","Corr.stat.pos","Corr.stat.neg","CV","Minfrac","Endogenous","Solvent","gen.plots","keep.singletons"}.
#' @param QC.id character vector specifying identifier in filename designating a
#'   Pooled QC sample.  Only the first value will be used.  Default is
#'   \code{"Pooled_QC_"}
#' @return data frame Peak.list.trimmed original Peak.list without all
#'   metabolite groups with coefficient of variation greater than user specified
#'   threshold; if dataset contains blanks, data.frame with all NA values is
#'   returned
#' @importFrom stats sd
#' @examples
#' library(LUMA)
#' file <- system.file('extdata','Search_Parameters.txt', package = "LUMA")
#' search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
#' test <- trim_cv(Peak.list = Peaklist_Pos$From_CAMERA, search.par = search.par)
#' nrow(Peaklist_Pos$From_CAMERA) -  nrow(test) #Equals 13
#'
#' test <- trim_cv(Peak.list = Peaklist_Pos$Combined_Isotopes_and_Adducts,
#' search.par = search.par)
#' nrow(Peaklist_Pos$Combined_Isotopes_and_Adducts) -  nrow(test) #Equals 9
trim_cv = function(Peak.list, search.par,QC.id) {

    #Set Default Values
    if(missing(QC.id))
      QC.id <- "Pooled_QC_"

    Peak.list.cv <- calc_cv(Peak.list,QC.id)
    CV.cutoff <- as.numeric(search.par[1, "CV"])
    Peak.list.trimmed <- Peak.list.cv[Peak.list.cv["X.CV"][[1]] < CV.cutoff, ]
    return(Peak.list.trimmed)
}

#' @title Trims by MinFrac
#'
#' @export
#' @description Removes metabolites with MinFrac smaller than the user specified
#'   threshold. The maximum MinFrac value is chosen from all features within a
#'   metabolite group.
#' @param object used for method dispatch. Can be any object. See usage for
#'   details
#' @param Peak.list data frame. Must have MinFrac column.  Should contain output
#'   columns from XCMS and CAMERA, and additional columns from IHL.search,
#'   Calc.MinFrac, CAMERA.parser, Calc.corr.stat and Combine.phenodata base
#'   functions.
#' @param search.par a single-row data frame with 11 variables containing
#'   user-defined search parameters. Must contain the columns
#'   \code{"ppm","rt","Voidrt","Corr.stat.pos","Corr.stat.neg","CV","Minfrac","Endogenous","Solvent","gen.plots","keep.singletons"}.
#' @return data frame Peak.list.trimmed original Peak.list containing all
#'   metabolite groups containing at least one feature that has MinFrac value
#'   greater than user specified threshold; if all MinFrac values are NA (i.e.
#'   dataset contains blanks), NULL is returned
#' @examples
#'   file <- system.file('extdata','Search_Parameters.txt', package = "LUMA")
#'   search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
#'
#'   method = "mz"
#'   class(method) = method
#'   test <- trim_minfrac(Peak.list = Peaklist_Pos$From_CAMERA_with_MinFrac,
#'   search.par = search.par, object = method)
#'   nrow(Peaklist_Pos$From_CAMERA_with_MinFrac) -  nrow(test) #equals 6
#'
#'   method = "monoMass"
#'   class(method) = method
#'   test <- trim_minfrac(Peak.list =
#'   Peaklist_Pos$Combined_Isotopes_and_Adducts, search.par = search.par, object
#'   = method)
#'   nrow(Peaklist_Pos$Combined_Isotopes_and_Adducts) - nrow(test) #equals 4
trim_minfrac = function(object, Peak.list, search.par) {
  UseMethod("trim_minfrac", object)
}

#' @rdname trim_minfrac
#' @export
trim_minfrac.mz = function(object, Peak.list, search.par) {

  MF <- Peak.list[["MinFrac"]]

  if(all(!is.na(MF))) {   # Safeguard in case this function is called on blanks; if BLANK = TRUE in Script info, then MF will contain all NA values

    MF.cutoff <- as.numeric(search.par[1, "Minfrac"])
    Peak.list.trimmed <- Peak.list[MF >= MF.cutoff, ]
    return(Peak.list.trimmed)
  }

}

#' @rdname trim_minfrac
#' @export
trim_minfrac.monoMass = function(object, Peak.list, search.par) {

  MF <- Peak.list[["MinFrac"]]


  if(all(!is.na(MF))) {   # Safeguard in case this function is called on blanks; if BLANK = TRUE in Script info, then MF will contain all NA values
    AllMF <- strsplit(MF, split = ";")
    AllMF <- lapply(AllMF, function(x) as.numeric(x))  #Convert character values to numeric values
    MaxMF <- lapply(AllMF, function(x) max(x))
    # MeanMF <- lapply(AllMF, function(x) mean(x)) #calculates mean values for minfrac trimming; not used
    MF.cutoff <- as.numeric(search.par[1, "Minfrac"])
    # temp <- data.frame(MinFrac = MF, Max.cutoff = MaxMF>=MF.cutoff, Mean.cutoff = MeanMF>=MF.cutoff, CV =
    # Peak.list[,'%CV'], Corr.stat = Peak.list[,'Correlation.stat'])
    # temp[which(temp$Max.cutoff!=temp$Mean.cutoff),] #compares mean thresholding to max thresholding
    Peak.list.trimmed <- Peak.list[MaxMF >= MF.cutoff, ]
    return(Peak.list.trimmed)
  }
}


#' @title Trims by retention time
#'
#' @export
#' @description Removes components with retention times smaller than the user
#'   specified threshold. See \code{remove_void_volume} and
#'   \code{CullVoidVolume} for examples that use this function.
#' @param object used for method dispatch. Can be any object. See usage for
#'   details
#' @param Peak.list data frame. Must have column \code{Correlation.stat}.
#'   Should contain output columns from XCMS and CAMERA, and additional columns
#'   from \code{match_Annotation}, \code{calc_minfrac}, \code{calc_corrstat} and
#'   \code{plot_metgroup} functions.
#' @param void.rt numeric retention time cutoff value corresponding to the LC
#'   void volume
#' @param rt.list numeric vector containing retention times for all features or
#'   compounds
#' @return NULL
trim_rt <- function(object, Peak.list, void.rt, rt.list) {
  UseMethod("trim_rt", object)
}

#' @rdname trim_rt
#' @export
trim_rt.mz <- function(object, Peak.list, void.rt, rt.list = Peak.list["rt"][[1]]) {
  drops <- Peak.list[rt.list < void.rt, "EIC_ID"][[1]]  #Creates a vector of features with rt in the void volume
  Peak.list.trimmed <- Peak.list[which(!(Peak.list$EIC_ID) %in% drops),]
  return(Peak.list.trimmed)
}

#' @rdname trim_rt
#' @export
trim_rt.monoMass  <- function(object, Peak.list, void.rt, rt.list = Peak.list["rt"][[1]]) {
  drops <- Peak.list[rt.list < void.rt, "metabolite_group"]  #Creates a vector of metabolite groups that contain at least one feature with rt in the void volume
  length(which(!unlist(Peak.list[, "metabolite_group"]) %in% unlist(drops)))
  met.list <- Peak.list["metabolite_group"]
  Peak.list.trimmed <- Peak.list[which(!unlist(met.list) %in% unlist(drops)), ]
  return(Peak.list.trimmed)
}
