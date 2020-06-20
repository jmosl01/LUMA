#' @title Wraps xcms
#'
#' @export
#' @description Run a simple XCMS workflow with user defined input parameters
#'   and return xcms objects
#' @param mzdatafiles a character vector of data files with full path names
#' @param XCMS.par a single-row data frame with 13 variables containing
#'   parameters for \code{xcms}. Must include the columns \code{c("Peakwidth1",
#'   "Peakwidth2","ppm","noise","snthresh","mzdiff","prefilter1","prefilter2",
#'   "center","gapInit","bw","mzwid","minfrac")}.
#' @param file.base character return from \code{gen_filebase}.
#' @return two \code{xcmsSet} objects \code{xset,xset4} without and with
#'   retention time alignment, peak grouping, and imputed missing values,
#'   respectively.
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'
#'     file <- system.file("extdata","Sample_Data.csv", package =
#'     "lcmsfishdata") # is case sensitive on Linux
#'     sample_data <- read.table(file, header = TRUE, sep = ",")
#'     mzdatafiles <- sample_data$CT.ID
#'
#'     file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode
#'     = "Positive", ion.id = c("Pos","Neg"))
#'     file2 <- system.file("extdata","Best_XCMS_parameters_positive.csv",
#'     package = "lcmsfishdata")  # is case sensitive on Linux
#'     XCMS.par <- read.table(file2, header = TRUE, sep = ",")
#'
#'     \dontrun{
#'
#'     #Runs XCMS This requires access to raw datafiles and won't work with
#'     #lcmsfishdata. Better to use your own data here.
#'     test <- wrap_xcms(mzdatafiles = mzdatafiles, XCMS.par = XCMS.par, file.base = file.base)
#'     }
#'
#' }
wrap_xcms = function(mzdatafiles, XCMS.par, file.base) {
  #added me >
  #mzdatafiles <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE) #This will cause xcms to run on EVERY file in mzML directory
  # file.base <- gen_filebase(mzdatafiles, BLANK, ion.id, IonMode) #Dont do this
  #added me <

  if(!requireNamespace("xcms", quietly = TRUE)) {
    stop("You must install xcms to use xcms_wrap! See installation instructions at:
         \n\nhttps://www.bioconductor.org/packages/release/bioc/html/xcms.html\n\n\n")
  } else {

    if(requireNamespace("BiocParallel", quietly = TRUE)) {

      xset <- xcms::xcmsSet(files = mzdatafiles, method = "centWave", peakwidth = c(XCMS.par$Peakwidth1, XCMS.par$Peakwidth2),
                      ppm = XCMS.par$ppm, noise = XCMS.par$noise, snthresh = XCMS.par$snthresh, mzdiff = XCMS.par$mzdiff, prefilter = c(XCMS.par$prefilter1,
                                                                                                                                        XCMS.par$prefilter2), mzCenterFun = "wMean", integrate = 1, fitgauss = FALSE, verbose.columns = FALSE,
                      BPPARAM = BiocParallel::SnowParam(workers = BiocParallel::snowWorkers(), type = "SOCK", stop.on.error = TRUE, progressbar = FALSE))

    } else {
      warning("Running xcms_wrap without BiocParallel takes a very long time! See installation instructions at:
         \n\nhttps://bioconductor.org/packages/release/bioc/html/BiocParallel.html\n\n\n")

    }

    pdf(file = paste(file.base, "RTDev Plot.pdf", sep = "_"))
    xset2 <- xcms::retcor(xset, method = "obiwarp", plottype = "deviation", distFunc = "cor_opt", profStep = 1, center = XCMS.par$center,
                    response = 1, gapInit = XCMS.par$gapInit, gapExtend = 2.7, factorDiag = 2, factorGap = 1, localAlignment = 0)
    dev.off()
    xset3 <- xcms::group(xset2, method = "density", bw = XCMS.par$bw, mzwid = XCMS.par$mzwid, minfrac = XCMS.par$minfrac,
                   minsamp = 1, max = 50)

    xset4 <- xcms::fillPeaks(xset3, BPPARAM = BiocParallel::SnowParam(workers = BiocParallel::snowWorkers(), type = "SOCK", stop.on.error = TRUE,
                                                  progressbar = FALSE))
    return(list(xset, xset4))
  }
}

#' @title Wraps CAMERA
#'
#' @export
#' @description Run a simple CAMERA workflow with user defined input parameters
#'   and return xsAnnotate objects
#' @param xcms.obj an \code{xcmsSet} object that has had peak picking, retention time
#'   alignment, peak grouping, and imputing missing values performed.
#' @param CAMERA.par a single-row data frame with 9 variables containing CAMERA
#'   parameters. Must include columns \code{c("perfwhm","sigma","minfrac","mzabs",
#'   "maxiso","corval_eic","corval_exp","pval","mzabs.1")}.
#' @param IonMode a character string defining the ionization mode.  Must be one
#'   of \code{c("Positive","Negative")}.
#' @return two grouped \code{xsannotate} objects \code{mz1setpos,anposGa}
#'   without and with annotated isotopes and ion adducts and fragments,
#'   respectively.
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr findIsotopes findAdducts getPeaklist
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#' file <- system.file("extdata","XCMS_objects_pos.Rdata", package =
#' "lcmsfishdata") # is case sensitive on Linux
#' load(file, envir = environment())
#' file2 <- system.file("extdata","Best_CAMERA_parameters_positive.csv", package
#' = "lcmsfishdata")  # is case sensitive on Linux
#' CAMERA.par <- read.table(file2, header = TRUE, sep = ",")
#'
#' \dontrun{
#'
#' #Runs CAMERA. This requires access to raw datafiles and won't work with
#' #lcmsfishdata. Better to use your own data here.
#' test <- wrap_camera(xcms.obj = xset4, CAMERA.par = CAMERA.par, IonMode = "Positive")
#' }
#'
#' }
wrap_camera = function(xcms.obj, CAMERA.par, IonMode) {
    best.perfwhm <- CAMERA.par$perfwhm
    best.sigma <- CAMERA.par$sigma
    best.mzabs.iso <- CAMERA.par$mzabs
    best.minfrac <- CAMERA.par$minfrac
    best.maxiso <- CAMERA.par$maxiso
    best.corval_exp <- CAMERA.par$corval_exp
    best.corval_eic <- CAMERA.par$corval_eic
    best.pval <- CAMERA.par$pval
    best.mzabs.add <- CAMERA.par$mzabs.1
    #me >
    graph_method <- "lpc"
    CAMERA_IonMode <- tolower(IonMode)
    #me <


    cat("\n\n")
    mz1setpos <- xsAnnotate(xs = xcms.obj, sample = NA)
    cat("\n\n")
    mz1setpos <- groupFWHM(object = mz1setpos, perfwhm = best.perfwhm, sigma = best.sigma, intval = "into")
    cat("\n\n")
    mz1setposIso <- findIsotopes(mz1setpos, maxcharge = 2, maxiso = best.maxiso, ppm = 3, mzabs = best.mzabs.iso,
        intval = "into", minfrac = best.minfrac, filter = TRUE)
    cat("\n\n")
    mz1setposIso.new <- groupCorr(object = mz1setposIso, cor_eic_th = best.corval_eic, cor_exp_th = best.corval_exp,
        graphMethod = graph_method, pval = best.pval, calcIso = TRUE, calcCiS = TRUE, calcCaS = TRUE, psg_list = NULL,
        xraw = NULL)
    cat("\n\n")
    anposGa <- findAdducts(mz1setposIso.new, polarity = CAMERA_IonMode, ppm = 3, mzabs = best.mzabs.add, multiplier = 2,
        rules = rules)
    cat("\n\n")
    peakGa <- getPeaklist(object = anposGa)
    cat("\n\n")
    EIC_ID <- row.names(peakGa)
    peak_data <- cbind(EIC_ID, peakGa)
    return(list(mz1setpos, anposGa))
}

#' @title Writes xlsx table
#'
#' @export
#' @description Write xlsx table output from \code{plot_metgroup}.  Essentially
#'   a wrapper for \code{openxlsx::saveWorkbook}
#' @param file.base character return from \code{gen_filebase}.
#' @param validate.sheets list of sheets to write to xlsx file.  Currently must
#'   be of length 2
#' @param myname character string to append to file.base to name xlsx file
#' @param mysheets character vector to name sheets in xlsx file. Currently must
#'   be of length 2. Default is \code{c("clear","muddy")}
#' @return \code{Workbook} from the \code{openxlsx} package
#' @importFrom openxlsx createWorkbook addWorksheet writeDataTable saveWorkbook
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'
#'   file <- system.file("extdata","Sample_Class.txt", package = "lcmsfishdata") # is case sensitive on Linux
#'   Sample.df <- read.table(file, header = TRUE, sep = "\t")
#'   file2 <- system.file("extdata","CAMERA_objects_Pos.Rdata", package = "lcmsfishdata") # is case sensitive on Linux
#'   load(file2, envir = environment())
#'   Peak.list <- lcmsfishdata::Peaklist_Pos[["input_parsed"]]
#'   file3 <- system.file("extdata","Sample_Data.csv", package = "lcmsfishdata") # is case sensitive on Linux
#'   sample_data <- read.table(file3, header = TRUE, sep = ",")
#'   mzdatafiles <- sample_data$CT.ID
#'
#'   file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode
#'   = "Positive", ion.id = c("Pos","Neg"))
#'
#'   mylist <- plot_metgroup(CAMERA.obj = anposGa, Sample.df = Sample.df,
#'   Peak.list = Peak.list, center = 2, BLANK = FALSE, gen.plots = FALSE,
#'   IonMode = "Positive", file.base = file.base, QC.id = "Pooled_QC", maxlabel
#'   = 10)
#'   class(mylist) ##is list
#'   length(mylist) ## with 2 elements
#'   validate.sheets <- mylist[[2]]
#'   myname <- "Validate_Metabolite_Groups"
#'   test <- write_xlsx(validate.sheets = validate.sheets, file.base = file.base, myname = myname)
#'   class(test) #returns Workbook object from the openxlsx package
#'   file.remove(paste(file.base, paste(myname,".xlsx", sep = ""), sep = "_"))
#'   }
write_xlsx <- function(validate.sheets,file.base,myname,mysheets) {

  #Set Default values
  if(missing(mysheets))
    mysheets <- c("clear","muddy")


  if(!is.list(validate.sheets) || length(validate.sheets) != 2)
    stop("validate.sheets must be a list with exactly two objects.", call. = FALSE)
  if(!is.vector(mysheets) || length(mysheets) != 2)
    stop("mysheets must be a vector of length 2!", call. = FALSE)

  wb = createWorkbook()
  S1 = addWorksheet(wb, mysheets[1])
  S2 = addWorksheet(wb, mysheets[2])
  df1 <- data.frame(validate.sheets$clear)
  df2 <- data.frame(validate.sheets$muddy)
  writeDataTable(wb, sheet = S1,x = df1)
  writeDataTable(wb, sheet = S2,x = df2)
  saveWorkbook(wb,
               file = paste(file.base,
                                paste(myname,".xlsx", sep = ""), sep = "_"),
               overwrite = TRUE)
  return(wb)
}
