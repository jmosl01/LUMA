#' @title Wraps xcms
#'
#' @export
#' @description Run XCMS with user defined input parameters and return xcms objects
#' @param mzdatafiles a character vector of data files with full path names
#' @param XCMS.par a single-row data frame with 13 variables containing XCMS parameters. The column names must be c('Peakwidth1','Peakwidth2','ppm','noise','snthresh','mzdiff','prefilter1','prefilter2','center','gapInit','bw','mzwid','minfrac')
#' @param file.base a single-row data frame with 13 variables containing XCMS parameters. The column names must be c('Peakwidth1','Peakwidth2','ppm','noise','snthresh','mzdiff','prefilter1','prefilter2','center','gapInit','bw','mzwid','minfrac')
#' @return two XCMS objects xset and xset4 without and with retention time alignment, peak grouping, and imputing missing values
#' @import xcms
#' @importFrom BiocParallel SnowParam snowWorkers
wrap_xcms = function(mzdatafiles, XCMS.par, file.base) {
  #added me >
  #mzdatafiles <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE) #This will cause xcms to run on EVERY file in mzML directory
  # file.base <- gen_filebase(mzdatafiles, BLANK, ion.id, ion.mode) #Dont do this
  #added me <
    xset <- xcmsSet(files = mzdatafiles, method = "centWave", peakwidth = c(XCMS.par$Peakwidth1, XCMS.par$Peakwidth2),
        ppm = XCMS.par$ppm, noise = XCMS.par$noise, snthresh = XCMS.par$snthresh, mzdiff = XCMS.par$mzdiff, prefilter = c(XCMS.par$prefilter1,
            XCMS.par$prefilter2), mzCenterFun = "wMean", integrate = 1, fitgauss = FALSE, verbose.columns = FALSE,
        BPPARAM = SnowParam(workers = snowWorkers(), type = "SOCK", stop.on.error = TRUE, progressbar = TRUE))
    pdf(file = paste(file.base, "RTDev Plot.pdf", sep = "_"))
    xset2 <- retcor(xset, method = "obiwarp", plottype = "deviation", distFunc = "cor_opt", profStep = 1, center = XCMS.par$center,
        response = 1, gapInit = XCMS.par$gapInit, gapExtend = 2.7, factorDiag = 2, factorGap = 1, localAlignment = 0)
    dev.off()
    xset3 <- group(xset2, method = "density", bw = XCMS.par$bw, mzwid = XCMS.par$mzwid, minfrac = XCMS.par$minfrac,
        minsamp = 1, max = 50)

    xset4 <- fillPeaks(xset3, BPPARAM = SnowParam(workers = snowWorkers(), type = "SOCK", stop.on.error = TRUE,
        progressbar = TRUE))
    return(list(xset, xset4))
}

#' @title Wraps CAMERA
#'
#' @export
#' @description Run CAMERA with user defined input parameters and return xsAnnotate objects
#' @param xcms.obj an xcms object that has had peak picking, retention time alignment, peak grouping, and imputing missing values performed
#' @param CAMERA.par a single-row data frame with 9 variables containing CAMERA parameters. The column names must be c('perfwhm','sigma','minfrac','mzabs','maxiso','corval_eic','corval_exp','pval','mzabs.1')
#' @param ion.mode a character string defining the ionization mode.  Must be either 'positive' or 'negative'
#' @return two grouped xsannotate objects mz1setpos and anposGa without and with annotated isotopes and ion adducts and fragments
#' @import CAMERA
wrap_camera = function(xcms.obj, CAMERA.par, ion.mode) {
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
    CAMERA.ion.mode <- tolower(ion.mode)
    #me <
    mz1setpos <- xsAnnotate(xs = xcms.obj, sample = NA)
    mz1setpos <- groupFWHM(object = mz1setpos, perfwhm = best.perfwhm, sigma = best.sigma, intval = "into")
    mz1setposIso <- findIsotopes(mz1setpos, maxcharge = 2, maxiso = best.maxiso, ppm = 3, mzabs = best.mzabs.iso,
        intval = "into", minfrac = best.minfrac, filter = TRUE)
    mz1setposIso.new <- groupCorr(object = mz1setposIso, cor_eic_th = best.corval_eic, cor_exp_th = best.corval_exp,
        graphMethod = graph_method, pval = best.pval, calcIso = TRUE, calcCiS = TRUE, calcCaS = TRUE, psg_list = NULL,
        xraw = NULL)
    anposGa <- findAdducts(mz1setposIso.new, polarity = CAMERA.ion.mode, ppm = 3, mzabs = best.mzabs.add, multiplier = 2,
        rules = rules)
    peakGa <- getPeaklist(object = anposGa)
    EIC_ID <- row.names(peakGa)
    peak_data <- cbind(EIC_ID, peakGa)
    return(list(mz1setpos, anposGa))
}

#' @title Writes xlsx table
#'
#' @export
#' @description Write xlsx table output from validate_metgroup.  Essentially a wrapper for xlsx::saveWorkbook
#' @param file.base character return from gen_filebase function
#' @param validate.sheets list of sheets to write to xlsx file.  Currently must be of length 2
#' @param myname name to append to file.base to create file name for xlsx
#' @param mysheets character vector to name sheets in xlsx file. Currently must be of length 2
#' @return class jobjRef object
#' @importFrom openxlsx createWorkbook addWorksheet writeDataTable saveWorkbook
write_xlsx <- function(validate.sheets,file.base,myname,mysheets) {
  if(!is.list(validate.sheets) || length(validate.sheets) != 2)
    stop("validate.sheets must be a list with exactly two objects.", call. = FALSE)
  if(missing(mysheets))
    mysheets <- c("clear","muddy")
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
