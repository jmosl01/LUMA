context("test-wrapper-functions")

test_that("xlsx file is correctly written", {

  if(require(lcmsfishdata, quietly = TRUE)) {

    file <- system.file("extdata/Sample_Class.txt", package = "lcmsfishdata")
    Sample.df <- read.table(file, header = TRUE, sep = "\t")
    file2 <- system.file("extdata/CAMERA_objects_Pos.Rdata", package = "lcmsfishdata")
    load(file2, envir = environment())
    Peak.list <- lcmsfishdata::Peaklist_Pos[["input_parsed"]]
    file3 <- system.file("extdata/Sample_Data.csv", package = "lcmsfishdata")
    sample_data <- read.table(file3, header = TRUE, sep = ",")
    mzdatafiles <- sample_data$CT.ID

    file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode
                              = "Positive", ion.id = c("Pos","Neg"))

    mylist <- plot_metgroup(CAMERA.obj = anposGa, Sample.df = Sample.df,
                            Peak.list = Peak.list, center = 2, BLANK = FALSE, gen.plots = FALSE,
                            IonMode = "Positive", file.base = file.base, QC.id = "Pooled_QC", maxlabel
                            = 10)
    expect_equal("list",class(mylist)) ##is list
    expect_equal(2,length(mylist)) ## with 2 elements
    validate.sheets <- mylist[[2]]
    myname <- "Validate_Metabolite_Groups"
    test <- write_xlsx(validate.sheets = validate.sheets, file.base = file.base, myname = myname)
    expect_equal("Workbook",class(test)[1]) #returns Workbook object
    expect_equal("openxlsx", attr(class(test), "package")) #from openxlsx package
    expect_true(file.remove(paste(file.base, paste(myname,".xlsx", sep = ""), sep = "_"))) #File was written to working directory
  }


})

test_that("CAMERA run is setup correctly", {

  if(require(lcmsfishdata, quietly = TRUE)) {
    file <- system.file("extdata/XCMS_objects_pos.Rdata", package = "lcmsfishdata")
    load(file, envir = environment())
    file2 <- system.file("extdata/Best_CAMERA_parameters_positive.csv", package = "lcmsfishdata")
    CAMERA.par <- read.table(file2, header = TRUE, sep = ",")

    expect_true(exists("xset", envir = environment())) #Loads raw xcms file
    expect_true(exists("xset4", envir = environment())) #Loads corrected xcms file
    test <- c("perfwhm","sigma","minfrac","mzabs","maxiso","corval_eic","corval_exp","pval","mzabs.1")
    expect_true(all(test %in% names(CAMERA.par))) #CAMERA parameters file has all the right column headers
    #Runs CAMERA. This requires access to raw datafiles and won't work with lcmsfishdata. Better to use your own data here.
    expect_error(wrap_camera(xcms.obj = xset4, CAMERA.par = CAMERA.par, IonMode = "Positive"))

  }


})

test_that("XCMS run is setup correctly", {

  if(require(lcmsfishdata, quietly = TRUE)) {

    file <- system.file("extdata/Sample_Data.csv", package = "lcmsfishdata")
    sample_data <- read.table(file, header = TRUE, sep = ",")
    mzdatafiles <- sample_data$CT.ID

    file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode
                              = "Positive", ion.id = c("Pos","Neg"))
    file2 <- system.file("extdata/Best_XCMS_parameters_positive.csv", package
                         ="lcmsfishdata")
    XCMS.par <- read.table(file2, header = TRUE, sep = ",")


    test <- c("Peakwidth1","Peakwidth2","ppm","noise","snthresh","mzdiff","prefilter1","prefilter2","center","gapInit","bw","mzwid","minfrac")

    expect_true(all(test %in% names(XCMS.par))) #XCMS parameters file has all the right column headers
    #Runs XCMS This requires access to raw datafiles and won't work with lcmsfishdata. Better to use your own data here.
    expect_error(wrap_xcms(mzdatafiles = mzdatafiles, XCMS.par = XCMS.par, file.base = file.base))

  }


})
