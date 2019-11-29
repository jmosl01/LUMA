#' @title Input Files Dialog Box
#'
#' @export
#' @description Creates dialog box for user to select input files when initiating LUMA workflow
#' @param WorkingDir What is the working directory for the LUMA workflow
#' @param multiple should dialog box with multiple fields be used; caution this is buggy and doesn't always work on all platforms; default if FALSE
#' @importFrom svDialogs ok_cancel_box dlg_form dlg_open dlg_dir dlg_list
#' @return named list
InputFiles_dlg <- function(WorkingDir,multiple = FALSE) {

  #Set initial values
  XCMS.par <- NULL

  #Set default values
  if(missing(multiple)) multiple = FALSE

  #Start Dialog Box; if multiple, then use dlg_form; otherwise, use individual dlg functions for each query
  if(multiple) {

    #Temporarirly change to WorkingDir
    mydir <- getwd()

    setwd(WorkingDir)

    #using dlgForm
    form <- list(
      "SampleClass:SFL" = "Select the Sample Class file or copy and paste into console",
      "SampleData:SFL" = "Select the Sample Data file or copy and paste into console",
      "AnnotatedLibrary:SFL" = "Select the Annotated Library file or copy and paste into console",
      "SearchPar:SFL" = "Select the Search Parameters file or copy and paste into console",
      "XCMSObj:CHK" = FALSE,
      "CAMERAObj:CHK" = FALSE,
      "Adducts:SFL" = "Select the Adducts file or copy and paste into console"
    )
    myresults <- dlg_form(form, title = "Input Files")$res

    XCMSObj <- myresults$XCMSObj
    CAMERAObj <- myresults$CAMERAObj

  } else {

    #Using individual dlg functions for each query
    SampleClass <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the Sample Class file or copy and paste into console")$res
    SampleData <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the Sample Data file or copy and paste into console")$res
    AnnotatedLibrary <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the Annotated Library file or copy and paste into console")$res
    SearchPar <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the Search Parameters file or copy and paste into console")$res
    XCMSObj <- ok_cancel_box(message = "If you are using existing saved XCMS objects, click OK. Otherwise click Cancel.")
    CAMERAObj <- ok_cancel_box(message = "If you are using existing saved CAMERA objects, click OK. Otherwise click Cancel.")
    Adducts <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the Adducts file or copy and paste into console")$res

  }

  if(XCMSObj) {

    XCMSObj <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the saved XCMS Objects file or copy and paste into console")$res
    form <- list("XCMSCenter:NUM" = "Enter the value for the XCMS center parameter; must be a positive integer")

    xcsmresult <- dlg_form(form, title = "XCMS Parameter Input")$res
    XCMS.par <- data.frame(center = as.numeric(xcsmresult$XCMSCenter))

  } else {

    XCMSObj <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the XCMS parameters file or copy and paste into console")$res

  }

  if(CAMERAObj) {

    CAMERAObj <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the saved CAMERA Objects file or copy and paste into console")$res

  } else {

    CAMERAObj <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the CAMERA parameters file or copy and paste into console")$res

  }


  if(multiple) {

    myresults <- lapply(myresults, function(x) paste(WorkingDir, x, sep = "/"))

    myresults$XCMSObj <- XCMSObj

    myresults$CAMERAObj <- CAMERAObj

    myresults$XCMS.par <- XCMS.par

    setwd(mydir)

  } else {

    myresults <- list(SampleClass = SampleClass,
                      SampleData = SampleData,
                      AnnotatedLibrary = AnnotatedLibrary,
                      SearchPar = SearchPar,
                      XCMSObj = XCMSObj,
                      CAMERAObj = CAMERAObj,
                      Adducts = Adducts,
                      XCMS.par = XCMS.par)

  }


  return(myresults)

}


#' @title Script Info Dialog Box
#'
#' @export
#' @description Creates dialog box for user to set basic workflow (used to be called script) information
#' @param multiple should dialog box with multiple fields be used; default if FALSE
#' @return named list
ScriptInfo_dlg <- function(multiple = FALSE) {

  #Set default values
  if(missing(multiple)) multiple = FALSE

  #Start Dialog Box; if multiple, then use dlg_form; otherwise, use individual dlg functions for each query

  if(multiple) {

    #Using dlgForm
    form <- list(
      "WorkingDir:DIR" = "Select your input directory or copy and paste into console",
      "Are you using the LUMA recommended data directory:CHK" = FALSE,
      "BLANK:CHK" = FALSE,
      "IonMode:CB" = c("Positive","Negative")
    )
    myresults <- dlg_form(form, title = "Script Info")$res

    while(length(myresults) == 0) {
      myresults <- dlg_form(form, title = "Script Info")$res
    }

    #Convert escaped backslashes to forward slashes in directory
    myresults$WorkingDir <- gsub("\\\\", "/", myresults$WorkingDir)

    #Set names so I can access access elements of dialog results by name
    names(myresults)[1:2] <- c("WorkingDir","DataDir")

  } else {

    #Using individual dlg functions for each query
    WorkingDir <- dlg_dir(default = getwd(), title = "Select your input directory or copy and paste into console")$res
    DataDir <- ok_cancel_box(message = "If you are using the LUMA recommended data directory, click OK. Otherwise click Cancel.")
    BLANK <- ok_cancel_box(message = "If blank samples are being processed, click OK. Otherwise click Cancel.")
    IonMode <- dlg_list(choices = c("Positive","Negative"))$res

    myresults <- list(WorkingDir = WorkingDir, DataDir = DataDir, BLANK = BLANK, IonMode = IonMode)



  }

  if(myresults[[2]]) {
    DataDir <- dlg_dir(default = getwd(), title = "Data Directory")$res
    mylist <- list.dirs(DataDir)
    if(length(mylist) != 4) {

      stop("Your Data Directory must only contain three subdirectories called 'Samples', 'Blanks', and 'PooledQCs'")

    } else {

      mylist <- mylist[-1]
      SamplesDir <- mylist[grep("Samples", mylist)]
      BlanksDir <- mylist[grep("Blanks", mylist)]
      PooledQCsDir <- mylist[grep("PooledQCs", mylist)]
      mydir <- list(SamplesDir = SamplesDir, BlanksDir = BlanksDir, PooledQCsDir = PooledQCsDir)
    }

  }  else {

    form <- list(
      "SampleDir:DIR" = "Select your Samples directory (or copy and paste into console)",
      "BlanksDir:DIR" = "Select your Blanks directory (or copy and paste into console)",
      "PooledQCsDir:DIR" = "Select your Pooled QCs directory (or copy and paste into console)"
    )
    mydir <- dlg_form(form, title = "Data Directories")$res

  }

  myresults$DataDir <- mydir

  return(myresults)

}

#' @title CAMERA Files Dialog Box
#'
#' @export
#' @description Creates dialog box for user to select saved CAMERA objects when running LUMA workflow
#' @param multiple should dialog box with multiple fields be used; caution this is buggy and doesn't always work on all platforms; default if FALSE
#' @importFrom svDialogs ok_cancel_box dlg_form dlg_open dlg_dir dlg_list
#' @return named list
CAMERAFiles_dlg <- function(multiple = FALSE) {

  #Set default values
  if(missing(multiple)) multiple = FALSE


  #Start Dialog Box; if multiple, then use dlg_form; otherwise, use individual dlg functions for each query
  if(multiple) {


    #using dlgForm
    form <- list(
      "WorkingDir:DIR" = "Select your working directory or copy and paste into console",
      "CAMERA_Pos_Obj:CHK" = FALSE,
      "CAMERA_Neg_Obj:CHK" = FALSE
    )
    myresults <- dlg_form(form, title = "Input Files")$res

    WorkingDir <- myresults$WorkingDir
    CAMERAPosObj <- myresults$CAMERA_Pos_Obj
    CAMERANegObj <- myresults$CAMERA_Neg_Obj

  } else {

    #Using individual dlg functions for each query
    WorkingDir <- dlg_dir(default = getwd(), title = "Select your working directory or copy and paste into console")$res
    CAMERAPosObj <- ok_cancel_box(message = "If you are using existing saved CAMERA objects in positive mode, click OK. Otherwise click Cancel.")
    CAMERANegObj <- ok_cancel_box(message = "If you are using existing saved CAMERA objects in negative mode, click OK. Otherwise click Cancel.")

  }

  if(CAMERAPosObj) {

    CAMERAPosObj <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the saved CAMERA Objects file for positive mode or copy and paste into console")$res

  } else {

    stop("You cannot combine ion modes with plotting if you don't have saved CAMERA objects.")

  }

  if(CAMERANegObj) {

    CAMERANegObj <- dlg_open(default = paste(WorkingDir, "*.*", sep = "/" ), title = "Select the saved CAMERA Objects file for negative mode or copy and paste into console")$res

  } else {

    stop("You cannot combine ion modes with plotting if you don't have saved CAMERA objects.")

  }


  if(multiple) {

    myresults <- lapply(myresults, function(x) paste(WorkingDir, x, sep = "/"))

    myresults$CAMERAPosObj <- CAMERAPosObj

    myresults$CAMERANegObj <- CAMERANegObj


  } else {

    myresults <- list(CAMERAPosObj = CAMERAPosObj, CAMERANegObj = CAMERANegObj)

  }


  return(myresults)

}


