##Initialize global variables
xset <- NULL
xset4 <- NULL
mz1setpos <- NULL
anposGa <- NULL

## internal constructor utility functions ----
.get_DataFiles = function(mzdatapath,ion.mode,BLANK,ion.id,blanks.dir) {
  ## Selects the datafiles to use for data processing
  mzdatafiles <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE)
  if(ion.mode == "Positive" && BLANK == TRUE){
    mzdatafiles <- subset(mzdatafiles, subset = grepl(ion.id[1], mzdatafiles, ignore.case = TRUE))
    mzdatafiles <- mzdatafiles[c(grep(blanks.dir,mzdatafiles, ignore.case = TRUE))]
  } else {
    if(ion.mode == "Negative" && BLANK == TRUE){
      mzdatafiles <- subset(mzdatafiles, subset = grepl(ion.id[2], mzdatafiles, ignore.case = TRUE))
      mzdatafiles <- mzdatafiles[c(grep(blanks.dir,mzdatafiles))]
    } else {
      if(ion.mode == "Positive" && BLANK == FALSE){
        mzdatafiles <- subset(mzdatafiles, subset = grepl(ion.id[1], mzdatafiles, ignore.case = TRUE))
        mzdatafiles <- mzdatafiles[-c(grep(blanks.dir,mzdatafiles))]
      } else {
        if(ion.mode == "Negative" && BLANK == FALSE){
          mzdatafiles <- subset(mzdatafiles, subset = grepl(ion.id[2], mzdatafiles, ignore.case = TRUE))
          mzdatafiles <- mzdatafiles[-c(grep(blanks.dir,mzdatafiles))]
        } else {
          stop("Ion mode must be Positive or Negative.\nBe sure to specify whether to analyze blanks by setting BLANK to a logical. \nSee LUMA vignette for more details.\n\n")
        }
      }

    }

  }
  return(mzdatafiles)
}

.set_PreProcessFileNames = function(ion.mode,BLANK) {
  if(ion.mode == "Positive" && BLANK == TRUE){
    XCMS.file <- "XCMS_objects_Blanks_Pos"
    CAMERA.file <- "CAMERA_objects_Blanks_Pos"
  } else {
    if(ion.mode == "Negative" && BLANK == TRUE){
      XCMS.file <- "XCMS_objects_Blanks_Neg"
      CAMERA.file <- "CAMERA_objects_Blanks_Neg"
    } else {
      if (ion.mode == "Positive") {
        XCMS.file <- "XCMS_objects_Pos"
        CAMERA.file <- "CAMERA_objects_Pos"
      } else {
        XCMS.file <- "XCMS_objects_Neg"
        CAMERA.file <- "CAMERA_objects_Neg"
      }
    }
  }
  return(list(XCMS.file,CAMERA.file))
}

.xcmsSanityCheck = function(XCMS.obj) {
  if(length(XCMS.obj@filled) == 0) {
    warning("LUMA works best on xcms data that has been filled.\n\n")
    } else {
      xset4 <<- xset4 <- XCMS.obj
      save(xset4, file = XCMS.file)
      return(XCMS.obj)
    }
  }

.CAMERASanityCheck = function(CAMERA.obj) {
  if(length(CAMERA.obj@annoGrp) == 0) {
    warning("LUMA works best on CAMERA data that has been annotated.\n\n")
  } else {
    anposGa <<- anposGa <- CAMERA.obj
    save(anposGa, file = CAMERA.file)
    return(CAMERA.obj)
  }
}

.get_rules = function(ion.mode, files) {
  ## Reads in the adduct rules list for CAMERA
  if(ion.mode == "Positive"){
    rules <- read.csv(file = files[1])
    CAMERA.ion.mode <- "positive"
  } else {
    if(ion.mode == "Negative"){
      rules <- read.csv(file = files[2])
      CAMERA.ion.mode <- "negative"
    }
  }
  return(rules)
}

.get_Peaklist = function(object,convert.rt) {
  if(missing(convert.rt))
    convert.rt = TRUE
  peakGa <- getPeaklist(object)
  EIC_ID<-row.names(peakGa)
  peak_data <- cbind(EIC_ID, peakGa)
  ## Converts retention times to min from sec in Peaklist -----
  if(convert.rt) {
    rt.list<-peak_data["rt"]
    rt.min<-apply((as.matrix(rt.list)),1, function(x) x/60)
    peak_data["rt"] <- rt.min
  }
  return(peak_data)
}

.PreProcess_Files = function(XCMS.file,CAMERA.file,mytable) {
  if(file.exists(XCMS.file)) {
    cat("Reading in XCMS files.\n\n")
    load(file = XCMS.file)
    if(file.exists(CAMERA.file)){
      cat("Reading in CAMERA files.\n\n")
      load(file = CAMERA.file)
      ## Converts retention times to min from sec in Peaklist -----
      peak_data <- .get_Peaklist(anposGa)
      write_tbl(mydf = peak_data,
                peak.db = peak_db,
                myname = "From CAMERA")
    } else {
      cat("Running CAMERA Only!")
      # Runs CAMERA on datafiles --------------------
      ## Code to run CAMERA on XCMS object that has been rt corrected, grouped, and peaks filled
      time.CAMERA <- system.time({
        myresults <- wrap_camera(xset4 = xset4,
                                 CAMERA.par = CAMERA.par,
                                 ion.mode = CAMERA.ion.mode)
        mz1setpos <<- mz1setpos <- myresults[[1]]
        anposGa <<- anposGa <- myresults[[2]]
        peak_data <- .get_Peaklist(anposGa)
      })
      cat(paste("PreProcessing with CAMERA took ",round(print(time.CAMERA[3]))," seconds of elapsed time.\n\n",sep = ""))
      # Section END

      # Saves XCMS and CAMERA objects for re-analysis and peaklist for data processing ----
        save(xset,xset4,file=XCMS.file)
        save(mz1setpos,anposGa,file=CAMERA.file)

      ## Converts retention times to min from sec in Peaklist -----
      peak_data <- .get_Peaklist(anposGa)
      write_tbl(mydf = peak_data,
                peak.db = peak_db,
                myname = "From CAMERA")
    }

  } else {
    if(!file.exists(XCMS.file) & !file.exists(CAMERA.file)) {
      cat("Running XCMS and CAMERA: Be Patient!")

      # Runs XCMS on datafiles --------------------------------------
      time.XCMS <- system.time({
        myresults <- wrap_xcms(mzdatafiles = mzdatafiles,
                               XCMS.par = XCMS.par)
        xset <<- xset <- myresults[[1]]
        xset4 <<- xset4 <- myresults[[2]]
      })
      cat(paste("PreProcessing with XCMS took ",round(print(time.XCMS[3]))," seconds of elapsed time.\n\n",sep = ""))
      ## Section End
      # Runs CAMERA on datafiles --------------------
      ## Code to run CAMERA on XCMS object that has been rt corrected, grouped, and peaks filled
      time.CAMERA <- system.time({
        myresults <- wrap_camera(xset4 = xset4,
                                 CAMERA.par = CAMERA.par,
                                 ion.mode = CAMERA.ion.mode)
        mz1setpos <<- mz1setpos <- myresults[[1]]
        anposGa <<- anposGa <- myresults[[2]]
        peak_data <- .get_Peaklist(anposGa)
      })
      cat(paste("PreProcessing with CAMERA took ",round(print(time.CAMERA[3]))," seconds of elapsed time.\n\n",sep = ""))
      # Section END

      # Saves XCMS and CAMERA objects for re-analysis and peaklist for data processing ----
        save(xset,xset4,file=XCMS.file)
        save(mz1setpos,anposGa,file=CAMERA.file)

      ## Converts retention times to min from sec in Peaklist -----
      peak_data <- .get_Peaklist(anposGa)

      write_tbl(mydf = peak_data,
                peak.db = peak_db,
                myname = mytable)
    } else {
      stop("You must run XCMS before CAMERA! \nPlease remove CAMERA objects from your script directory.\n\n")
    }
    # Section END
  }
  return(xset4)
}

## END
## Other module utility functions ----
.gen_res = function(ion.mode,search.par,Peak.list,Sample.df,BLANK) {
  if (ion.mode == "Positive") {
    cor.stat <- as.numeric(search.par[1, "Corr.stat.pos"])
  } else {
    if (ion.mode == "Negative") {
      cor.stat <- as.numeric(search.par[1, "Corr.stat.neg"])
    }
  }
  Peaklist_corstat <- Peak.list[which(Peak.list$Correlation.stat >= cor.stat), ]
  if (BLANK == FALSE) {
    sexes <- unique(paste(Sample.df$Sex, "_", sep = ""))
    #Flags all of the sample columns and the metabolite group data
    res <- lapply(colnames(Peaklist_corstat),
                  function(ch) unique(grep(paste("Pooled_QC_", paste(strsplit(sexes,"(?<=.[_]$)", perl = TRUE), collapse = "|"),
                                                 "metabolite_group", sep = "|"), ch)))
  } else {
    if (ion.mode == "Positive" && BLANK == TRUE) {
      #Flags all of the sample columns and the metabolite group data
      res <- lapply(colnames(Peaklist_corstat),
                    function(ch) unique(grep("_Pos|metabolite_group", ch, ignore.case = TRUE)))
    } else {
      if (ion.mode == "Negative" && BLANK == TRUE) {
        #Flags all of the sample columns and the metabolite group data
        res <- lapply(colnames(Peaklist_corstat),
                      function(ch) unique(grep("_Neg|metabolite_group", ch, ignore.case = TRUE)))
      }
    }
  }
  return(list(Peaklist_corstat,res))
}

##END
