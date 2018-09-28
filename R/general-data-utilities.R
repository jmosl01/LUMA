#Generates a parsing vector
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

## internal constructor utility functions
.get_DataFiles = function(mzdatapath,ion.mode,BLANK,ion.id,blanks,dir) {
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
          msg <- "Ion mode must be Positive or Negative.\nBe sure to specify whether to analyze blanks by setting BLANK to a logical. \nSee LUMA vignette for more details.\n\n"
          .LUMAmsg <<- c(.LUMAmsg,msg)
          cat(msg)
          stop()
        }
      }

    }

  }
  return(mzdatafiles)
}

.get_rules = function(ion.mode, files) {
  ## Reads in the adduct rules list for CAMERA
  if(ion.mode == "Positive"){
    rules <- read.csv(file = files[1])
    CAMERA.ion.mode <<- "positive"
  } else {
    if(ion.mode == "Negative"){
      rules <- read.csv(file = files[2])
      CAMERA.ion.mode <<- "negative"
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

