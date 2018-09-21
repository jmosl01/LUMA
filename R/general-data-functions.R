#' @title Combine phenotype data
#'
#' @export
#' @description Combine phenotype data for each metabolite group with summed intensity values
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Peak.list data frame. Must have Correlation.stat, metabolite_group, and mono_mass columns.  Should contain output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac, CAMERA.parser and EIC.plotter functions.
#' @param Summed.list data frame containing metabolite group as first column and the rest summed intensities for each sample
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param BLANK a logical indicating whether blanks are being evaluated
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @return data frame Peak.list.summed with the munged phenotype columns up front followed by QC and sample columns
#' @importFrom plyr ddply
#' @importFrom dplyr '%>%' mutate_if summarise bind_cols
#' @importFrom utils str txtProgressBar setTxtProgressBar
#' @importFrom stringr str_count
combine_phenodata = function(Sample.df, Peak.list, Summed.list, search.par, BLANK, ion.mode) {

  mylist <- .gen_res(ion.mode,search.par,Peak.list,Sample.df,BLANK)
  Peaklist_corstat <- mylist[[1]]
  res <- mylist[[2]]

  # Creates a data frame with only the pheno data columns combined for each metabolite group
  pheno.list <- Peaklist_corstat[sapply(res, function(x) length(x) < 1)]  #Extracts all of the pheno columns for combining by metabolite group
  pheno.list[,"metabolite_group"] <- Peaklist_corstat[,"metabolite_group"]
  pheno.list <- pheno.list %>% mutate_if(is.factor, as.character)
  attributes(pheno.list)
  str(pheno.list)
  mylist <- sort(unique(Peaklist_corstat$metabolite_group))
  new.pheno.list <- as.data.frame(mylist)
  # i = 11 ##for debugging purposes
  total = length(colnames(pheno.list))
  cat("Combining metadata for isotopes and adducts.\n\n\n")
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  mypaste = function(pheno.list, i) {
    summarise(pheno.list, X = paste0(pheno.list[, i], collapse = ";"))
  }
  for (i in 1:total) {
    if (is.character(pheno.list[, i])) {
      # works new.pheno.list[,colnames(pheno.list)[i]] <- ddply(pheno.list, .(metabolite_group), summarise,
      # X=paste0(MS.ID, collapse = ';'))[2]
      new.pheno.list[, colnames(pheno.list)[i]] <- ddply(pheno.list, ~metabolite_group, function(x) mypaste(x,
                                                                                                            i))[2]
    } else if (!is.na(.isWhole(pheno.list[, i]))) {
      new.pheno.list[, colnames(pheno.list)[i]] <- ddply(pheno.list, ~metabolite_group, function(x) mypaste(x,
                                                                                                            i))[2]
    }
    setTxtProgressBar(pb, i)
  }

  # Create original pheno list to be replaced once combined with the summed results
  Iso.only.list <- subset(pheno.list, pheno.list$Correlation.stat == 1)  #extracts only the prime feature for each metabolite group across all columns

  # Replace pheno columns with combined phenotype data
  new.pheno.list <- new.pheno.list[, -which(colnames(new.pheno.list) %in% "metabolite_group")]
  nm <- intersect(colnames(new.pheno.list), colnames(Iso.only.list))
  colnames(new.pheno.list) <- c("metabolite_group", colnames(new.pheno.list)[-1])

  new.pheno.list[["MS.ID"]][match(Iso.only.list$metabolite_group, new.pheno.list$metabolite_group)]
  Iso.only.list[nm] <- lapply(nm, function(x) new.pheno.list[[x]][match(Iso.only.list$metabolite_group, new.pheno.list$metabolite_group)])

  # Combine munged pheno columns with summed data
  df1 <- Iso.only.list[order(Iso.only.list$metabolite_group), ]
  df2 <- Summed.list[order(Summed.list$metabolite_group), ]
  Peak.list.summed <- bind_cols(df1, df2)
  colnames(Peak.list.summed)
  temp <- as.character(Peak.list.summed$MS.ID)
  str(temp)

  #Drop Singletons
  no.features <- str_count(temp, ";") + 1
  drop.singletons = !search.par[1, "keep.singletons"]
  if (drop.singletons & BLANK == FALSE) {
    drops <- grep(1, no.features, value = FALSE)
    temp <- Peak.list.summed[drops, ]
    Peak.list.summed <- Peak.list.summed[-drops, ]
  }

  #Combine monomolecular mass values into a single value
  myMonoMass <- Peak.list.summed$mono_mass
  allMonoMass <- strsplit(myMonoMass, split = ";")
  allMonoMass <- lapply(allMonoMass, function(x) as.numeric(x))
  meanMonoMass <- lapply(allMonoMass, function(x) mean(x))
  Peak.list.summed$mono_mass <- unlist(meanMonoMass)

  #Combine rt values into a single value
  myRT <- Peak.list.summed$rt
  allRT <- strsplit(myRT, split = ";")
  allRT <- lapply(allRT, function(x) as.numeric(x))
  meanRT <- lapply(allRT, function(x) mean(x))
  Peak.list.summed$meanRT <- unlist(meanRT)

  return(Peak.list.summed)
}

#' @title Generates LUMA data table from XCMS or CAMERA object
#'
#' @export
#' @description Generates LUMA data table from xsAnnotate object generated from CAMERA, with the option to convert retention times to min
#' @param object an xcmsSet or xsAnnotate object
#' @param convert.rt logical indicating whether to convert retention times to min
#' @return data frame of peak intensities across all samples with metadata columns from XCMS and CAMERA
LUMA_getPeaklist = function(object,convert.rt) {
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

##Generates a parsing vector
.gen_res <- function(ion.mode,search.par,Peak.list,Sample.df,BLANK) {
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
