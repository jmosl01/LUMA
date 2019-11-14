#' @title Removes void volume from Peak.list
#'
#' @export
#' @description Removes features or compounds found in the void volume of the chromatographic run from the Peak.list
#' @param Peak.list data frame. Must have Correlation.stat column.  Should contain output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac, Calc.corr.stat and EIC.plotter functions.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param method which method to apply to trim by retention time.  See trim_rt for details
#' @param ... Arguments to pass to trim_rt
#' @return data frame Peak.list.trimmed original Peak.list without all metabolite groups containing at least one feature in the void volume
remove_void_volume = function(Peak.list, search.par, method,...) {
    void.rt <- as.numeric(search.par[1, "Voidrt"])
    class(method) <- method
    Peak.list.trimmed <- trim_rt(method,Peak.list,void.rt,...)
    return(Peak.list.trimmed)
}

#' @title Trims by CV
#'
#' @export
#' @description Removes metabolites with coefficient of variation greater than the user specified threshold, calculated from the QC samples
#' @param Peak.list data frame. Must have QC sample columns that contain the string 'Pooled_QC_'.  Should contain output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac, CAMERA.parser, Calc.corr.stat and Combine.phenodata base functions.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @return data frame Peak.list.trimmed original Peak.list without all metabolite groups with coefficient of variation greater than user specified threshold; if dataset contains blanks, data.frame with all NA values is returned
#' @importFrom stats sd
trim_cv = function(Peak.list, search.par) {
    res <- lapply(colnames(Peak.list), function(ch) grep("Pooled_QC_", ch))
    QC.list <- Peak.list[sapply(res, function(x) length(x) > 0)]
    QCsd <- apply((as.matrix(QC.list)), 1, sd)
    QCmean <- rowMeans(QC.list)
    RSD <- QCsd/QCmean
    Peak.list[, "%CV"] <- RSD
    CV.cutoff <- as.numeric(search.par[1, "CV"])
    Peak.list.trimmed <- Peak.list[RSD < CV.cutoff, ]
    return(Peak.list.trimmed)
}

#' @title Trims by MinFrac
#'
#' @export
#' @description Removes metabolites with MinFrac smaller than the user specified threshold. The maximum MinFrac value is chosen from all features within a metabolite group.
#' @param Peak.list data frame. Must have MinFrac column.  Should contain output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac, CAMERA.parser, Calc.corr.stat and Combine.phenodata base functions.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @return data frame Peak.list.trimmed original Peak.list containing all metabolite groups containing at least one feature that has MinFrac value greater than user specified threshold; if all MinFrac values are NA (i.e. dataset contains blanks), NULL is returned
trim_minfrac = function(Peak.list, search.par) {
    MF <- Peak.list[, "MinFrac"]
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
#' @description Removes components with retention times smaller than the user specified threshold.
#' @param object used for method dispatch. Can be any object. See usage for details
#' @param Peak.list data frame. Must have Correlation.stat column.  Should contain output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac, Calc.corr.stat and EIC.plotter functions.
#' @param void.rt numeric retention time cutoff value corresponding to the LC void volume
#' @param rt.list numeric vector containing retention times for all features or compounds
#' @return NULL
trim_rt <- function(object, Peak.list, void.rt, rt.list) {
  UseMethod("trim_rt", object)
}

#' @rdname trim_rt
#' @export
trim_rt.mz <- function(object, Peak.list, void.rt, rt.list = Peak.list["rt"]) {
  drops <- Peak.list[rt.list < void.rt, "EIC_ID"]  #Creates a vector of features with rt in the void volume
  Peak.list.trimmed <- Peak.list[which(!(Peak.list$EIC_ID) %in% drops),]
  return(Peak.list.trimmed)
}

#' @rdname trim_rt
#' @export
trim_rt.monoMass  <- function(object, Peak.list, void.rt, rt.list = Peak.list["rt"]) {
  drops <- Peak.list[rt.list < void.rt, "metabolite_group"]  #Creates a vector of metabolite groups that contain at least one feature with rt in the void volume
  length(which(!unlist(Peak.list[, "metabolite_group"]) %in% unlist(drops)))
  met.list <- Peak.list["metabolite_group"]
  Peak.list.trimmed <- Peak.list[which(!unlist(met.list) %in% unlist(drops)), ]
  return(Peak.list.trimmed)
}
