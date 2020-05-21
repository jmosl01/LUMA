#' @title Formats Peak.list for SIMCA
#'
#' @export
#' @description Formats Peak.list for import into SIMCA
#' @param Peak.list data frame containing combined ion mode peaklist with ion mode duplicates removed.  Alternatively can be retrieved from databases.  Default is NULL
#' @param Sample.df data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Sample.data data frame with phenotype data as columns and a row for each study sample.  First column must be a unique sample identifier with the header 'CT-ID'.  Phenotype columns may vary, but must include two columns called 'Plate Number' and 'Plate Position' for determining run order.
#'    'Plate Number' must be numeric and is equivalent to batch number.  'Plate Position' must be alphanumeric and corresponds to row(alpha) and column(numeric) positions on, e.g. a 96-well plate
#' @param tbl.id character Table name to draw from database. Default is NULL
#' @param QC.id character vector specifying identifier in filename designating a Pooled QC sample.  Only the first value will be used.  Default is 'Pooled_QC_'
#' @param ... Arguments to pass parameters to database functions
#' @return NULL testing
#' @importFrom stats setNames
#' @importFrom gtools mixedorder
#' @importFrom plyr rbind.fill
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'   file <- system.file('extdata/Sample_Class.txt', package = "lcmsfishdata")
#'   Sample.df <- read.table(file, header = TRUE, sep = "\t")
#'  file2 <- system.file('extdata/Sample_Data.csv', package = "lcmsfishdata")
#'   Sample.data <- read.table(file2, header = TRUE, sep = ",")
#'  Peak.list <- Peaklist_db$Peaklist_Normalized
#'   test <- format_simca(Peak.list = Peak.list, Sample.df = Sample.df, Sample.data = Sample.data)
#' }
format_simca = function(Peak.list = NULL, Sample.df, Sample.data, tbl.id = NULL, QC.id = "Pooled_QC_", ...) {



  ##-----------------------------------------------------------------------------------------
  ## Initialize Peaklist from database
  ##-----------------------------------------------------------------------------------------
  if (missing(Peak.list))
    Peak.list = NULL
  if (missing(tbl.id))
    tbl.id = NULL
  if (is.null(tbl.id) && is.null(Peak.list)) {
    stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
  }
  if (is.null(Peak.list)) {
    Peak.list <- read_tbl(tbl.id, ...)
  }

  ##-----------------------------------------------------------------------------------------
  ## Format metabolite data.
  ##-----------------------------------------------------------------------------------------

  # Extract peak intensity data from Peak.list using the sexes as a search string.
  sexes <- unique(paste(Sample.df$Sex, "_", sep = ""))
  samples <- vector(mode = "character", length = length(colnames(Peak.list)))
  for (i in 1:length(sexes)) {
    rows_loop <- grep(sexes[i], colnames(Peak.list))
    samples[rows_loop] <- sexes[i]
  }
  res <- samples %in% sexes
  sample.peaks <- Peak.list[, res]


  ## Creates a new column for grouping by class based on user input
  groups <- paste(Sample.df$Sex, Sample.df$Class, sep = ";")  ## Generate search string for all classes
  groups <- strsplit(groups, split = ";")
  names(groups) <- paste(Sample.df$Sex, Sample.df$Class, sep = "_")
  group <- vector(mode = "character", length = length(colnames(sample.peaks)))
  for (i in 1:length(groups)) {
    rows_loop <- intersect(grep(groups[[i]][1], colnames(sample.peaks)),
                           grep(groups[[i]][2], colnames(sample.peaks))
    )
    group[rows_loop] <- names(groups)[i]
  }
  group <- unlist(group)

  # modify sample data to include the user defined exposure class
  Sample.data <- Sample.data[order(Sample.data[,1]), ]
  Sample.data <- setNames(cbind.data.frame(Sample.data[, 1], group, Sample.data[-1]), c(colnames(Sample.data[1]),
                                                                                        "Exposure Class", colnames(Sample.data[-1])))
  sample.ID <- suppressWarnings(as.numeric(sub("\\D*(\\d{6}).*", "\\1", colnames(sample.peaks))))  #pulls out the 6-digit numeric sample codes into a vector for matching against the Sample data spreadsheet
  if (all(is.na(sample.ID))) {
    ## Alternatively uses a unique alphanumeric ID provided by user
    sample.ID <- sub("\\D*(\\d{6}).*", "\\1", colnames(sample.peaks))
  }
  sample.peaks <- cbind.data.frame(Peak.list[, "MS.ID"], sample.peaks)
  # Create a new ID for each metabolite and record how many features were combined into that metabolite for
  # metadata
  temp <- as.character(sample.peaks[, 1])

  no.features <- str_count(temp, ";") + 1
  MetID <- gsub("^(.*?);.*", "\\1", temp)
  Mettag <- rep("Unidentified", length = length(MetID), mode = "character")
  Mettag[which(grepl("Annotated", temp, fixed = TRUE))] <- "Annotated"

  my_ion <- Peak.list$Ion.Mode.1
  my_mass <- round(Peak.list$mono_mass, digits = 5)
  my_rt <- round(Peak.list$meanRT, digits = 2)
  Metbase <- paste(my_ion,my_mass,my_rt,sep = "_")

  MetID <- paste(Metbase, Mettag, sep = "_")
  # End new ID generation

  tsample.peaks <- setNames(data.frame(t(sample.peaks[, -1])), MetID)
  # add columns for sex, sample class, and all phenotype data for each sample
  tsample.peaks <- cbind.data.frame(X = sample.ID, tsample.peaks)
  colnames(tsample.peaks)[1] = colnames(Sample.data)[1]
  tsample.peaks <- merge(Sample.data, tsample.peaks, by = colnames(tsample.peaks)[1])

  mylist <- split(tsample.peaks, f = tsample.peaks$Plate.Number)
  new.list <- lapply(mylist, function(X) X[mixedorder(X[, "Plate.Position"]), ])
  tsample.peaks <- rbind.fill(new.list)

  # make a dataframe of metabolite intensity values for pooled QC samples then transpose so the colnames are the
  # same as metabolite names: Pos_###_Annotated
  Pooled.QC <- QC.id  ## Generate search string for all sexes
  for (i in 1:length(Pooled.QC)) {
    rows_loop <- grep(Pooled.QC[i], colnames(Peak.list))
    samples[rows_loop] <- Pooled.QC[i]
  }
  QC <- samples %in% Pooled.QC
  QC.peaks <- Peak.list[, QC]
  QC.peaks <- cbind.data.frame(Peak.list[, "MS.ID"], QC.peaks)
  tQC.peaks <- setNames(data.frame(t(QC.peaks[, -1])), MetID)
  # add columns for sex, sample class, and all phenotype data for Pooled QCs, and fill with NA
  dummy.data <- setNames(data.frame(matrix(ncol = length(Sample.data), nrow = nrow(tQC.peaks))), colnames(Sample.data))
  tQC.peaks <- cbind.data.frame(dummy.data, tQC.peaks)


  # make a dataframe of all metadata columns then transpose so the colnames are the same as the metabolite
  # names: Pos_###_Annotated
  log.list <- list()
  log.list[[1]] <- res
  log.list[[2]] <- QC
  meta.peaks <- Peak.list[, !Reduce("|", log.list)]
  temp <- str_count(meta.peaks$Duplicate_EIC, ";")
  Ion.duplicate <- vector(mode = "logical", length = length(temp))
  Ion.duplicate[which(temp > 0)] = TRUE  #add flag for the metabolites detected in both ion modes
  endo <- meta.peaks[, "Endogenous_flag"]
  meta.peaks$Endogenous_flag <- as.logical(meta.peaks$Endogenous_flag)  #convert endogenous flag to a logical
  meta.peaks <- cbind.data.frame(meta.peaks[, 1:11], no.features, meta.peaks[, 12:ncol(meta.peaks)], Ion.duplicate)
  tmeta.peaks <- setNames(data.frame(t(meta.peaks[, -1])), MetID)
  dummy.data <- setNames(data.frame(matrix(ncol = length(Sample.data), nrow = nrow(tmeta.peaks))), colnames(Sample.data))
  tmeta.peaks <- cbind.data.frame(dummy.data, tmeta.peaks)
  # pheno.peaks <- Peak.list[,-grep('[MF]_',colnames(Peak.list))] colnames(pheno.peaks)

  # combine the three dataframes, first sample data, then QC data, lastly metadata
  SIMCA.data <- rbind.data.frame(tsample.peaks, tQC.peaks, tmeta.peaks)
  SIMCA.data <- cbind.data.frame(X = rownames(SIMCA.data), SIMCA.data)
  SIMCA.data[, 1] <- gsub("^X", "\\1", SIMCA.data[, 1])
  SIMCA.data <- setNames(cbind.data.frame(SIMCA.data, row.names = NULL), c("", colnames(SIMCA.data)[-1]))



}

#' @title Formatting of metabolomic data for MetaboAnalystR.
#'
#' @export
#' @description This function initializes objects that will hold the metabolite
#'   data, formats peak intensity data into one of the formats acceptable by
#'   MetaboAnalystR, and sets the metabolite data object.
#' @param mSetObj NULL
#' @param Peak.list data frame containing combined ion mode peaklist with ion
#'   mode duplicates removed.
#' @param Sample.df data frame with class info as columns.  Must contain a
#'   separate row entry for each unique sex/class combination. Must contain the
#'   columns 'Sex','Class','n','Endogenous'.
#' @param Sample.data data frame with phenotype data as columns and a row for
#'   each study sample.  First column must be a unique sample identifier with
#'   the header 'CT-ID'.  Phenotype columns may vary, but must include two
#'   columns called 'Plate Number' and 'Plate Position' for determining run
#'   order.
#' @param tbl.id character Table name to draw from database. Default is NULL
#' @param ... Arguments to pass parameters to database functions
#' @return mSetObj
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#' file <- system.file("extdata/Sample_Class.txt", package = "LUMA")
#' Sample.df <- read.table(file, header = TRUE, sep = "\t")
#' file2 <- system.file("extdata/Sample_Data.csv", package = "LUMA")
#' Sample.data <- read.table(file2, header = TRUE, sep = ",")
#' Peak.list <- Peaklist_db$Peaklist_Normalized
#' new_mSetObj <- format_MetabolomicData(mSetObj = mSetObj, Peak.list = Peak.list, Sample.df = Sample.df, Sample.data = Sample.data)
#' new_mSetObj$dataSet$orig.cls
#' }
format_MetabolomicData <- function(mSetObj, Peak.list, Sample.df, Sample.data, tbl.id, ...)
{
  UseMethod("format_MetabolomicData", mSetObj)
}

#' @rdname format_MetabolomicData
#' @export
format_MetabolomicData.pktable <- function(mSetObj, Peak.list, Sample.df, Sample.data, tbl.id, ...)
{

  ##-----------------------------------------------------------------------------------------
  ## Initialize Peaklist from database
  ##-----------------------------------------------------------------------------------------
  if (missing(Peak.list))
    Peak.list = NULL
  if (missing(tbl.id))
    tbl.id = NULL
  if (is.null(tbl.id) && is.null(Peak.list)) {
    stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
  }
  if (is.null(Peak.list)) {
    Peak.list <- read_tbl(tbl.id, ...)
  }




  ##-----------------------------------------------------------------------------------------
  ## Format metabolite data.
  ##-----------------------------------------------------------------------------------------

  # Extract peak intensity data from Peak.list using the sexes as a search string.
  sexes  <-  unique(paste(Sample.df$Sex, "_",  sep = ""))    ## Generate search string for all sexes
  samples  <-  vector(mode = "character", length = length(colnames(Peak.list)))
  for (i  in  1:length(sexes))
  {
    rows_loop  <-  grep(sexes[i], colnames(Peak.list))
    samples[rows_loop]  <-  sexes[i]
  }
  res  <-  samples %in% sexes
  sample.peaks  <- Peak.list[, res]

  # Generate sample IDs
  sample.ID  <-  suppressWarnings(as.numeric(sub("\\D*(\\d{6}).*", "\\1", colnames(sample.peaks))))
  if (all(is.na(sample.ID)))
  {
    ## Alternatively uses a unique alphanumeric ID provided by user
    sample.ID <- sub("\\D*(\\d{6}).*", "\\1", colnames(sample.peaks))
  }

  ## Code to generate treatment levels (sometimes confusingly called sample classes)
  # Creates a new column for grouping by class based on user input
  groups <- paste(Sample.df$Sex, Sample.df$Class, sep = ";")  ## Generate search string for all classes
  groups <- strsplit(groups, split = ";")
  names(groups) <- paste(Sample.df$Sex, Sample.df$Class, sep = "_")
  group <- vector(mode = "character", length = length(colnames(sample.peaks)))
  for (i in 1:length(groups)) {
    rows_loop <- intersect(grep(groups[[i]][1], colnames(sample.peaks)), grep(groups[[i]][2], colnames(sample.peaks)))
    group[rows_loop] <- names(groups)[i]
  }
  group <- unlist(group)

  # Modify sample data to include the user defined exposure class
  Sample.data <- Sample.data[order(Sample.data[,1]), ]
  Sample.data <- setNames(cbind.data.frame(Sample.data[, 1], group, Sample.data[-1]), c(colnames(Sample.data[1]),
                                                                                        "Exposure Class", colnames(Sample.data[-1])))

  # Extract only the treatment level column.
  treatment.Level <- Sample.data[, 2]

  # Transpose the peak intensity data frame.
  tsample.peaks <- t(sample.peaks)

  # Build the peak intensity table for the MetaboAnalystR package.
  tsample.peaks <- cbind.data.frame(sample.ID, treatment.Level, tsample.peaks)

  ##-----------------------------------------------------------------------------------------
  ## Set the metabolomic data object.
  ##-----------------------------------------------------------------------------------------
  mSetObj$dataSet$cls.type <-  "disc"
  mSetObj$dataSet$format <-  "rowu"
  dat  <-  tsample.peaks

  #From MetaboAnalyst::Read.TextData code
  smpl.nms <-dat[,1];

  all.nms <- colnames(dat);

  facA.lbl <- all.nms[2];

  cls.lbl <- facA <- dat[,2]; # default assign facA to cls.lbl in order for one-factor analysis
  conc <- dat[,-c(1,2)];
  var.nms <- colnames(conc);

  #assign the dimension names
  rownames(conc) <- smpl.nms
  colnames(conc) <- var.nms

  mSetObj$dataSet$type.cls.lbl <- class(cls.lbl);

  mSetObj$dataSet$orig.cls  <-  mSetObj$dataSet$cls  <-  cls.lbl;
  mSetObj$dataSet$orig  <-  conc;

  return(mSetObj)
}


#' @rdname format_MetabolomicData
#' @export
format_MetabolomicData.mass_all <- function(mSetObj, Peak.list, Sample.df, Sample.data, tbl.id, ...)
{

}

#' @rdname format_MetabolomicData
#' @export
format_MetabolomicData.default <- function(mSetObj, Peak.list, Sample.df, Sample.data, tbl.id, ...)
{
  warning(paste("format_MetabolomicData does not know how to handle object of class ",
                class(mSetObj),
                "and can only be used on classes pktable and mass_all"))
}


#' @title Printing metadata as CSV files for MetaboAnalystR.
#'
#' @export
#' @description This function prints the metabolite and sample metadata to CSV files.
#' @param mSetObj NULL
#' @param Peak.list data frame containing combined ion mode peaklist with ion mode duplicates removed.
#' @param Sample.df data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Sample.data data frame with phenotype data as columns and a row for each study sample.  First column must be a unique sample identifier with the header 'CT-ID'.  Phenotype columns may vary, but must include two columns called 'Plate Number' and 'Plate Position' for determining run order.
#' @importFrom utils write.csv
#' @return mSetObj
output_MetaData <- function(mSetObj, Peak.list, Sample.df, Sample.data)
{
  UseMethod("output_MetaData", mSetObj)
}

#' @rdname output_MetaData
#' @export
output_MetaData.pktable <- function(mSetObj, Peak.list, Sample.df, Sample.data)
{
  ##-----------------------------------------------------------------------------------------
  ## Output metabolite and sample metadata.
  ##-----------------------------------------------------------------------------------------

  # Determine location of peak intensity data from Peak.list using the sexes
  # as a search string.
  sexes  <-  unique(paste(Sample.df$Sex, "_",  sep = ""))    ## Generate search string for all sexes
  samples  <-  vector(mode = "character", length = length(colnames(Peak.list)))
  for (i  in  1:length(sexes))
  {
    rows_loop  <-  grep(sexes[i], colnames(Peak.list))
    samples[rows_loop]  <-  sexes[i]
  }
  res  <-  samples %in% sexes
  sample.peaks <- Peak.list[, res]

  # Apply the logical NOT operator "!" to "res" to identify the columns with the metadata.
  meta_res <- !res

  # Extract the metabolite metadata
  metabolite_metadata  <- Peak.list[, meta_res]

  # Write metabolite metadata to a CSV file
  write.csv(metabolite_metadata, "Metabolite_Metadata.csv")

  ## Code to generate the exposure class
  # Creates a new column for grouping by class based on user input
  groups <- paste(Sample.df$Sex, Sample.df$Class, sep = ";")  ## Generate search string for all classes
  groups <- strsplit(groups, split = ";")
  names(groups) <- paste(Sample.df$Sex, Sample.df$Class, sep = "_")
  group <- vector(mode = "character", length = length(colnames(sample.peaks)))
  for (i in 1:length(groups)) {
    rows_loop <- intersect(grep(groups[[i]][1], colnames(sample.peaks)), grep(groups[[i]][2], colnames(sample.peaks)))
    group[rows_loop] <- names(groups)[i]
  }
  group <- unlist(group)

  # Modify sample data to include the user defined exposure class
  Sample.data <- Sample.data[order(Sample.data[,1]), ]
  Sample.data <- setNames(cbind.data.frame(Sample.data[, 1], group, Sample.data[-1]), c(colnames(Sample.data[1]),
                                                                                        "Exposure Class", colnames(Sample.data[-1])))
  # Write the sample metadata to a CSV file.
  write.csv(Sample.data, "Sample_Metadata.csv")

}

#' @rdname output_MetaData
#' @export
output_MetaData.mass_all <- function(mSetObj, Peak.list, Sample.df, Sample.data)
  {

}

#' @rdname output_MetaData
#' @export
output_MetaData.default <- function(mSetObj, Peak.list, Sample.df, Sample.data)
{
  warning(paste("format_MetaData does not know how to handle object of class ",
                class(mSetObj),
                "and can only be used on classes pktable and mass_all"))
}
