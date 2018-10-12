#' @title Formatting of metabolite data for MetaboAnalystR.
#'
#' @description This function initializes objects that will hold the metabolite data, formats peak intensity data into one of the formats acceptable by MetaboAnalystR, and sets the metabolite data object.
#' @param Peak.list data frame containing combined ion mode peaklist with ion mode duplicates removed.
#' @param Sample.df data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Sample.data ata frame with phenotype data as columns and a row for each study sample.  First column must be a unique sample identifier with the header 'CT-ID'.  Phenotype columns may vary, but must include two columns called 'Plate Number' and 'Plate Position' for determining run order.
#' @return mSetObbj
#' @importFrom
format_MetaboAnalystR <- function(Peak.list, Sample.df, Sample.data)
{
  ##-----------------------------------------------------------------------------------------
  ## Initialize the metabolite object.
  ##-----------------------------------------------------------------------------------------
  dataSet <- list();
  dataSet$type <- data.type;
  dataSet$design.type <- "regular";    # one factor to two factor
  dataSet$cls.type <- "disc";          # default until specified otherwise
  dataSet$format <- "rowu";
  dataSet$paired <- paired;
  analSet <- list();
  analSet$type <- anal.type;

  mSetObj <- list();
  mSetObj$dataSet <- dataSet;
  mSetObj$analSet <- analSet;
  mSetObj$imgSet <- list();
  mSetObj$msgSet <- list();                                 # store various message during data processing
  mSetObj$msgSet$msg.vec <- vector(mode = "character");     # store error messages
  mSetObj$cmdSet <- vector(mode = "character");             # store R command


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
  sample.ID  <-  as.numeric(sub("\\D*(\\d{6}).*", "\\1", colnames(sample.peaks)))
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
  ## Set the metabolite data object.
  ##-----------------------------------------------------------------------------------------
  mSetObj$dataSet$cls.type <-  "disc"
  mSetObj$dataSet$format <-  "rowu"
  dat  <-  tsample.peaks
  mSetObj$dataSet$type.cls.lbl <- class(dat[, 1]);

  mSetObj$dataSet$orig.cls  <-  mSetObj$dataSet$cls  <-  as.factor(as.character(dat[, 2]));
  mSetObj$dataSet$orig.cls  <-  mSetObj$dataSet$cls  <-  as.factor(as.character(dat[, 1]));
  mSetObj$dataSet$orig  <-  dat[, -c(1,2)];

  return(mSetObj)
}
