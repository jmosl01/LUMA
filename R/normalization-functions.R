#' @title Replaces zeros and normalizes
#'
#' @export
#' @description Replaces zero values with half of the minimum value from all
#'   study samples and normalizes to unity
#' @param Peak.list a data frame containing combined ion mode peaklist with ion
#'   mode duplicates removed.  Alternatively can be retrieved from databases.
#'   Default is \code{NULL}
#' @param Sample.df a data frame with class info as columns.  Must contain a
#'   separate row entry for each unique sex/class combination. Must contain the
#'   columns \code{"Sex","Class","n","Endogenous"}.
#' @param tbl.id character vector of table names to draw from databases.  First
#'   value should be table name from positive ionization peak database, second
#'   should be table name from negative ionization peak database. Default is
#'   \code{NULL}
#' @param QC.id character vector specifying identifier in filename designating a
#'   Pooled QC sample.  Only the first value will be used.  Default is
#'   \code{"Pooled_QC_"}
#' @param mult numeric multiplier to be added to unit normalization.  Default is
#'   \code{1000}
#' @param ... Arguments to pass parameters to database functions
#' @return data frame containing normalized intensity matrix for all samples
#'   plus Pooled QCs
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'
#'   file <- system.file("extdata","Sample_Class.txt", package = "lcmsfishdata")
#'   Sample.df <- read.table(file, header = TRUE, sep = "\t")
#'   Peak.list <- lcmsfishdata::Peaklist_db$Peaklist_Combined_FINAL
#'   QC.ID <- "Pooled_QC"
#'   mult <- 1
#'
#'   test <- replace_zeros_normalize(Peak.list = Peak.list, Sample.df =
#'   Sample.df, QC.id = QC.ID, mult = mult)
#'   sum(test[,grep(QC.ID, colnames(test))]) # There are 11 Pooled QC samples,
#'   # each with unit normalized intensity values
#'
#'   mult2 <- 1000
#'   test2 <- replace_zeros_normalize(Peak.list = Peak.list, Sample.df =
#'   Sample.df, QC.id = QC.ID, mult = mult2)
#'   sum(test2[,grep(QC.ID, colnames(test))]) #Now they are each normalized to 1000
#'
#'   test_pos <- test[grep("Pos",test[["Ion.Mode"]]),grep(QC.ID, colnames(test))]
#'   test_neg <- test[grep("Neg",test[["Ion.Mode"]]),grep(QC.ID, colnames(test))]
#'
#'   sum(test_pos) == sum(test_neg) #Ion modes are given equal weighting
#'
#' }
replace_zeros_normalize = function(Peak.list = NULL, Sample.df, tbl.id = NULL, QC.id = "Pooled_QC_",
                                   mult = 1000, ...) {
  #Set default values
  if (missing(Peak.list))
    Peak.list = NULL
  if (missing(tbl.id))
    tbl.id = NULL
  if (missing(QC.id))
    QC.id = "Pooled_QC"
  if (missing(mult))
    mult = 1000
  if (is.null(tbl.id) && is.null(Peak.list)) {
    stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
  }
  if (is.null(Peak.list)) {
    Peak.list <- read_tbl(tbl.id, ...)
  }
  sexes <- unique(paste(Sample.df$Sex, "_", sep = ""))  ## Generate search string for all sexes
  samples <- vector(mode = "character", length = length(colnames(Peak.list)))
  for (i in 1:length(sexes)) {
    rows_loop <- grep(sexes[i], colnames(Peak.list))
    samples[rows_loop] <- sexes[i]
  }

  res <- samples %in% sexes
  sum.range.list <- Peak.list[, res]
  # return(sum.range.list) Checks out

  # Get all study samples for each ion mode separately Positive mode
  sum.range.pos <- split(sum.range.list, as.factor(Peak.list$Ion.Mode))$Pos
  sum.range.pos[sum.range.pos == 0] <- NA
  min.pos <- min(sum.range.pos, na.rm = TRUE)
  # Negative mode
  sum.range.neg <- split(sum.range.list, as.factor(Peak.list$Ion.Mode))$Neg
  sum.range.neg[sum.range.neg == 0] <- NA
  min.neg <- min(sum.range.neg, na.rm = TRUE)

  # Replace NA values with 0.5*min values from all samples Positive mode
  w.min.pos <- replace(sum.range.pos, is.na(sum.range.pos), 0.5 * min.pos)
  # Negative mode
  w.min.neg <- replace(sum.range.neg, is.na(sum.range.neg), 0.5 * min.neg)

  # Unit Normalize all samples by their sum Positive mode
  norm.pos <- normalize_unit(w.min.pos)
  normed <- lapply(norm.pos, function(x) x * mult/2)
  sum.range.pos <- as.data.frame(normed)
  # Negative mode
  norm.neg <- normalize_unit(w.min.neg)
  normed <- lapply(norm.neg, function(x) x * mult/2)
  sum.range.neg <- as.data.frame(normed)
  sum.range.list <- rbind(sum.range.pos, sum.range.neg)

  # Normalize the QCs separately
  sexes <- QC.id  ## Generate search string for pooled QCs
  samples <- vector(mode = "character", length = length(colnames(Peak.list)))
  for (i in 1:length(sexes)) {
    rows_loop <- grep(sexes[i], colnames(Peak.list))
    samples[rows_loop] <- sexes[i]
  }

  res <- samples %in% sexes
  QC.range.list <- Peak.list[, res]
  QC.range.pos <- split(QC.range.list, as.factor(Peak.list$Ion.Mode))$Pos
  QC.range.pos[QC.range.pos == 0] <- NA

  QC.range.neg <- split(QC.range.list, as.factor(Peak.list$Ion.Mode))$Neg
  QC.range.neg[QC.range.neg == 0] <- NA

  min.pos <- min(QC.range.pos, na.rm = TRUE)
  min.neg <- min(QC.range.neg, na.rm = TRUE)
  # Positive mode
  w.min.pos <- replace(QC.range.pos, is.na(QC.range.pos), 0.5 * min.pos)
  norm.pos <- normalize_unit(w.min.pos)
  normed <- lapply(norm.pos, function(x) x * mult/2)
  QC.range.pos <- as.data.frame(normed)

  # Negative mode
  w.min.neg <- replace(QC.range.neg, is.na(QC.range.neg), 0.5 * min.neg)
  norm.neg <- normalize_unit(w.min.neg)
  normed <- lapply(norm.neg, function(x) x * mult/2)
  QC.range.neg <- as.data.frame(normed)
  QC.range.list <- rbind(QC.range.pos, QC.range.neg)

  Norm.Peaklist <- Peak.list
  Norm.Peaklist[, colnames(sum.range.list)] <- sum.range.list
  Norm.Peaklist[, colnames(QC.range.list)] <- QC.range.list
  return(Norm.Peaklist)
}

#' @title normalize_unit
#'
#' @description normalizes each column to unity. See
#'   \code{replace_zeros_normalize} and \code{NormalizePeaklists} for examples
#'   that use this function.
#' @param x dataframe containing only numeric vectors representing samples with
#'   metabolite intensities.  Cannot contain zeros or NA values
#' @return An array with the same shape as \code{x}, but with the summary
#'   statistics swept out
normalize_unit <- function(x) {
  x <- sweep(x, 2, apply(x, 2, sum), "/")
  return(x)
}

