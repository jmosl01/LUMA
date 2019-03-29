#' @title Combines Peak.lists from two ion modes
#'
#' @export
#' @description Combines the two single adduct ion mode feature tables into one for EIC plotting
#' @param Peak.list a named list of data frames (one per ionization mode) containing intensity matrices across all study samples and Pooled QCs.
#' Names must be c('Positive','Negative').  Alternatively may use existing database connections by setting to NULL
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param tbl.id character vector of table names to draw from databases.  First value should be table containing compounds from positive ionization, second should be table containing compounds from negative ionization. Default is NULL
#' @param method which method to apply to search for duplicate entries.
#' See search_IonDup for details.
#' @param peak.db Formal class SQLiteConnection
#' @param ... Arguments to pass parameters to search_IonDup
#' @return data frame containing the intensity matrix for the peaklist with Duplicate IDs
#' @importFrom igraph clusters graph.adjacency
#' @importFrom stats ave
combine_ion_modes = function(Peak.list, search.par, tbl.id, method, peak.db, ...) {

  #Set default values
    if (missing(Peak.list))
        Peak.list = NULL
    if (missing(tbl.id))
        tbl.id = NULL
    if (is.null(Peak.list) && is.null(tbl.id)) {
        stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
    }
    if (is.null(Peak.list)) {
        df1 <- read_tbl(mytable = tbl.id[1], peak.db = peak.db)
        df2 <- read_tbl(mytable = tbl.id[2], peak.db = peak.db)
        Peak.list.pos <- df1
        Peak.list.neg <- df2
    } else {
        # Error check
        if (length(grep("Positive", names(Peak.list))) == 0 || length(grep("Negative", names(Peak.list))) == 0) {
            stop("Peaklist must be a named list. Try \n> names(Peak.list) = c(\"Positive\",\"Negative\")", call. = FALSE)
        }

        Peak.list.pos <- Peak.list["Positive"]
        Peak.list.neg <- Peak.list["Negative"]
    }
    Peak.list.pos[, "Ion Mode"] <- "Pos"
    Peak.list.neg[, "Ion Mode"] <- "Neg"

    class(method) <- method
    Peak.list.combined <- search_IonDup(method,Peak.list.pos,Peak.list.neg,search.par,...)
    return(Peak.list.combined)
}

#' @title Removes Ion Mode Duplicates
#'
#' @export
#' @description Removes the ion mode duplicates based on user-modified indices after visual inspection of EIC_plots
#' @param Peak.list a data frame containing combined ion mode peaklist with Duplicate IDs.
#' Alternatively can be retrieved from databases.  Default is NULL
#' @param Key.list list containing two numeric vectors \(one per ionization mode\) of ion mode duplicates to keep. Key (vector) corresponding to positive mode duplicates should be first.
#' @param tbl.id character string corresponding to table name to draw from database.
#' Default is NULL
#' @param ... Arguments to pass parameters to database functions
#' @return NULL testing
remove_ion_dup = function(Peak.list, Key.list, tbl.id, ...) {
    if (missing(Peak.list))
        Peak.list = NULL
    if (missing(tbl.id))
        tbl.id = NULL
    if (is.null(Peak.list) && is.null(tbl.id)) {
        stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
    }
    if (is.null(Peak.list)) {
        Peak.list <- read_tbl(tbl.id, ...)
    }

    # Get the unique list of Duplicate IDs
    x <- sapply(Peak.list$Duplicate_ID, function(x) sum(as.numeric(Peak.list$Duplicate_ID == x)))

    drops <- Peak.list$Duplicate_ID[x == 1] # Duplicate IDs which only appear once


    List.ID <- Peak.list$Duplicate_ID # List of all duplicate IDs

    res <- !List.ID %in% drops
    Un.ID <- unique(Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE)]) # List of unique duplicate IDs to get EICs


    # List of duplicate IDs for both positive and negative modes
    Dup.ID.Pos <- Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE) & sapply(Peak.list$`Ion Mode`, function(x) x ==
        "Pos")]
    Dup.ID.Neg <- Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE) & sapply(Peak.list$`Ion Mode`, function(x) x ==
        "Neg")]

    # Initialize the drop vectors
    x.pos <- rep(1, length.out = length(Dup.ID.Pos))  #needs to be as long as the number of positive plots in the PDF file
    x.neg <- rep(1, length.out = length(Dup.ID.Neg))  #needs to be as long as the number of negative plots in the PDF file
    EICs <- data.frame(Peak.list[Peak.list$Duplicate_ID %in% Un.ID, "EIC_ID"])

    # Key vectors of plots to keep
    if (length(Key.list) == 2) {
        key.pos <- Key.list[[1]]
        key.neg <- Key.list[[2]]
    } else {
        stop("Key.list must contain a list of length 2.  Each element should be a vector of duplicate IDs to keep.")
    }

    # 1 selects for removal, 0 for keeping; only keep the EICs that were not removed by the user
    key.pos <- as.numeric(row.names(key.pos))
    key.neg <- as.numeric(row.names(key.neg))
    drop.pos <- replace(x.pos, key.pos, 0)  #Creates the drop vector for positive plots
    drop.neg <- replace(x.neg, key.neg, 0)  #Creates the drop vector for negative plots

    drops <- c(drop.pos, drop.neg)
    drops <- sapply(drops, function(x) x > 0)
    EICs <- EICs[drops, ]

    ## The number of spectra to keep is equal to what you told the code
    Peak.list.trimmed <- Peak.list[!Peak.list$EIC_ID %in% EICs, ]
    return(Peak.list.trimmed)

}


#' @title Combines Peak.lists simply
#'
#' @export
#' @description Combines Peak.lists from both ionization modes \(positive and negative\) for a first look at class separations, e.g. with PCA.
#' @param Peak.list a named list of data frames \(one per ionization mode\) containing intensity matrices across all study samples and Pooled QCs.
#' Names must be c('Positive','Negative').  Alternatively may use existing database connections by setting to NULL and specifying database parameters with ...
#' @param tbl.id character vector of table names to draw from database.  First value should be table name for positive mode, second should be table name for negative mode. Default is NULL
#' @param ... arguments to pass to database functions
#' @return NULL testing
pre_combine_ion_modes = function(Peak.list = NULL, tbl.id, ...) {
  if (missing(tbl.id))
    tbl.id = NULL
  if (is.null(tbl.id) && is.null(Peak.list)) {
    stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
  }
  if (is.null(Peak.list)) {
    df1 <- read_tbl(tbl.id[1], ...)
    df2 <- read_tbl(tbl.id[2], ...)
    Peak.list.pos <- df1
    Peak.list.neg <- df2
  } else {
    # Error check
    if (length(grep("Positive", names(Peak.list))) == 0 || length(grep("Negative", names(Peak.list))) == 0) {
      stop("Peaklist must be a named list. Try \n> names(Peak.list) = c(\"Positive\",\"Negative\")", call. = FALSE)
    }

    Peak.list.pos <- Peak.list["Positive"]
    Peak.list.neg <- Peak.list["Negative"]
  }
  Peak.list.pos[, "Ion Mode"] <- "Pos"
  Peak.list.neg[, "Ion Mode"] <- "Neg"


  Col.names.pos <- colnames(Peak.list.pos)
  colnames(Peak.list.neg) <- Col.names.pos

  Peak.list.combined <- rbind(Peak.list.pos, Peak.list.neg)
  return(Peak.list.combined)
}