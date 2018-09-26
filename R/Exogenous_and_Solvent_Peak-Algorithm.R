#' @title Removes background components from Peak.list
#'
#' @export
#' @description Searches features or compounds in Peak.list that are present in process blanks and removes them.
#' Also flags compounds as exogenous, defined as mean abundance in endogenous samples below user-defined threshold.
#' @param Peak.list a named list of data frames (two per ionization mode) containing intensity matrices across all study samples and Pooled QCs and process blanks.  Names should be c('pos','neg','blanks_pos','blanks_neg').  Alternatively may use existing database connections by setting to NULL and specifying database parameters with ...
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param tbl.id character vector of table names to draw from databases.  First value should be table name from peak database, second second should be table name from solvent database. Default is NULL
#' @param method which method to use.  Can be monoMass or mz
#' @param lib.db character name of database to contain Solvent Library
#' @param ... Arguments to pass parameters to database and annotation functions
#' @return nested list a list for each ionization mode, each containing a list of two dataframes: the first contains the intensity matrix for the peaklist with solvent peaks removed, the second contains the intensity matrix for the solvent peaks
#' @importFrom plyr llply
#' @importFrom dplyr filter
remove_background_peaks = function(Peak.list, Sample.df, search.par, tbl.id, method, lib.db, ...) {
    if (missing(tbl.id))
        tbl.id = NULL
    if (missing(method))
      stop("Need to specify a method!", call. = FALSE)
    if (is.null(Peak.list) && is.null(tbl.id)) {
        stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
    }
    if (is.null(Peak.list)) {
        mydbs <- connect_lumadb(db.list, db.dir)
        nm = names(mydbs)[-1]
        nm1 = nm[-grep("blanks", nm)]  # Return all sample databases
        nm2 = nm[grep("blanks", nm)]  #Return all blank databases
        list1 <- llply(nm1, .fun = function(x) read_tbl(tbl.id[1], mydbs[[x]]))
        list2 <- llply(nm2, .fun = function(x) read_tbl(tbl.id[2], mydbs[[x]]))
        names(list1) <- gsub("_db", "\\1", nm1)
        names(list2) <- gsub("_db", "\\1", nm2)
        Peak.list <- list1
        Solv.list <- list2
    } else {
        nm = names(Peak.list)
        nm1 = nm[-grep("blanks", nm)]  # Return all sample databases
        nm2 = nm[grep("blanks", nm)]  #Return all blank databases
        Solv.list <- Peak.list[nm2]
        Peak.list <- Peak.list[nm1]
    }
    lib_db <- connect_libdb(lib.db = lib.db, db.dir = db.dir)
    peak_db <- connect_lumadb(db.list = db.list, db.dir = db.dir, new.db = new.db)

    class(method) <- method
    masterlist <- search_solv(method, Peak.list, Solv.list, Sample.df, search.par, ..., lib_db)
    return(masterlist)
}

#' @title Combines Peak.lists simply
#'
#' @export
#' @description Combines Peak.lists from both ionization modes \(positive and negative\) for a first look at class separations, e.g. with PCA.
#' @param Peak.list a named list of data frames \(one per ionization mode\) containing intensity matrices across all study samples and Pooled QCs.
#' Names must be c('Positive','Negative').  Alternatively may use existing database connections by setting to NULL and specifying database parameters with ...
#' @param tbl.id character vector of table names to draw from database.  First value should be table name for positive mode, second should be table name for negative mode. Default is NULL
#' @param ... Arguments to pass parameters to database functions
#' @return NULL testing
pre_combine_ion_modes = function(Peak.list, tbl.id, ...) {
    if (missing(tbl.id))
        tbl.id = NULL
    if (is.null(tbl.id) && is.null(Peak.list)) {
        stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
    }
    if (is.null(Peak.list)) {
        df1 <- read_tbl(tbl.id[1], peak.db)
        df2 <- read_tbl(tbl.id[2], peak.db)
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
