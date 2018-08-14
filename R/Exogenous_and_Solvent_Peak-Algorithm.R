#' @title Exogenous and background peak algorithm
#'
#' @export
#' @description Searches mz/rt features for exogenous features and removes background features.  Exogenous is defined as mean abundance in endogenous samples below user-defined threshold.  Background determined from process blanks.
#' @param Peak.list a named list of data frames (two per ionization mode) containing intensity matrices across all study samples and Pooled QCs and process blanks.  Names should be c('pos','neg','blanks_pos','blanks_neg').  Alternatively may use existing database connections by setting to NULL and specifying database parameters with ...
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param ion.id character vector specifying identifier in filename designating positive or negative ionization mode or both.  Positive identifier must come first. Default is c('Pos','Neg')
#' @param QC.id character vector specifying identifier in filename designating a Pooled QC sample.  Only the first value will be used.  Default is 'Pooled_QC_'
#' @param MB.id character vector specifying identifier in filename designating a Method Blank.  Only the first value will be used. Default is '^MB'
#' @param tbl.id character vector of table names to draw from databases.  First value should be table name from peak database, second second should be table name from solvent database. Default is NULL
#' @param db.id character vector specifying identifiers in database names designating sample \[1\] and blank \[2\] databases.  Default is c('Peaklist','Blanks')
#' @param ion.modes a character vector defining the ionization mode.  Must be either 'Positive', 'Negative' or both.
#' @param method which method to use.  Can be monoMass or mz
#' @param lib.db character name of database to contain Solvent Library
#' @param ... Arguments to pass parameters to database functions
#' @return nested list a list for each ionization mode, each containing a list of two dataframes: the first contains the intensity matrix for the peaklist with solvent peaks removed, the second contains the intensity matrix for the solvent peaks
#' @importFrom plyr llply
#' @importFrom dplyr filter
remove_background_peaks = function(Peak.list, Sample.df, search.par, ion.id, QC.id, MB.id, tbl.id, db.id,
    ion.modes, method, lib.db, ...) {
    if (missing(ion.id)) {
        ion.id = c("Pos", "Neg")
    }
    if (missing(QC.id))
        QC.id = "Pooled_QC_"
    if (missing(MB.id))
        MB.id = "^MB"
    if (missing(tbl.id))
        tbl.id = NULL
    if (missing(db.id))
        db.id = c("Peaklist", "Blanks")
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

    class(Peak.list) <- method
    masterlist <- search_solv(Peak.list, Solv.list, Sample.df, search.par, ion.id, QC.id, MB.id, tbl.id, db.id,
                              ion.modes,lib_db)
    return(masterlist)
}

#' @title Combines ion mode data simply
#'
#' @export
#' @description Combines peaklists from both ionization modes \(positive and negative\) for a first look at class separations
#' @param Peak.list a named list of data frames \(one per ionization mode\) containing intensity matrices across all study samples and Pooled QCs.  Names must be c('Positive','Negative').  Alternatively may use existing database connections by setting to NULL and specifying database parameters with ...
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

search_solv <- function(Peak.list, ...) {
  UseMethod("search_solv", Peak.list)
}

search_solv.mz <- function(Peak.list, Solv.list, Sample.df, search.par, ion.id, QC.id, MB.id, tbl.id, db.id,
                           ion.modes,lib_db) {
  endo.groups <- as.matrix(paste(Sample.df[which(Sample.df[, "Endogenous"] == TRUE), "Sex"], Sample.df[which(Sample.df[,
                                                                                                                       "Endogenous"] == TRUE), "Class"], sep = "_"))
  list.length = length(ion.modes)
  mylist <- NULL
  masterlist <- NULL
  for (i in 1:list.length) {
    ## Selects the filebase to use for data processing
    cur.ion = ion.modes[i]
    file.base = paste(db.id[1], ion.id[i], sep = "_")
    blank.base = paste(db.id[2], ion.id[i], sep = "_")
    cur.Peaklist <- Peak.list[[i]]
    cur.Solvlist <- Solv.list[[i]]

    res <- as.numeric(unlist(apply(endo.groups, 1, function(X) grep(X, colnames(cur.Peaklist)))))
    Endo.list <- cur.Peaklist[, res]
    search.list <- cur.Peaklist %>% select(EIC_ID, mz, rt) %>% dplyr::collect()

    # Calculate the mean of QC values for each compound
    res <- lapply(colnames(cur.Peaklist), function(ch) grep(QC.id, ch))
    QC.list <- cur.Peaklist[sapply(res, function(x) length(x) > 0)]
    QCmean <- rowMeans(QC.list)
    cur.Peaklist[, "Endogenous_flag"] <- NULL
    Endo.mean <- apply(Endo.list, 1, mean)
    Endo.lim <- as.numeric(search.par[1, "Endogenous"])
    Endo.flag <- Endo.mean > Endo.lim
    cur.Peaklist[, "Endogenous_flag"] <- Endo.flag

    # Calculate the mean of MB values for each compound
    res <- lapply(colnames(cur.Solvlist), function(ch) unique(grep(MB.id, ch)))  #Flags the EIC_ID and Isotope cluster columns
    Solvent.new.list <- cur.Solvlist[sapply(res, function(x) length(x) > 0)]  #Extracts the EIC_ID and Isotope cluster columns
    MBmean <- rowMeans(Solvent.new.list)
    cur.Solvlist$mean <- MBmean
    if (cur.ion == "Positive") {
      cat("Removing Background Compounds in Positive mode:")
      copy_to(lib_db, cur.Solvlist, name = "Pos_list", temporary = FALSE, overwrite = TRUE)
      IHL <- tbl(lib_db, "Pos_list")
    } else {
      if (cur.ion == "Negative") {
        cat("Removing Background Compounds in Negative mode:")
        copy_to(lib_db, cur.Solvlist, name = "Neg_list", temporary = FALSE, overwrite = TRUE)
        IHL <- tbl(lib_db, "Neg_list")
      } else {
        stop("Ionization mode must be Positive or Negative. It is case sensitive!", call. = FALSE)
      }
    }
    # calculates the min and max range values for searching against the solvent list
    d.mz <- search.list$mz * as.numeric(search.par[1, "ppm"])/10^6
    d.rt <- as.numeric(search.par[1, "rt"])

    search.list$mz.min <- search.list$mz - d.mz
    search.list$mz.max <- search.list$mz + d.mz
    search.list$rt.min <- search.list$rt - d.rt
    search.list$rt.max <- search.list$rt + d.rt
    search.list$sample.mean = QCmean

    bin = vector(mode = "logical", length = search.list %>% nrow)
    Solvent.ratio <- as.numeric(search.par[1, "Solvent"])
    ## attempts to match all peaks against the In House Library ! Try to build a query using *apply functions to
    ## search all of the search list at once; should speed ! things up considerably
    cnt = 1
    total <- nrow(search.list)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    for (j in 1:total) {
      mz.min = search.list$mz.min[j]
      mz.max = search.list$mz.max[j]
      rt.min = search.list$rt.min[j]
      rt.max = search.list$rt.max[j]
      test.list <- IHL %>% filter(between(mz, mz.min, mz.max)) %>% filter(between(rt, rt.min, rt.max)) %>%
        dplyr::collect()
      if (nrow(test.list) == 0) {
      } else {
        if (nrow(test.list) == 1) {
          if (search.list$sample.mean[j]/test.list$mean <= Solvent.ratio) {
            temp <- test.list[, c("MS.ID", "Formula", "Name", "Annotated.adduct", "Conf.Level", "FISh.Coverage",
                                  "isotopes", "adduct")]
            bin[j] = TRUE
            cnt = cnt + 1
          }

        } else {
          if (nrow(test.list) >= 1) {
            if (search.list$sample.mean[j]/max(test.list$mean) <= Solvent.ratio) {
              temp <- test.list[, c("MS.ID", "Formula", "Name", "Annotated.adduct", "Conf.Level", "FISh.Coverage",
                                    "isotopes", "adduct")]
              bin[j] = TRUE
              cnt = cnt + 1
            }
          }
        }
      }
      setTxtProgressBar(pb, j)
    }
    cur.Peaklist$Solv_flag <- bin

    length(which(cur.Peaklist$Solv_flag %in% FALSE))  # The number of features matched between Peaklist and Endolist
    cur.Peaklist.trimmed <- cur.Peaklist[which(cur.Peaklist$Solv_flag %in% FALSE), ]
    Solv.list.trimmed <- cur.Peaklist[which(cur.Peaklist$Solv_flag %in% TRUE), ]
    mylist <- list(cur.Peaklist.trimmed, Solv.list.trimmed)
    names(mylist) <- c("Peaklist", "Solvent_Peaks")
    masterlist[[i]] <- mylist
  }
  names(masterlist) <- ion.modes
  return(masterlist)
}

search_solv.monoMass <- function(Peak.list, Solv.list, Sample.df, search.par, ion.id, QC.id, MB.id, tbl.id, db.id,
                                 ion.modes,lib_db) {
  endo.groups <- as.matrix(paste(Sample.df[which(Sample.df[, "Endogenous"] == TRUE), "Sex"], Sample.df[which(Sample.df[,
                                                                                                                       "Endogenous"] == TRUE), "Class"], sep = "_"))
  list.length = length(ion.modes)
  mylist <- NULL
  masterlist <- NULL
  for (i in 1:list.length) {
    ## Selects the filebase to use for data processing
    cur.ion = ion.modes[i]
    file.base = paste(db.id[1], ion.id[i], sep = "_")
    blank.base = paste(db.id[2], ion.id[i], sep = "_")
    cur.Peaklist <- Peak.list[[i]]
    cur.Solvlist <- Solv.list[[i]]

    res <- as.numeric(unlist(apply(endo.groups, 1, function(X) grep(X, colnames(cur.Peaklist)))))
    Endo.list <- cur.Peaklist[, res]
    search.list <- cur.Peaklist %>%
      select(EIC_ID, mono_mass, meanRT) %>%
      dplyr::collect()

    # Calculate the mean of QC values for each compound
    res <- lapply(colnames(cur.Peaklist), function(ch) grep(QC.id, ch))
    QC.list <- cur.Peaklist[sapply(res, function(x) length(x) > 0)]
    QCmean <- rowMeans(QC.list)
    cur.Peaklist[, "Endogenous_flag"] <- NULL
    Endo.mean <- apply(Endo.list, 1, mean)
    Endo.lim <- as.numeric(search.par[1, "Endogenous"])
    Endo.flag <- Endo.mean > Endo.lim
    cur.Peaklist[, "Endogenous_flag"] <- Endo.flag

    # Calculate the mean of MB values for each compound
    res <- lapply(colnames(cur.Solvlist), function(ch) unique(grep(MB.id, ch)))  #Flags the EIC_ID and Isotope cluster columns
    Solvent.new.list <- cur.Solvlist[sapply(res, function(x) length(x) > 0)]  #Extracts the EIC_ID and Isotope cluster columns
    MBmean <- rowMeans(Solvent.new.list)
    cur.Solvlist$mean <- MBmean
    if (cur.ion == "Positive") {
      cat("Removing Background Compounds in Positive mode:")
      copy_to(lib_db, cur.Solvlist, name = "Pos_list", temporary = FALSE, overwrite = TRUE)
      IHL <- tbl(lib_db, "Pos_list")
    } else {
      if (cur.ion == "Negative") {
        cat("Removing Background Compounds in Negative mode:")
        copy_to(lib_db, cur.Solvlist, name = "Neg_list", temporary = FALSE, overwrite = TRUE)
        IHL <- tbl(lib_db, "Neg_list")
      } else {
        stop("Ionization mode must be Positive or Negative. It is case sensitive!", call. = FALSE)
      }
    }
    # calculates the min and max range values for searching against the solvent list
    d.mz <- search.list$mono_mass * as.numeric(search.par[1, "ppm"])/10^6
    d.rt <- as.numeric(search.par[1, "rt"])

    search.list$mz.min <- search.list$mono_mass - d.mz
    search.list$mz.max <- search.list$mono_mass + d.mz
    search.list$rt.min <- search.list$meanRT - d.rt
    search.list$rt.max <- search.list$meanRT + d.rt
    search.list$sample.mean = QCmean

    bin = vector(mode = "logical", length = search.list %>% nrow)
    Solvent.ratio <- as.numeric(search.par[1, "Solvent"])
    ## attempts to match all peaks against the In House Library ! Try to build a query using *apply functions to
    ## search all of the search list at once; should speed ! things up considerably
    cnt <- 1
    total <- nrow(search.list)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    for (j in 1:total) {
      mz.min = search.list$mz.min[j]
      mz.max = search.list$mz.max[j]
      rt.min = search.list$rt.min[j]
      rt.max = search.list$rt.max[j]
      test.list <- IHL %>%
        filter(between(mono_mass, mz.min, mz.max)) %>%
        filter(between(meanRT, rt.min, rt.max)) %>%
        dplyr::collect()
      if (nrow(test.list) == 0) {
      } else {
        if (nrow(test.list) == 1) {
          if (search.list$sample.mean[j]/test.list$mean <= Solvent.ratio) {
            temp <- test.list[, c("MS.ID", "Formula", "Name", "Annotated.adduct", "Conf.Level", "FISh.Coverage",
                                  "isotopes", "adduct")]
            bin[j] = TRUE
            cnt = cnt + 1
          }

        } else {
          if (nrow(test.list) >= 1) {
            if (search.list$sample.mean[j]/max(test.list$mean) <= Solvent.ratio) {
              temp <- test.list[, c("MS.ID", "Formula", "Name", "Annotated.adduct", "Conf.Level", "FISh.Coverage",
                                    "isotopes", "adduct")]
              bin[j] = TRUE
              cnt = cnt + 1
            }
          }
        }
      }
      setTxtProgressBar(pb, j)
    }
    cur.Peaklist$Solv_flag <- bin

    length(which(cur.Peaklist$Solv_flag %in% FALSE))  # The number of features matched between Peaklist and Endolist
    cur.Peaklist.trimmed <- cur.Peaklist[which(cur.Peaklist$Solv_flag %in% FALSE), ]
    Solv.list.trimmed <- cur.Peaklist[which(cur.Peaklist$Solv_flag %in% TRUE), ]
    mylist <- list(cur.Peaklist.trimmed, Solv.list.trimmed)
    names(mylist) <- c("Peaklist", "Solvent_Peaks")
    masterlist[[i]] <- mylist
  }
  names(masterlist) <- ion.modes
  return(masterlist)

}
