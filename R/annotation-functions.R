#' @title IHL.search
#'
#' @export
#' @description Compare CAMERA isotope annotation with user-defined annotation library
#' @param Peak.list a table of class 'tbl_dbi', 'tbl_sql', 'tbl_lazy', or 'tbl' with samples as columns.  Should contain all output columns from XCMS and CAMERA, both metadata and sample data. Retention times must be in min.
#' @param Annotated.library a data frame with annotated metabolite entries. Must contain columns called 'Name', 'Formula', Molecular.Weight' and 'RT..Min.'.  Can contain additional info as separate columns.
#' @param rules a data frame containing the rule list used by CAMERA to annotate ion adducts and fragments.  Must contain the columns 'name','nmol','charge','massdiff','oidscore','quasi','ips'.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @param lib_db RSQLite connection
#' @return data frame containing the original table with added columns 'Name','MS.ID','Formula','Annotated.adduct' and any additional info columns from Annotated.library
#' @importFrom dplyr '%>%' select copy_to tbl between
#' @importFrom glue glue_collapse
#' @importFrom utils str txtProgressBar
search_IHL = function(Peak.list, Annotated.library, rules, search.par, ion.mode, lib_db) {
  search.list <- Peak.list %>% select(EIC_ID, mz, rt) %>% dplyr::collect()

  ## Creates full adduct list for all compounds in the Annotated library
  IHL <- Annotated.library[rep(seq_len(nrow(Annotated.library)), each = nrow(rules)), ]
  x <- rules$nmol
  IHL.temp <- sweep(IHL, 1, x, "*")
  x <- rules$massdiff
  if (ion.mode == "Positive") {
    IHL.temp <- sweep(IHL.temp, 1, x, "+")
    Ion.Mode <- "Pos"
    bin <- paste(Ion.Mode, "_", search.list$EIC_ID, "_", sep = "")

  } else {
    if (ion.mode == "Negative") {
      IHL.temp <- sweep(IHL.temp, 1, x, "+")
      IHL.temp <- sweep(IHL.temp, 1, -1, "*")
      Ion.Mode <- "Neg"
      bin <- paste(Ion.Mode, "_", search.list$EIC_ID, "_", sep = "")

    } else {
      stop("You must include the ionization mode!")
    }
  }
  x <- rules$charge
  IHL.adduct.data <- sweep(IHL.temp, 1, x, "/")
  IHL[, "mz"] <- IHL.adduct.data$Molecular.Weight
  IHL[, "adduct"] <- rules$name
  copy_to(lib_db, IHL, name = paste("Annotated Library", Ion.Mode, sep = "_"), temporary = FALSE, overwrite = TRUE)
  rm(IHL, IHL.adduct.data, IHL.temp)

  IHL <- tbl(lib_db, paste("Annotated Library", Ion.Mode, sep = "_"))
  d.mz <- search.list$mz * as.numeric(search.par[1, "ppm"])/10^6
  d.rt <- as.numeric(search.par[1, "rt"])

  # calculates the min and max range values for searching the in house library
  search.list$mz.min <- search.list$mz - d.mz
  search.list$mz.max <- search.list$mz + d.mz
  search.list$rt.min <- search.list$rt - d.rt
  search.list$rt.max <- search.list$rt + d.rt
  search.list$MS.ID = NA
  search.list$Formula = NA
  search.list$Name = NA
  search.list$Annotated.adduct = NA
  search.list$Conf.Level = NA
  search.list$FISh.Coverage = NA

  ## attempts to match all peaks against the In House Library ! Try to build a query using *apply functions to
  ## search all of the search list at once; should speed ! things up considerably
  ## i = 27 Used for debugging purposes
  total = nrow(search.list)
  cat("Annotating features against the In House Library.\n\n\n")
  pb = txtProgressBar(min = 0, max = total, style = 3)
  cnt = 1
  for (i in 1:nrow(search.list)) {
    mz.min = search.list$mz.min[i]
    mz.max = search.list$mz.max[i]
    rt.min = search.list$rt.min[i]
    rt.max = search.list$rt.max[i]
    test.list <- IHL %>% dplyr::filter(between(mz, mz.min, mz.max)) %>% dplyr::filter(between(RT..min., rt.min,
                                                                                              rt.max)) %>% dplyr::collect()
    if (nrow(test.list) == 0) {
      search.list$MS.ID[i] = paste(bin[i], "Unidentified", sep = "")
    } else {
      if (nrow(test.list) >= 1) { ##Matches all features against the In House Library
        #Prints the results of the match to console
        # cat("\n\n\nFeature annotated for match number ", cnt, ".\nMatch results below.\n\n\n", sep = "")
        # str(test.list)
        search.list$MS.ID[i] = paste(bin[i], "Annotated", sep = "")
        search.list$Name[i] = glue_collapse(test.list$Name, sep = ";", width = Inf, last = " or ")
        search.list$Formula[i] = glue_collapse(unique(gsub(" ", "", test.list$Formula)), sep = ";", width = Inf,
                                                last = " or ")
        search.list$Annotated.adduct[i] = glue_collapse(test.list$adduct, sep = ";", width = Inf, last = " or ")
        search.list$Conf.Level[i] = glue_collapse(test.list$Levels, sep = ";", width = Inf, last = " or ")
        search.list$FISh.Coverage[i] = glue_collapse(test.list$Fish.Coverage, sep = ";", width = Inf,
                                                      last = " or ")
        cnt = cnt + 1
      }
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)
  named.peak.list <- search.list[, 8:ncol(search.list)]
  colnames(named.peak.list)
  temp <- as.data.frame(Peak.list)
  colnames(Peak.list)
  if (nrow(temp) == nrow(named.peak.list)) {
    temp <- cbind(Peak.list, named.peak.list)
  }
  colnames(temp)
  return(temp)
}

#' @title Searches for Ion mode duplicates
#'
#' @export
#' @description Searches data tables from both ionization modes for duplicate entries
#' @param object used for method dispatch. Can be any object. Class must be "mz" or "monoMass"
#' @param Peak.list.pos Positive ionization mode data table
#' @param Peak.list.neg Negative ionization mode data table
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @return data frame containing the original table with added columns 'Name','MS.ID','Formula','Annotated.adduct' and any additional info columns from Annotated.library
search_IonDup <- function(object, Peak.list.pos,Peak.list.neg,search.par) {
  UseMethod("search_IonDup", object)
}

#' @method search_IonDup mz
#' @export
search_IonDup.mz  <- function(object,Peak.list.pos,Peak.list.neg,search.par) {
  ## Trim the feature table down to just those columns necessary for duplicate matching
  col.names <- c("Ion Mode", "EIC_ID", "mono_mass", "rt")
  mono.pos <- subset(Peak.list.pos, select = paste(col.names))

  mono.neg <- subset(Peak.list.neg, select = paste(col.names))

  mono.comb <- rbind(mono.pos, mono.neg)

  # get the distance matrices for mz and rt
  system.time({
    d.mz <- as.matrix(dist(mono.comb$mono_mass))
    dim(d.mz)
    v = as.numeric(as.character(mono.comb$mono_mass))
    v2 <- rep(v, each = dim(d.mz)[1])
    d.amu <- d.mz
    d.ppm <- d.amu/v2
    v = rep(10^6, length = ncol(d.ppm))
    v2 <- rep(v, each = dim(d.ppm)[1])
    d.mz <- d.ppm * v2
  })
  system.time(d.rt <- as.matrix(dist(mono.comb$rt)))

  # build the adjacency matrix
  m <- d.mz <= as.numeric(search.par[1, "ppm"]) & d.rt <= as.numeric(search.par[1, "rt"])

  # obtain the connected features
  g <- graph.adjacency(m)
  z <- clusters(g)$membership
  mono.comb$Duplicate_EIC <- ave(as.character(mono.comb$EIC_ID), z, FUN = function(s) paste(s, collapse = ","))
  mono.comb$Duplicate_ID <- z
  bin <- row.names(mono.comb)


  # Combine peaklists and add duplicate flag and ID
  colnames(Peak.list.pos)
  colnames(Peak.list.neg)
  colnames(Peak.list.neg) <- colnames(Peak.list.pos)
  Peak.list.combined <- as.data.frame(rbind(Peak.list.pos, Peak.list.neg))

  Peak.list.combined[, "Duplicate_ID"] <- mono.comb$Duplicate_ID
  Peak.list.combined[, "Duplicate_EIC"] <- mono.comb$Duplicate_EIC
  return(Peak.list.combined)
}

#' @method search_IonDup monoMass
#' @export
search_IonDup.monoMass  <- function(object,Peak.list.pos,Peak.list.neg,search.par) {
  ## Trim the feature table down to just those columns necessary for duplicate matching
  col.names <- c("Ion Mode", "EIC_ID", "mono_mass", "meanRT")
  mono.pos <- subset(Peak.list.pos, select = paste(col.names))

  mono.neg <- subset(Peak.list.neg, select = paste(col.names))

  mono.comb <- rbind(mono.pos, mono.neg)

  # get the distance matrices for mz and rt
  system.time({
    d.mz <- as.matrix(dist(mono.comb$mono_mass))
    dim(d.mz)
    v = as.numeric(as.character(mono.comb$mono_mass))
    v2 <- rep(v, each = dim(d.mz)[1])
    d.amu <- d.mz
    d.ppm <- d.amu/v2
    v = rep(10^6, length = ncol(d.ppm))
    v2 <- rep(v, each = dim(d.ppm)[1])
    d.mz <- d.ppm * v2
  })
  system.time(d.rt <- as.matrix(dist(mono.comb$meanRT)))

  # build the adjacency matrix
  m <- d.mz <= as.numeric(search.par[1, "ppm"]) & d.rt <= as.numeric(search.par[1, "rt"])

  # obtain the connected features
  g <- graph.adjacency(m)
  z <- clusters(g)$membership
  mono.comb$Duplicate_EIC <- ave(as.character(mono.comb$EIC_ID), z, FUN = function(s) paste(s, collapse = ","))
  mono.comb$Duplicate_ID <- z
  bin <- row.names(mono.comb)


  # Combine peaklists and add duplicate flag and ID
  colnames(Peak.list.pos)
  colnames(Peak.list.neg)
  colnames(Peak.list.neg) <- colnames(Peak.list.pos)
  Peak.list.combined <- as.data.frame(rbind(Peak.list.pos, Peak.list.neg))

  Peak.list.combined[, "Duplicate_ID"] <- mono.comb$Duplicate_ID
  Peak.list.combined[, "Duplicate_EIC"] <- mono.comb$Duplicate_EIC
  return(Peak.list.combined)
}

#' @title Searches for background components
#'
#' @export
#' @description Searches data table from samples for components coming from background
#' @param Peak.list data table containing sample data
#' @param Solv.list data table containing blank data
#' @param ... arguments to be passed to other functions
#' @return data frame containing the original table with added columns 'Name','MS.ID','Formula','Annotated.adduct' and any additional info columns from Annotated.library
search_solv <- function(Peak.list, Solv.list, search.par, ...) {
  UseMethod("search_solv", Peak.list)
}

#' @method search_solv mz
#' @export
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

#' @method search_solv monoMass
#' @export
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
