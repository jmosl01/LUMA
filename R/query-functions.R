#' @title Matches Peak.list annotations against Annotated Library
#'
#' @export
#' @description Compares isotope and adduct annotations within Peak.list (inherited from CAMERA) to user-defined annotation library and returns annotation results
#' @param Peak.list a table of class 'tbl_dbi', 'tbl_sql', 'tbl_lazy', or 'tbl' with samples as columns.  Should contain all output columns from XCMS and CAMERA, both metadata and sample data. Retention times must be in min.
#' @param Annotated.library a data frame with annotated metabolite entries. Must contain the first four columns called 'Name', 'Formula', Molecular.Weight' and 'RT..Min.', respectively.  Can contain additional columns.
#' @param rules a data frame containing the rule list used by CAMERA to annotate ion adducts and fragments.  Must contain the columns 'name','nmol','charge','massdiff','oidscore','quasi','ips'.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @param lib_db RSQLite connection
#' @return data frame containing the original table with added columns 'Name','MS.ID','Formula','Annotated.adduct' and any additional info columns from Annotated.library
#' @importFrom dplyr '%>%' select copy_to tbl between
#' @importFrom utils str txtProgressBar
match_Annotation = function(Peak.list, Annotated.library, rules, search.par, ion.mode, lib_db) {

  myresults <- .gen_IHL(Peak.list, Annotated.library, rules, ion.mode, lib_db)
  search.list <- myresults[[1]]
  bin <- myresults[[2]]
  myion.mode <- myresults[[3]]

  ## Creates delta mz and rt for all compounds in the Search.list
  d.mz <- search.list[["mz"]] * search.par[["ppm"]][1]/10^6
  d.rt <- search.par[["rt"]][1]

  # calculates the min and max range values for searching the in house library
  mz.min <- search.list[["mz"]] - d.mz
  mz.max <- search.list[["mz"]] + d.mz
  rt.min <- search.list[["rt"]] - d.rt
  rt.max <- search.list[["rt"]] + d.rt

  # Sets up the match results columns including all phenodata
  search.list["MS.ID"] = NA_character_
  search.list["Name"] = NA_character_
  search.list["Formula"] = NA_character_
  search.list["Molecular.Weight"] = NA_character_
  search.list["RT..Min."] = NA_character_
  search.list["Annotated.adduct"] = NA_character_
  search.list["match.mz"] = NA_character_
  search.list["match.adduct"] = NA_character_
  search.list[,seq(from = 12, to = 12 - 1 + length(colnames(Library.phenodata)), by = 1)] = NA_character_
  colnames(search.list) <- c(colnames(search.list)[1:11],colnames(Library.phenodata))

  ## attempts to match all peaks against the In House Library ! Try to build a query using *apply functions to
  ## search all of the search list at once; should speed ! things up considerably
  # i = 13 # Used for debugging purposes
  total = nrow(search.list)
  cat("Annotating features against the In House Library.\n\n\n")
  pb = txtProgressBar(min = 0, max = total, style = 3)
  cnt = 1
  for (i in 1:nrow(search.list)) {
    m.min = mz.min[[i]]
    m.max = mz.max[[i]]
    r.min = rt.min[[i]]
    r.max = rt.max[[i]]

    test.list <- tbl(lib_db, paste("Annotated Library", myion.mode, sep = "_")) %>%
      dplyr::filter(between(mz, m.min, m.max)) %>%
      dplyr::filter(between(RT..Min., r.min, r.max)) %>%
      dplyr::collect()

    colnames(test.list)[colnames(test.list)=="mz"] <- "match.mz"
    colnames(test.list)[colnames(test.list)=="adduct"] <- "match.adduct"

    if (nrow(test.list) == 0) {
      search.list[i,"MS.ID"] = paste(bin[[i]], "Unidentified", sep = "_")
    } else {
      # if(nrow(test.list) == 1) {
      #   search.list[i,"MS.ID"] = paste(bin[[i]], "Annotated", sep = "_")
      #   search.list[i,colnames(test.list)] <- test.list[[1]]
      # } else {
      if (nrow(test.list) >= 1) { ##Matches all features against the In House Library
        search.list[i,"MS.ID"] = paste(bin[[i]], "Annotated", sep = "_")
        # search.list[i,"Name"] = glue_collapse(test.list[,"Name"], sep = ";", width = Inf)
        # search.list[i,"Formula"] = glue_collapse(unique(gsub(" ", "", test.list[,"Formula"])), sep = ";", width = Inf)
        # search.list[i,"Molecular.Weight"] = glue_collapse(unique(test.list[,"Molecular.Weight"]), sep = ";", width = Inf)
        # search.list[i,"RT..Min."] = glue_collapse(unique(test.list[,"RT..Min."]), sep = ";", width = Inf)
        # search.list[i,"Annotated.adduct"] = glue_collapse(test.list[,"adduct"], sep = ";", width = Inf)

        total_loop = length(colnames(test.list))

        for (j in 1:total_loop) {
          search.list[i,colnames(test.list)[j]] <- ddply(test.list, ~Ion.Mode, function(x) .mypaste(x,j))[2]
        }


        # for(j in seq(from = 10, to = length(colnames(search.list)), by = 1)) {
        #   k = j - 5
        # search.list[i,j] = glue_collapse(test.list[,k], sep = ";", width = Inf)
        # }
        cnt = cnt + 1
        # }

      }

    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  named.peak.list <- search.list[, 4:ncol(search.list)]
  colnames(named.peak.list)
  temp <- as.data.frame(Peak.list)
  colnames(Peak.list)
  if (nrow(temp) == nrow(named.peak.list)) {
    temp <- cbind(Peak.list, named.peak.list)
  }
  colnames(temp)
  return(temp)
}

#' @title Searches Peak.list for ion mode duplicates
#'
#' @export
#' @description Searches Peak.list with combined ionization mode data tables for duplicate entries
#' @param object used for method dispatch. Can be any object. See usage for details
#' @param Peak.list.pos Positive ionization mode data table
#' @param Peak.list.neg Negative ionization mode data table
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param col.names character vector of column names to include when searching for duplicate entries.
#' Default is to include the ion mode, unique EIC_ID, monomolecular mass and retention time
#' @return data frame containing the original Peak.list with added columns "Duplicate_ID" and "Duplicate_EIC"
search_IonDup <- function(object, Peak.list.pos,Peak.list.neg,search.par,col.names) {
  UseMethod("search_IonDup", object)
}

#' @rdname search_IonDup
#' @export
search_IonDup.mz  <- function(object,Peak.list.pos,Peak.list.neg,search.par,col.names) {

  # Set default values
  if(missing(col.names))
    col.names <- c("Ion Mode", "EIC_ID", "mz", "rt")

  mono.pos <- subset(Peak.list.pos, select = paste(col.names))

  mono.neg <- subset(Peak.list.neg, select = paste(col.names))

  mono.comb <- rbind(mono.pos, mono.neg)

  # get the distance matrices for exact mass and rt
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

#' @rdname search_IonDup
#' @export
search_IonDup.monoMass  <- function(object,Peak.list.pos,Peak.list.neg,search.par,col.names) {

  # Set default values
  if(missing(col.names))
    col.names <- c("Ion Mode", "EIC_ID", "mono_mass", "meanRT")

  mono.pos <- subset(Peak.list.pos, select = paste(col.names))

  mono.neg <- subset(Peak.list.neg, select = paste(col.names))

  mono.comb <- rbind(mono.pos, mono.neg)

  # get the distance matrices for exact mass and rt
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

#' @title Finds background components in Peak.list
#'
#' @export
#' @description Find components within Peak.list that are also present in the process blanks below a user-defined sample:blank ratio
#' @param object used for method dispatch. Can be any object. See usage for details
#' @param Peak.list data table containing sample data
#' @param Solv.list data table containing blank data
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param ion.id character vector specifying identifier in filename designating positive or negative ionization mode or both.  Positive identifier must come first. Default is c('Pos','Neg')
#' @param QC.id character vector specifying identifier in filename designating a Pooled QC sample.  Only the first value will be used.
#' Default is 'Pooled_QC_'
#' @param MB.id character vector specifying identifier in filename designating a Method Blank.
#' Only the first value will be used. Default is '^MB'
#' @param db.id character vector specifying identifiers in database names designating sample \[1\] and blank \[2\] databases.
#' Default is c('Peaklist','Blanks')
#' @param ion.modes a character vector defining the ionization mode.
#' Must be either 'Positive', 'Negative' or both. Default is c('Positive','Negative')
#' @param lib_db RSQLite connection
#' @return nested list a list for each ionization mode, each containing a list of two dataframes: the first contains the intensity matrix for the peaklist with solvent peaks removed, the second contains the intensity matrix for the solvent peaks
find_Background <- function(object, Peak.list, Solv.list, Sample.df, search.par, lib_db,ion.id,QC.id,MB.id,db.id,ion.modes) {
  UseMethod("find_Background", object)
}

#' @rdname  find_Background
#' @export
find_Background.mz <- function(object, Peak.list, Solv.list, Sample.df, search.par, lib_db,ion.id,QC.id,MB.id,db.id,ion.modes) {
  #Set default values
  if(missing(ion.id))
    ion.id <- c("Pos","Neg")
  if(missing(QC.id))
    QC.id <- "Pooled_QC_"
  if(missing(MB.id))
    MB.id <- "^MB"
  if(missing(db.id))
    db.id <- c("Peaklist","Blanks")
  if(missing(ion.modes))
    ion.modes <- c("Positive","Negative")

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
      cat("\nRemoving Background Compounds in Positive mode:\n\n")
      copy_to(lib_db, cur.Solvlist, name = "Pos_list", temporary = FALSE, overwrite = TRUE)
      IHL <- tbl(lib_db, "Pos_list")
    } else {
      if (cur.ion == "Negative") {
        cat("\nRemoving Background Compounds in Negative mode:\n\n")
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
            # temp <- test.list[, c("MS.ID", "Formula", "Name", "Annotated.adduct", "Conf.Level", "FISh.Coverage",
            #                       "isotopes", "adduct")]
            bin[j] = TRUE
            cnt = cnt + 1
          }

        } else {
          if (nrow(test.list) >= 1) {
            if (search.list$sample.mean[j]/max(test.list$mean) <= Solvent.ratio) {
              # temp <- test.list[, c("MS.ID", "Formula", "Name", "Annotated.adduct", "Conf.Level", "FISh.Coverage",
              #                       "isotopes", "adduct")]
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

#' @rdname find_Background
#' @export
find_Background.monoMass <- function(object, Peak.list, Solv.list, Sample.df, search.par, lib_db,ion.id,QC.id,MB.id,db.id,ion.modes) {
  #Set default values
  if(missing(ion.id))
    ion.id <- c("Pos","Neg")
  if(missing(QC.id))
    QC.id <- "Pooled_QC_"
  if(missing(MB.id))
    MB.id <- "^MB"
  if(missing(db.id))
    db.id <- c("Peaklist","Blanks")
  if(missing(ion.modes))
    ion.modes <- c("Positive","Negative")

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

    # Calculate the mean of endo values for each compound
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
      cat("\n\nRemoving Background Compounds in Positive mode:\n\n")
      copy_to(lib_db, cur.Solvlist, name = "Pos_list", temporary = FALSE, overwrite = TRUE)
      IHL <- tbl(lib_db, "Pos_list")
    } else {
      if (cur.ion == "Negative") {
        cat("\n\nRemoving Background Compounds in Negative mode:\n\n")
        copy_to(lib_db, cur.Solvlist, name = "Neg_list", temporary = FALSE, overwrite = TRUE)
        IHL <- tbl(lib_db, "Neg_list")
      } else {
        stop("Ionization mode must be Positive or Negative. It is case sensitive!\n\n", call. = FALSE)
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
            # temp <- test.list[, c("MS.ID", "Formula", "Name", "Annotated.adduct", "Conf.Level", "FISh.Coverage",
            #                       "isotopes", "adduct")]
            bin[j] = TRUE
            cnt = cnt + 1
          }

        } else {
          if (nrow(test.list) >= 1) {
            if (search.list$sample.mean[j]/max(test.list$mean) <= Solvent.ratio) {
              # temp <- test.list[, c("MS.ID", "Formula", "Name", "Annotated.adduct", "Conf.Level", "FISh.Coverage",
              #                       "isotopes", "adduct")]
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
