#' @title Searches metabolite data against In House Library (IHL)
#'
#' @export
#' @description Searches metabolite data against user-provided IHL. Useful for comparing CAMERA isotope and ion adduct annotations against known information from user's experience.
#' @param Peak.list a table of class 'tbl_dbi', 'tbl_sql', 'tbl_lazy', or 'tbl' with samples as columns.  Should contain all output columns from XCMS and CAMERA, both metadata and sample data. Retention times must be in min.
#' @param Annotated.library a data frame with annotated metabolite entries. Must contain columns called 'Name', 'Formula', Molecular.Weight' and 'RT..Min.'.  Can contain additional info as separate columns.
#' @param rules a data frame containing the rule list used by CAMERA to annotate ion adducts and fragments.  Must contain the columns 'name','nmol','charge','massdiff','oidscore','quasi','ips'.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @param lib_db RSQLite connection
#' @param mycols a character vector containing all the names to include in the search output
#' @param ... parameters to pass to database functions
#' @return data frame containing the original table with added columns 'Name','MS.ID','Formula','Annotated.adduct' and any additional info columns from Annotated.library
#' @importFrom dplyr '%>%' select copy_to tbl between
#' @importFrom glue collapse
#' @importFrom utils str txtProgressBar
search_IHL = function(Peak.list, Annotated.library, rules, search.par, ion.mode, lib_db, mycols, ...) {
  if (ion.mode == "Positive") {
    Ion.Mode <- "Pos"
    bin <- paste(Ion.Mode, "_", search.list$EIC_ID, "_", sep = "")

  } else {
    if (ion.mode == "Negative") {
      Ion.Mode <- "Neg"
      bin <- paste(Ion.Mode, "_", search.list$EIC_ID, "_", sep = "")

    } else {
      stop("You must include the ionization mode!")
    }
  }

  search.list <- get_features(myname,
                              peak.db,
                              asdf=TRUE)
  IHL <- gen_adducts(Annotated.library,
                     rules,
                     ion.mode)
  write_tbl(mytbl = IHL,
            peak.db = lib_db,
            myname = "Annotated Library")
  new.search.list <- calc_ranges(search.list,
                                 search.par,
                                 mycols)
  search.list <- new.search.list
  ## attempts to match all peaks against the In House Library
  ## Try to build a query using *apply functions to
  ## search all of the search list at once; should speed things up considerably
  i = 27 ##Used for debugging purposes
  total = nrow(search.list)
  pb = txtProgressBar(title = "Annotating Features.", min = 0, max = total, width = NA)
  counter = 1
  for (i in 1:nrow(search.list)) {
    mz.min = search.list$mz.min[i]
    mz.max = search.list$mz.max[i]
    rt.min = search.list$rt.min[i]
    rt.max = search.list$rt.max[i]
    test.list <- IHL %>%
                      dplyr::filter(between(mz, mz.min, mz.max)) %>%
                            dplyr::filter(between(RT..min., rt.min, rt.max)) %>%
                                  dplyr::collect()
    if (nrow(test.list) == 0) {
      search.list$MS.ID[i] = paste(bin[i], "Unidentified", sep = "")
    } else {
      if (nrow(test.list) >= 1) {
        cat("\n\n\nFeature annotated for match number ", counter, ".\nMatch results below.\n\n\n", sep = "")
        str(test.list)
        search.list$MS.ID[i] = paste(bin[i], "Annotated", sep = "")
        search.list$Name[i] = glue::collapse(test.list$Name, sep = ";", width = Inf, last = " or ")
        search.list$Formula[i] = glue::collapse(unique(gsub(" ", "", test.list$Formula)), sep = ";", width = Inf,
                                                last = " or ")
        search.list$Annotated.adduct[i] = glue::collapse(test.list$adduct, sep = ";", width = Inf, last = " or ")
        search.list$Conf.Level[i] = glue::collapse(test.list$Levels, sep = ";", width = Inf, last = " or ")
        search.list$FISh.Coverage[i] = glue::collapse(test.list$Fish.Coverage, sep = ";", width = Inf,
                                                      last = " or ")
        counter = counter + 1
      }
    }

    setTxtProgressBar(pb, i, title = paste("Annotating Features: ", round(i/total * 100, 0), "% done", sep = ""))
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

#' @title Calculates search ranges
#'
#' @description Calculates the min and max ranges for searching against the In House Library (IHL)
#' @param search.list tbl containing the list of mz and rt values to search against the IHL. Must contain only three columns, "EIC_ID", "mz", and "rt"
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param mycols a character vector containing all the names to include in the search output
#' @return NULL testing
calc_ranges = function(search.list,search.par,mycols) {
  if(ncol(search.list) != 3) {
    ##error check
    stop("search.list must contain only three columns, \"EIC_ID\", \"mz\", and \"rt\"!")
  }
  mylist <- search.list
  d.mz <- search.list$mz * as.numeric(search.par[1, "ppm"])/10^6
  d.rt <- as.numeric(search.par[1, "rt"])

  # calculates the min and max range values for searching the in house library
  mylist$mz.min <- search.list$mz - d.mz
  mylist$mz.max <- search.list$mz + d.mz
  mylist$rt.min <- search.list$rt - d.rt
  mylist$rt.max <- search.list$rt + d.rt
  n <- length(mycols) + 7
  mylist[,8:n] <- NA
  colnames(mylist) <- c(colnames(search.list),"mz.min","mz.max","rt.min","rt.max",mycols)
  return(mylist)
}

#' @title Generates adducts for in house library
#'
#' @description Generates a list of adducts for all metabolites in the supplied annotation library based on the rules input
#' @param Annotated.library a data frame with annotated metabolite entries. Must contain columns called 'Name', 'Formula', Molecular.Weight' and 'RT..Min.'.  Can contain additional info as separate columns.
#' @param rules a data frame containing the rule list used by CAMERA to annotate ion adducts and fragments.  Must contain the columns 'name','nmol','charge','massdiff','oidscore','quasi','ips'.
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @return NULL testing
gen_adducts = function(Annotated.library,rules,ion.mode) {
  ## Creates full adduct list for all compounds in the Annotated library
  IHL <- Annotated.library[rep(seq_len(nrow(Annotated.library)), each = nrow(rules)), ]
  x <- rules$nmol
  IHL.temp <- sweep(IHL, 1, x, "*")
  x <- rules$massdiff
  if (ion.mode == "Positive") {
    IHL.temp <- sweep(IHL.temp, 1, x, "+")

  } else {
    if (ion.mode == "Negative") {
      IHL.temp <- sweep(IHL.temp, 1, x, "+")
      IHL.temp <- sweep(IHL.temp, 1, -1, "*")

    } else {
      stop("You must include the ionization mode! Try\n ion.mode = c(\"Positive\",\"Negative\"")
    }
  }
  x <- rules$charge
  IHL.adduct.data <- sweep(IHL.temp, 1, x, "/")
  IHL[, "mz"] <- IHL.adduct.data$Molecular.Weight
  IHL[, "adduct"] <- rules$name
  return(IHL)
}
