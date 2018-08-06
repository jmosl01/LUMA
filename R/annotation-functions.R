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
  bin
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
  ## search all of the search list at once; should speed ! things up considerably i = 27 ##Used for debugging
  ## purposes
  total = nrow(search.list)
  pb = txtProgressBar(title = "Annotating Features.", min = 0, max = total, width = NA)
  counter = 1
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
      if (nrow(test.list) >= 1) {
        cat("\n\n\nFeature annotated for match number ", counter, ".\nMatch results below.\n\n\n", sep = "")
        str(test.list)
        search.list$MS.ID[i] = paste(bin[i], "Annotated", sep = "")
        search.list$Name[i] = glue_collapse(test.list$Name, sep = ";", width = Inf, last = " or ")
        search.list$Formula[i] = glue_collapse(unique(gsub(" ", "", test.list$Formula)), sep = ";", width = Inf,
                                                last = " or ")
        search.list$Annotated.adduct[i] = glue_collapse(test.list$adduct, sep = ";", width = Inf, last = " or ")
        search.list$Conf.Level[i] = glue_collapse(test.list$Levels, sep = ";", width = Inf, last = " or ")
        search.list$FISh.Coverage[i] = glue_collapse(test.list$Fish.Coverage, sep = ";", width = Inf,
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

