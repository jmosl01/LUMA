#' @title Calculate correlation matrices for metabolite groups
#'
#' @description Calculates the correlation matrices for metabolite groups based on the best feature within the group that belongs to the primary metabolite
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Peak.list a data frame from CAMERA that has been parsed.  Should contain all output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac and CAMERA.parser.
#' @param get.mg numerical vector of metabolite groups that have more than one feature
#' @param BLANK a logical indicating whether blanks are being evaluated
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @return a table of class 'tbl_df',tbl' or 'data.frame' with variables as columns.  Has all of the columns as the original data frame with one additional column 'Correlation.stat'
#' @importFrom stats cor dist
#' @importFrom utils setTxtProgressBar
#' @importFrom flashClust flashClust
calc_corrstat = function(Sample.df, Peak.list, get.mg, BLANK, ion.mode) {

    ## Error check
    Peak.list.pspec <- Peak.list[which(Peak.list$metabolite_group %in% get.mg), ]
    if (sum(ifelse(colnames(Peak.list) == colnames(Peak.list.pspec), 0, 1)) != 0) {
        stop("Something is wrong with your metabolite groups. Check out from Matlab script!")
    }
    nrow(Peak.list)
    nrow(Peak.list.pspec)

    if (BLANK == FALSE) {
        sexes <- unique(paste(Sample.df$Sex, "_", sep = ""))
        sample.cols <- grep(paste("Pooled_QC_", paste(strsplit(sexes, "(?<=.[_])", perl = TRUE), collapse = "|"),
            sep = "|"), colnames(Peak.list))
    } else {
        if (ion.mode == "Positive" && BLANK == TRUE) {
            sample.cols <- grep("_Pos", colnames(Peak.list), ignore.case = TRUE)
        } else {
            if (ion.mode == "Negative" && BLANK == TRUE) {
                sample.cols <- grep("_Neg", colnames(Peak.list), ignore.case = TRUE)
            }
        }
    }

    corr.group <- as.matrix(unique(Peak.list.pspec[, "metabolite_group"]))
    corr.group <- sort(corr.group)
    length(corr.group)
    ## Error Check
    if (sum(ifelse(corr.group != get.mg, 1, 0)) != 0) {
        stop("Something is wrong with your metabolite groups. Check out from Matlab script!")
    }
    corr.stat = vector(mode = "numeric", length = nrow(Peak.list.pspec))
    corr.df <- data.frame(corr.stat = corr.stat, row.names = as.character(Peak.list.pspec$EIC_ID))
    rownames(corr.df)
    colnames(corr.df)
    corr.df
    total = length(get.mg)
    # i = 17 For debugging purposes
    cat("Generating Correlation Matrices.\n\n\n")
    pb = txtProgressBar(min = 0, max = total, width = NA)
    for (i in 1:length(corr.group)) {
        my.df <- Peak.list[which(Peak.list$metabolite_group %in% get.mg[i]), ]
        colnames(my.df)
        if (nrow(my.df) == 1) {
        } else {
            if (nrow(my.df) == 2) {
                my.mat <- as.matrix(my.df[, sample.cols])
                test.mat <- as.matrix(t(my.mat))
                dimnames(test.mat) <- list(colnames(my.df[, sample.cols]), my.df$EIC_ID)
                colnames(test.mat)
                res <- cor(test.mat)
                d <- dist(res)
                sim.by.hclust <- flashClust(d)
                attributes(d)
                labels(d)
                # plot(sim.by.hclust)
                attributes(sim.by.hclust)
                sim.by.hclust$merge
                sim.by.hclust$labels
                labels(d)[sim.by.hclust$order]
                corr.stat <- res[which(rowSums(res) %in% max(rowSums(res))), ]
                new.df <- as.data.frame(corr.stat, row.names = names(corr.stat))
                j <- rownames(new.df)
                corr.df[j, ] <- new.df
                corr.df[j, ]
                corr.df
            } else {
                my.mat <- as.matrix(my.df[, sample.cols])
                test.mat <- as.matrix(t(my.mat))
                dimnames(test.mat) <- list(colnames(my.df[, sample.cols]), my.df$EIC_ID)
                colnames(test.mat)
                res <- cor(test.mat)
                d <- dist(res)
                sim.by.hclust <- flashClust(d)
                attributes(d)
                labels(d)
                # plot(sim.by.hclust)
                attributes(sim.by.hclust)
                sim.by.hclust$merge
                sim.by.hclust$labels
                labels(d)[sim.by.hclust$order]
                corr.stat <- res[which(rowSums(res) %in% max(rowSums(res))), ]
                new.df <- as.data.frame(corr.stat, row.names = names(corr.stat))
                j <- rownames(new.df)
                corr.df[j, ] <- new.df
                corr.df[j, ]
                corr.df
                # ordered.hclust <- reorder(sim.by.hclust, wts = my.df$monoisotopic_flg, agglo.FUN = 'mean')
                # ordered.hclust$value ordered.hclust$labels

            }
        }
        setTxtProgressBar(pb, i)
    }
    Peak.list.pspec[, "Correlation.stat"] <- corr.df
    close(pb)
    return(Peak.list.pspec)
}

#' @title Calculate Minfrac for sample classes
#'
#' @export
#' @description Calculates the fraction of samples within each class that contains a given feature.  Only the minimum fraction across all sample classes is returned.
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param xset4 an xcms object that has had peak picking, retention time alignment, peak grouping, and imputing missing values performed
#' @param BLANK a logical indicating whether blanks are being evaluated
#' @param Peak.list a table of class 'tbl_df',tbl' or 'data.frame' with variables as columns.  Should contain all output columns from XCMS and CAMERA, and additional columns from IHL.search.
#' @return data frame containing the original table with one additional column 'Minfrac' at the end, followed by the CAMERA columns 'isotopes','adduct','pcgroup'
#' @importFrom xcms peakTable
#' @importFrom utils read.table write.table str head
#' @importFrom stats variable.names
calc_minfrac = function(Sample.df, xset4, BLANK, Peak.list) {
    peakSN <- peakTable(xset4, filebase=NULL, value="sn") #writes the SN peak table to file
    SN.list <- data.frame(X = rownames(peakSN),peakSN)

    if (BLANK == TRUE) {
    } else {
        sexes <- unique(paste(Sample.df$Sex, "_", sep = ""))
        samples <- vector(mode = "character", length = length(colnames(SN.list)))
        for (i in 1:length(sexes)) {
            rows_loop <- grep(sexes[i], colnames(SN.list))
            samples[rows_loop] <- sexes[i]
        }
        sum.range.list <- SN.list[, samples %in% sexes]
        t.list<-as.data.frame(t(sum.range.list)) #transposes the SN list

        ## Grabs the unique sample ID number for each sample, if present in the filename
        bin <- variable.names(t.list, full = TRUE)
        colnames(sum.range.list)
        sample.ID <- sub("\\D*(\\d{6}).*", "\\1", colnames(sum.range.list))
        colnames(sum.range.list)

        ## Creates a new column for grouping by class based on user input
        groups <- paste(Sample.df$Sex, Sample.df$Class, sep = "_")
        group <- vector(mode = "character", length = length(colnames(sum.range.list)))
        for (i in 1:length(groups)) {
            rows_loop <- grep(groups[i], colnames(sum.range.list))
            group[rows_loop] <- groups[i]
        }
        class.n <- unique(Sample.df$n)

        if (length(class.n) > 1) {
            t.list <- cbind(group, t.list)  #combines the two dataframes together
            sn.list <- split(t.list, as.factor(group))  #Splits the data frame by a factor
            sn.list.split <- split(sn.list, sapply(sn.list, function(x) nrow(x)))

            # Create matrix to contain all of the minfrac values for each group from each list in sn.list
            min.frac.fulltable <- matrix(nrow = length(bin), ncol = sum(sapply(sn.list, function(x) length(x))))
            # j = 1 This and following line used for debugging purposes k = 1
            for (j in 1:length(sn.list)) {
                na_count <- sapply(as.data.frame(sn.list[j]), function(x, y) sum(length(which(is.na(x)))), y = group)  #counts the number of na values per feature across each group
                na_count <- data.frame(na_count)  #puts the list in a data frame
                for (i in 1:length(groups)) {
                  rows_loop <- grep(groups[i], rownames(na_count))
                  group[rows_loop] <- groups[i]
                }
                na_count <- cbind(group, na_count)  #combines the grouping variable with the na count values
                na.count.list <- split(na_count, as.factor(group))  #Splits the data frame into lists of data frames by the grouping variable
                str(na.count.list)

                all_count <- sapply(as.data.frame(sn.list[j]), function(x, y) length(x), y = group)  #counts the number of na values per feature across each group
                all_count <- data.frame(all_count)  #puts the list in a data frame
                str(all_count)
                for (i in 1:length(groups)) {
                  rows_loop <- grep(groups[i], rownames(all_count))
                  group[rows_loop] <- groups[i]
                }
                all_count <- cbind(group, all_count)  #combines the grouping variable with the na count values
                all.count.list <- split(all_count, as.factor(group))  #Splits the data frame into lists of data frames by the grouping variable
                str(all.count.list)

                ## The next section creates a new data matrix with the number of rows and columns equal to the number of
                ## features and groups, respectively.  The new table is called 'minfrac.fulltable_loop'
                minfrac.fulltable_loop <- matrix(nrow = length(bin), ncol = length(all.count.list))
                colnames(minfrac.fulltable_loop) <- c(names(na.count.list))  ###Names the columns for the new table as indicated
                rownames(minfrac.fulltable_loop) <- (bin)  ###Pulls out and names the table rows based on the metabolite label (I do this as a QA to make sure the loop is correctly linking the ANOVA output and metbaolite IDs

                ## this section calculates the minfrac values for each exposure class in a for loop and populates the
                ## minfrac.fulltable_loop with these values
                group.names <- names(na.count.list)
                for (i in names(na.count.list)) {
                  na_loop <- na.count.list[i]
                  na_loop <- do.call("rbind", na_loop)
                  all_loop <- all.count.list[i]
                  all_loop <- do.call("rbind", all_loop)
                  min_frac_loop <- rapply(as.data.frame(na_loop$na_count), function(x, y) x/y, y = all_loop$all_count)  ##calculates the min frac values for each feature by exposure class
                  min_frac_loop <- data.frame(min_frac_loop)  ##changes the format of the min frac results into a data frame
                  colnames(min_frac_loop) <- paste(i, "", "min_frac")
                  min_frac_loop <- min_frac_loop[-1, ]  ##Removes the value corresponding to the grouping variable
                  minfrac.fulltable_loop[, i] = min_frac_loop  ##Pastes the min frac results from an exposure class into the next col of the min frac full table.
                }
                ## Code to calculate min_frac for one of the exposure classes; used to QA the loop output
                i = names(na.count.list[1])
                group.names <- names(na.count.list)
                test.na <- na.count.list[i]
                test.na <- do.call("rbind", test.na)
                test.all <- all.count.list[i]
                test.all <- do.call("rbind", test.all)
                min_frac <- rapply(as.data.frame(test.na$na_count), function(x, y) x/y, y = as.data.frame(test.all$all_count))  #Need to find a way to get rid of the subset operator $
                min_frac <- data.frame(min_frac)
                str(min_frac)
                colnames(min_frac) <- paste(i, "", "min_frac")
                min_frac <- min_frac[-1, ]
                head(min_frac)

                # Bring loop-based minfrac values out of the loop by placing them in the minfrac master table
                min.frac.fulltable[, k:sum(k - 1, ncol(minfrac.fulltable_loop))] <- minfrac.fulltable_loop

                # Update the counter to keep track of where to copy the results in the master table
                k <- k + ncol(minfrac.fulltable_loop)
            }
        } else {
            if (length(class.n) == 1) {
                t.list <- cbind(group, t.list)  #combines the two dataframes together
                sn.list <- split(t.list, as.factor(group))  #Splits the data frame by a factor

                na_count <- sapply(as.data.frame(sn.list), function(x, y) sum(length(which(is.na(x)))), y = group)  #counts the number of na values per feature across each group
                na_count <- data.frame(na_count)  #puts the list in a data frame
                # group<-unlist(strsplit(gsub('(([MF])_([RN])S_T([0482]))|.', '\\1', rownames(na_count)), '\\s+'))
                # #generates a character vector of grouping variables for each na count value
                for (i in 1:length(groups)) {
                  rows_loop <- grep(groups[i], rownames(na_count))
                  group[rows_loop] <- groups[i]
                }
                na_count <- cbind(group, na_count)  #combines the grouping variable with the na count values
                na.count.list <- split(na_count, as.factor(group))  #Splits the data frame into lists of data frames by the grouping variable
                all_count <- sapply(as.data.frame(sn.list), function(x, y) length(x), y = group)  #counts the number of na values per feature across each group
                all_count <- data.frame(all_count)  #puts the list in a data frame
                for (i in 1:length(groups)) {
                  rows_loop <- grep(groups[i], rownames(all_count))
                  group[rows_loop] <- groups[i]
                }
                all_count <- cbind(group, all_count)  #combines the grouping variable with the na count values
                all.count.list <- split(all_count, as.factor(group))  #Splits the data frame into lists of data frames by the grouping variable

                ## The next section creates a new data matrix with the number of rows and columns equal to the number of
                ## features and groups, respectively.  The new table is called 'min.frac.fulltable'
                min.frac.fulltable <- matrix(nrow = length(bin), ncol = length(all.count.list))
                colnames(min.frac.fulltable) <- c(names(na.count.list))  ###Names the columns for the new table as indicated
                rownames(min.frac.fulltable) <- (bin)  ###Pulls out and names the table rows based on the metabolite label (I do this as a QA to make sure the loop is correctly linking the ANOVA output and metbaolite IDs
                head(min.frac.fulltable)  ##A double-check to make sure the new table was created properly

                # this section calculates the minfrac values for each exposure class in a for loop and populates the
                # min.frac.fulltable with these values
                group.names <- names(na.count.list)
                for (i in names(na.count.list)) {
                  na_loop <- na.count.list[i]
                  na_loop <- do.call("rbind", na_loop)
                  all_loop <- all.count.list[i]
                  all_loop <- do.call("rbind", all_loop)
                  min_frac_loop <- rapply(as.data.frame(na_loop$na_count), function(x, y) x/y, y = all_loop$all_count)  ##calculates the min frac values for each feature by exposure class
                  min_frac_loop <- data.frame(min_frac_loop)  ##changes the format of the min frac results into a data frame
                  colnames(min_frac_loop) <- paste(i, "", "min_frac")
                  min_frac_loop <- min_frac_loop[-1, ]  ##Removes the value corresponding to the grouping variable
                  min.frac.fulltable[, i] = min_frac_loop  ##Pastes the min frac results from an exposure class into the next row of the min frac full table.
                }
                # Code to calculate min_frac for one of the exposure classes; used to QA the loop output
                i = names(na.count.list[2])
                group.names <- names(na.count.list)
                test.na <- na.count.list[i]
                test.na <- do.call("rbind", test.na)
                test.all <- all.count.list[i]
                test.all <- do.call("rbind", test.all)
                min_frac <- rapply(as.data.frame(test.na$na_count), function(x, y) x/y, y = as.data.frame(test.all$all_count))  #Need to find a way to get rid of the subset operator $
                min_frac <- data.frame(min_frac)
                colnames(min_frac) <- paste(i, "", "min_frac")
                min_frac <- min_frac[-1, ]
                head(min_frac)

            } else print("Please include the number of samples in each class in the Sample Class file!")
        }
        ## Code to calculate MinFrac as the minimum fraction of samples with a feature across all exposure classes and
        ## sexes
        MinFrac.table <- matrix(nrow = length(bin), ncol = 1)  #This creates a new table to fill with one MinFrac value for each feature
        for (i in 1:length(bin)) {
            minfrac_loop <- 1 - min(min.frac.fulltable[i, ], na.rm = TRUE)
            MinFrac.table[i, ] <- minfrac_loop
        }
    }

    CAMERA.list <- Peak.list[, c("isotopes", "adduct", "pcgroup")]  #Extracts all of the CAMERA columns to a separate dataframe
    drops <- c("isotopes", "adduct", "pcgroup")
    stripped.list <- Peak.list[, !(names(Peak.list) %in% drops)]  #Removes the CAMERA columns
    stripped.list[, "MinFrac"] <- MinFrac.table  #Attaches the MInFrac column to the stripped down peak.list
    new.peak.list <- cbind(stripped.list, CAMERA.list)
    raw <- new.peak.list
    return(raw)
}

#' @title Combine phenotype data
#'
#' @export
#' @description Combine phenotype data for each metabolite group with summed intensity values
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Peak.list data frame. Must have Correlation.stat column.  Should contain output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac, CAMERA.parser and EIC.plotter functions.
#' @param Summed.list data frame containing metabolite group as first column and the rest summed intensities for each sample
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param BLANK a logical indicating whether blanks are being evaluated
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @return data frame Peak.list.summed with the munged phenotype columns up front followed by QC and sample columns
#' @importFrom plyr ddply
#' @importFrom dplyr '%>%' mutate_if summarise bind_cols
#' @importFrom utils str txtProgressBar setTxtProgressBar
#' @importFrom stringr str_count
combine_phenodata = function(Sample.df, Peak.list, Summed.list, search.par, BLANK, ion.mode) {

  mylist <- gen_res(ion.mode,search.par,Peak.list,Sample.df,BLANK)
  Peaklist_corstat <- mylist[[1]]
  res <- mylist[[2]]

    # Creates a data frame with only the pheno data columns combined for each metabolite group
    pheno.list <- Peaklist_corstat[sapply(res, function(x) length(x) < 1)]  #Extracts all of the pheno columns for combining by metabolite group
    pheno.list[,"metabolite_group"] <- Peaklist_corstat[,"metabolite_group"]
    pheno.list <- pheno.list %>% mutate_if(is.factor, as.character)
    attributes(pheno.list)
    str(pheno.list)
    mylist <- sort(unique(Peaklist_corstat$metabolite_group))
    new.pheno.list <- as.data.frame(mylist)
    # i = 11 ##for debugging purposes
    total = length(colnames(pheno.list))
    cat("Combining metadata for isotopes and adducts.\n\n\n")
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    mypaste = function(pheno.list, i) {
        summarise(pheno.list, X = paste0(pheno.list[, i], collapse = ";"))
    }
    for (i in 1:total) {
        if (is.character(pheno.list[, i])) {
            # works new.pheno.list[,colnames(pheno.list)[i]] <- ddply(pheno.list, .(metabolite_group), summarise,
            # X=paste0(MS.ID, collapse = ';'))[2]
            new.pheno.list[, colnames(pheno.list)[i]] <- ddply(pheno.list, ~metabolite_group, function(x) mypaste(x,
                i))[2]
        } else if (!is.na(is.whole(pheno.list[, i]))) {
            new.pheno.list[, colnames(pheno.list)[i]] <- ddply(pheno.list, ~metabolite_group, function(x) mypaste(x,
                i))[2]
        }
        setTxtProgressBar(pb, i)
    }

    # Create original pheno list to be replaced once combined with the summed results
    Iso.only.list <- subset(pheno.list, pheno.list$Correlation.stat == 1)  #extracts only the prime feature for each metabolite group across all columns

    # Replace pheno columns with combined phenotype data
    new.pheno.list <- new.pheno.list[, -which(colnames(new.pheno.list) %in% "metabolite_group")]
    nm <- intersect(colnames(new.pheno.list), colnames(Iso.only.list))
    colnames(new.pheno.list) <- c("metabolite_group", colnames(new.pheno.list)[-1])

    new.pheno.list[["MS.ID"]][match(Iso.only.list$metabolite_group, new.pheno.list$metabolite_group)]
    Iso.only.list[nm] <- lapply(nm, function(x) new.pheno.list[[x]][match(Iso.only.list$metabolite_group, new.pheno.list$metabolite_group)])

    # Combine munged pheno columns with summed data
    df1 <- Iso.only.list[order(Iso.only.list$metabolite_group), ]
    df2 <- Summed.list[order(Summed.list$metabolite_group), ]
    Peak.list.summed <- bind_cols(df1, df2)
    colnames(Peak.list.summed)
    temp <- as.character(Peak.list.summed$MS.ID)
    str(temp)
    no.features <- str_count(temp, ";") + 1
    drop.singletons = !search.par[1, "keep.singletons"]
    if (drop.singletons & BLANK == FALSE) {
        drops <- grep(1, no.features, value = FALSE)
        temp <- Peak.list.summed[drops, ]
        Peak.list.summed <- Peak.list.summed[-drops, ]
    }
    return(Peak.list.summed)
}


#' @title Sum Features into Metabolites
#'
#' @export
#' @description Sums all features belonging to the same metabolite into a single intensity value per metabolite group per sample
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Peak.list data frame. Must have Correlation.stat column.  Should contain output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac, CAMERA.parser and EIC.plotter functions.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param BLANK a logical indicating whether blanks are being evaluated
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @return data table data frame sum.range.list with the first column containing metabolite group and the rest containing sample and QC columns
#' @importFrom data.table as.data.table
sum_features = function(Sample.df, Peak.list, search.par, BLANK, ion.mode) {
    mylist <- gen_res(ion.mode,search.par,Peak.list,Sample.df,BLANK)
    Peaklist_corstat <- mylist[[1]]
    res <- mylist[[2]]
    sum.range.list <- Peaklist_corstat[sapply(res, function(x) length(x) > 0)]  #Extracts all of the sample columns for summing by metabolite group
    DT <- as.data.table(sum.range.list)  #Puts the summing columns in data table format
    sum.range.list <- DT[, lapply(.SD, sum), by = metabolite_group]  #sums features within each sample by metabolite group
    return(sum.range.list)
}

is.whole <- function(a) {
  (is.numeric(a) && floor(a) == a) || (is.complex(a) && floor(Re(a)) == Re(a) && floor(Im(a)) == Im(a))
}

gen_res <- function(ion.mode,search.par,Peak.list,Sample.df,BLANK) {
  if (ion.mode == "Positive") {
    cor.stat <- as.numeric(search.par[1, "Corr.stat.pos"])
  } else {
    if (ion.mode == "Negative") {
      cor.stat <- as.numeric(search.par[1, "Corr.stat.neg"])
    }
  }
  Peaklist_corstat <- Peak.list[which(Peak.list$Correlation.stat >= cor.stat), ]
  if (BLANK == FALSE) {
    sexes <- unique(paste(Sample.df$Sex, "_", sep = ""))
    paste(sexes, sep = "|")
    res <- lapply(colnames(Peaklist_corstat), function(ch) unique(grep(paste("Pooled_QC_", paste(strsplit(sexes,
                                                                                                          "(?<=.[_])", perl = TRUE), collapse = "|"), "metabolite_group", sep = "|"), ch)))  #Flags all of the sample columns and the metabolite group data
  } else {
    if (ion.mode == "Positive" && BLANK == TRUE) {
      res <- lapply(colnames(Peaklist_corstat), function(ch) unique(grep("_Pos|metabolite_group", ch, ignore.case = TRUE)))  #Flags all of the sample columns and the metabolite group data
    } else {
      if (ion.mode == "Negative" && BLANK == TRUE) {
        res <- lapply(colnames(Peaklist_corstat), function(ch) unique(grep("_Neg|metabolite_group", ch,
                                                                           ignore.case = TRUE)))  #Flags all of the sample columns and the metabolite group data
      }
    }
  }
  return(list(Peaklist_corstat,res))
}
