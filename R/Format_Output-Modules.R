#' @title Formats Peak.list for SIMCA
#'
#' @export
#' @description Formats Peak.list for import into SIMCA
#' @param Peak.list data frame containing combined ion mode peaklist with ion mode duplicates removed.  Alternatively can be retrieved from databases.  Default is NULL
#' @param Sample.df data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Sample.data data frame with phenotype data as columns and a row for each study sample.  First column must be a unique sample identifier with the header 'CT-ID'.  Phenotype columns may vary, but must include two columns called 'Plate Number' and 'Plate Position' for determining run order.
#'    'Plate Number' must be numeric and is equivalent to batch number.  'Plate Position' must be alphanumeric and corresponds to row(alpha) and column(numeric) positions on, e.g. a 96-well plate
#' @param tbl.id character vector of table names to draw from databases.  First value should be table name from positive ionization peak database, second should be table name from negative ionization peak database. Default is NULL
#' @param QC.id character vector specifying identifier in filename designating a Pooled QC sample.  Only the first value will be used.  Default is 'Pooled_QC_'
#' @param ... Arguments to pass parameters to database functions
#' @return NULL testing
#' @importFrom stats setNames
#' @importFrom gtools mixedorder
#' @importFrom plyr rbind.fill
format_simca = function(Peak.list = NULL, Sample.df, Sample.data, tbl.id = NULL, QC.id = "Pooled_QC_", ...) {
    if (missing(Peak.list))
        Peak.list = NULL
    if (missing(tbl.id))
        tbl.id = NULL
    if (is.null(tbl.id) && is.null(Peak.list)) {
        stop("Need to specify tbl.id if using databases to retrieve Peak.list!", call. = FALSE)
    }
    if (is.null(Peak.list)) {
        Peak.list <- read_tbl(tbl.id, peak.db)
    }
    sexes <- unique(paste(Sample.df$Sex, "_", sep = ""))  ## Generate search string for all sexes
    samples <- vector(mode = "character", length = length(colnames(Peak.list)))
    for (i in 1:length(sexes)) {
        rows_loop <- grep(sexes[i], colnames(Peak.list))
        samples[rows_loop] <- sexes[i]
    }
    res <- samples %in% sexes
    sample.peaks <- Peak.list[, res]


  ## Creates a new column for grouping by class based on user input
    groups <- paste(Sample.df$Sex, Sample.df$Class, sep = ";")  ## Generate search string for all classes
    groups <- strsplit(groups, split = ";")
    names(groups) <- paste(Sample.df$Sex, Sample.df$Class, sep = "_")
    group <- vector(mode = "character", length = length(colnames(sample.peaks)))
    for (i in 1:length(groups)) {
        rows_loop <- intersect(grep(groups[[i]][1], colnames(sample.peaks)),
                               grep(groups[[i]][2], colnames(sample.peaks))
                               )
        group[rows_loop] <- names(groups)[i]
    }
    group <- unlist(group)

    # modify sample data to include the user defined exposure class
    Sample.data <- Sample.data[order(Sample.data[,1]), ]
    Sample.data <- setNames(cbind.data.frame(Sample.data[, 1], group, Sample.data[-1]), c(colnames(Sample.data[1]),
        "Exposure Class", colnames(Sample.data[-1])))
    sample.ID <- as.numeric(sub("\\D*(\\d{6}).*", "\\1", colnames(sample.peaks)))  #pulls out the 6-digit numeric sample codes into a vector for matching against the Sample data spreadsheet
    if (all(is.na(sample.ID))) {
        ## Alternatively uses a unique alphanumeric ID provided by user
        sample.ID <- sub("\\D*(\\d{6}).*", "\\1", colnames(sample.peaks))
    }
    sample.peaks <- cbind.data.frame(Peak.list[, "MS.ID"], sample.peaks)
    # Create a new ID for each metabolite and record how many features were combined into that metabolite for
    # metadata
    temp <- as.character(sample.peaks[, 1])

    no.features <- str_count(temp, ";") + 1
    MetID <- gsub("^(.*?);.*", "\\1", temp)
    Mettag <- rep("Unidentified", length = length(MetID), mode = "character")
    Mettag[which(grepl("Annotated", temp, fixed = TRUE))] <- "Annotated"

    my_ion <- Peak.list$Ion.Mode
    my_mass <- round(Peak.list$mono_mass, digits = 5)
    my_rt <- round(Peak.list$meanRT, digits = 2)
    Metbase <- paste(my_ion,my_mass,my_rt,sep = "_")

    MetID <- paste(Metbase, Mettag, sep = "_")
    # End new ID generation

    tsample.peaks <- setNames(data.frame(t(sample.peaks[, -1])), MetID)
    # add columns for sex, sample class, and all phenotype data for each sample
    tsample.peaks <- cbind.data.frame(X = sample.ID, tsample.peaks)
    colnames(tsample.peaks)[1] = colnames(Sample.data)[1]
    tsample.peaks <- merge(Sample.data, tsample.peaks, by = colnames(tsample.peaks)[1])

    mylist <- split(tsample.peaks, f = tsample.peaks$Plate.Number)
    new.list <- lapply(mylist, function(X) X[mixedorder(X[, "Plate.Position"]), ])
    tsample.peaks <- rbind.fill(new.list)

    # make a dataframe of metabolite intensity values for pooled QC samples then transpose so the colnames are the
    # same as metabolite names: Pos_###_Annotated
    Pooled.QC <- QC.id  ## Generate search string for all sexes
    for (i in 1:length(Pooled.QC)) {
        rows_loop <- grep(Pooled.QC[i], colnames(Peak.list))
        samples[rows_loop] <- Pooled.QC[i]
    }
    QC <- samples %in% Pooled.QC
    QC.peaks <- Peak.list[, QC]
    QC.peaks <- cbind.data.frame(Peak.list[, "MS.ID"], QC.peaks)
    tQC.peaks <- setNames(data.frame(t(QC.peaks[, -1])), MetID)
    # add columns for sex, sample class, and all phenotype data for Pooled QCs, and fill with NA
    dummy.data <- setNames(data.frame(matrix(ncol = length(Sample.data), nrow = nrow(tQC.peaks))), colnames(Sample.data))
    tQC.peaks <- cbind.data.frame(dummy.data, tQC.peaks)


    # make a dataframe of all metadata columns then transpose so the colnames are the same as the metabolite
    # names: Pos_###_Annotated
    log.list <- list()
    log.list[[1]] <- res
    log.list[[2]] <- QC
    meta.peaks <- Peak.list[, !Reduce("|", log.list)]
    temp <- str_count(meta.peaks$Duplicate_EIC, ";")
    Ion.duplicate <- vector(mode = "logical", length = length(temp))
    Ion.duplicate[which(temp > 0)] = TRUE  #add flag for the metabolites detected in both ion modes
    endo <- meta.peaks[, "Endogenous_flag"]
    meta.peaks$Endogenous_flag <- as.logical(meta.peaks$Endogenous_flag)  #convert endogenous flag to a logical
    meta.peaks <- cbind.data.frame(meta.peaks[, 1:11], no.features, meta.peaks[, 12:ncol(meta.peaks)], Ion.duplicate)
    tmeta.peaks <- setNames(data.frame(t(meta.peaks[, -1])), MetID)
    dummy.data <- setNames(data.frame(matrix(ncol = length(Sample.data), nrow = nrow(tmeta.peaks))), colnames(Sample.data))
    tmeta.peaks <- cbind.data.frame(dummy.data, tmeta.peaks)
    # pheno.peaks <- Peak.list[,-grep('[MF]_',colnames(Peak.list))] colnames(pheno.peaks)

    # combine the three dataframes, first sample data, then QC data, lastly metadata
    SIMCA.data <- rbind.data.frame(tsample.peaks, tQC.peaks, tmeta.peaks)
    SIMCA.data <- cbind.data.frame(X = rownames(SIMCA.data), SIMCA.data)
    SIMCA.data[, 1] <- gsub("^X", "\\1", SIMCA.data[, 1])
    SIMCA.data <- setNames(cbind.data.frame(SIMCA.data, row.names = NULL), c("", colnames(SIMCA.data)[-1]))



}
