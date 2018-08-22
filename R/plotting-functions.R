#' @title Loop-based EIC Plotter
#'
#' @export
#' @description Plots EICs and Pseudospectra from (parsed) metabolite groups using plotEICs-methods and plotPsSpectrum-methods from CAMERA. Also plots correlation matrices and clustered dendrograms for each (parsed) metabolite group.
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Peak.list a data frame from CAMERA that has been parsed.  Should contain all output columns from XCMS and CAMERA, and additional columns from IHL.search, Calc.MinFrac and CAMERA.parser.
#' @param center numeric value indicating which sample to pick for plotting purposes
#' @param BLANK a logical indicating whether blanks are being evaluated. Default is FALSE
#' @param gen.plots a logical indicating whether to create plots for metabolite groups.  Default is FALSE
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @param anposGa xsannotate object with annotated isotopes and ion adducts and fragments
#' @param file.base character string used to name graphical output.  Will be appended with '_CorrPlots.pdf'
#' @param QC.id character identifier for pooled QC samples. Default is 'Pooled_QC'
#' @param maxlabel numeric How many m/z labels to print
#' @return List of length 2.  1st element is a data frame with all columns as the original data frame with one additional column 'Correlation.stat'.  2nd element is a data frame specifically used to validate CAMERA results.
#' @importFrom  CAMERA plotEICs
#' @importFrom CAMERA plotPsSpectrum
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom graphics layout par plot text
#' @importFrom Hmisc rcorr
#' @importFrom corrplot corrplot corrMatOrder
#' @importFrom dplyr group_by
plot_metgroup = function(anposGa, Sample.df, Peak.list, center, BLANK, gen.plots, ion.mode, file.base, QC.id, maxlabel) {
    if (missing(BLANK))
        BLANK = FALSE
    if (missing(gen.plots))
        gen.plots = FALSE
    if (missing(QC.id))
        QC.id = "Pooled_QC"
    if (missing(maxlabel))
        maxlabel = 10


    # Change the psspectra list in the xsAnnotate object to be the unique list of metabolite_groups
    X <- split(Peak.list$EIC_ID, as.numeric(Peak.list$metabolite_group))
    names(X) <- sort(unique(Peak.list$metabolite_group))
    new_anposGa <- anposGa
    new_anposGa@pspectra <- X

    ## Need to set a numeric value in this slot, else plotEICs will give error
    new_anposGa@sample <- c(1:nrow(new_anposGa@xcmsSet@phenoData))

    # Gets the EICs and Psspectra for each metabolite group that has more than one feature
    pspec.length <- sapply(new_anposGa@pspectra, function(x) length(x))
    get.mg <- which(pspec.length > 1)
    names(get.mg) <- NULL

    Peak.list.pspec <- calc_corrstat(Sample.df, Peak.list, get.mg, BLANK, ion.mode)
    validate.list <- validate_metgroup(Peak.list.pspec)

    Peak.list.new <- list(Peak.list.pspec, validate.list)
    # ! Very important!! Needs a sample index for each new psspectrum; without it, it will only find files ! for
    # the first n spectra, where n is the original number of psspectra. ! Find out how CAMERA performs automatic
    # selection, then replicate here for the new psspectra ! For now, I am using the center sample for every
    # psspectra
    new_anposGa@psSamples <- rep(center, length(X))  #Tells CAMERA to select all psspectra from the center file
    # new_anposGa@psSamples <- rep(-1, length(X)) #Experimental new_anposGa@sample <- as.numeric(NA) #doesn't work
    # for some reason

    # New code to plot correlation matrix plots for each metabolite group containing more than one feature
    if (gen.plots && !BLANK) {
        sexes <- unique(paste(Sample.df$Sex, "_", sep = ""))
        sample.cols <- grep(paste(QC.id, paste(strsplit(sexes, "(?<=.[_])", perl = TRUE), collapse = "|"), sep = "|"),
            colnames(Peak.list))
        pdf(file = paste(file.base, "CorrPlots.pdf", sep = "_"))

        par(mar = c(5, 4, 4, 2) + 0.1)


        total = length(get.mg)
        i = 8  #For Debugging purposes
        for (i in 1:length(get.mg)) {
            layout(matrix(c(2, 1, 1, 4, 3, 3, 4, 3, 3), nrow = 3, ncol = 3, byrow = TRUE))
            ## Code to plot Upper Correlation matrix
            my.df <- Peak.list[which(Peak.list$metabolite_group %in% get.mg[i]), ]
            if (nrow(my.df) == 1) {
            } else {
                if (nrow(my.df) == 2) {
                  my.mat <- as.matrix(my.df[, grep(paste(QC.id, paste(strsplit(sexes, "(?<=.[_])", perl = TRUE),
                    collapse = "|"), sep = "|"), colnames(my.df))])
                  test.mat <- as.matrix(t(my.mat))
                  dimnames(test.mat) <- list(colnames(my.df[, sample.cols]), my.df$EIC_ID)
                  plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
                  text(x = 0.5, y = 0.5, paste("Dendrograms can't be plotted for two features.\n", "Use the correlation plot to the left to \n",
                    "Tell if two features are from the same metabolite"), cex = 1.6, col = "black")
                  res2 <- rcorr(test.mat)
                  M <- res2$r
                  P <- res2$P
                  # Insignificant correlations are leaved blank
                  order.hc2 <- corrplot::corrMatOrder(M, order = "hclust", hclust.method = "complete")
                  col.new = rainbow(maxlabel)[order.hc2]
                  corrplot(M, type = "upper", order = "hclust", hclust.method = "complete", p.mat = P, sig.level = 0.01, insig = "blank",
                    tl.col = col.new)

                  # Plots the EICs and Pseudo-spectra for each metabolite group containing more than one feature in a for loop
                  EIC.plots <- plotEICs(new_anposGa, pspec = get.mg[i], maxlabel = maxlabel, sleep = 0)  ## Plot the EICs

                  plotPsSpectrum(new_anposGa, pspec = get.mg[i], maxlabel = maxlabel, sleep = 0)  #Plot the Mass Spectra

                } else {
                  my.mat <- as.matrix(my.df[, grep(paste(QC.id, paste(strsplit(sexes, "(?<=.[_])", perl = TRUE),
                    collapse = "|"), sep = "|"), colnames(my.df))])
                  test.mat <- as.matrix(t(my.mat))
                  dimnames(test.mat) <- list(colnames(my.df[, sample.cols]), my.df$EIC_ID)
                  colnames(test.mat)
                  res <- cor(test.mat)
                  d <- dist(res)
                  sim.by.hclust <- flashClust(d)
                  attributes(d)
                  labels(d)
                  plot(sim.by.hclust)
                  res2 <- rcorr(test.mat)
                  M <- res2$r
                  P <- res2$P
                  # Insignificant correlation are crossed corrplot(M, type='upper', order='hclust', p.mat = P,
                  # sig.level = 0.01, insig = 'pch', pch = 3) Insignificant correlations are leaved blank
                  order.hc2 <- corrplot::corrMatOrder(M, order = "hclust", hclust.method = "complete")
                  col.new = rainbow(maxlabel)[order.hc2]
                  corrplot(M, type = "upper", order = "hclust", hclust.method = "complete", p.mat = P, sig.level = 0.01, insig = "blank",
                    tl.col = col.new)
                  # Plots the EICs and Pseudo-spectra for each metabolite group containing more than one feature in a for loop
                  EIC.plots <- plotEICs(new_anposGa, pspec = get.mg[i], maxlabel = maxlabel, sleep = 0)  ## Plot the EICs

                  plotPsSpectrum(new_anposGa, pspec = get.mg[i], maxlabel = maxlabel, sleep = 0)  #Plot the Mass Spectra
                }
            }
        }
        dev.off()
    }

    return(Peak.list.new)  #Use the second element in the list to validate CAMERA

}

#' @title Loop-based EIC Plotter for Ion Mode Duplicates
#'
#' @export
#' @description Plots EICs and Pseudospectra from (parsed) metabolite groups using plotEICs-methods and plotPsSpectrum-methods from CAMERA. Also plots correlation matrices and clustered dendrograms for each (parsed) metabolite group.
#' @param Peak.list a data frame containing combined ion mode peaklist with Duplicate IDs.  Alternatively can be retrieved from databases.  Default is NULL
#' @param gen.plots a logical indicating whether to create plots for metabolite groups.  Default is FALSE
#' @param anposGa xsannotate object with annotated isotopes and ion adducts and fragments for positive mode.
#' @param xpos.cor xcmsSet object shoud have grouping, retention time correction and fillPeaks applied.  Default is to look for this in anposGa
#' @param annegGa xsannotate object with annotated isotopes and ion adducts and fragments for negative mode.
#' @param xneg.cor xcmsSet object shoud have grouping, retention time correction and fillPeaks applied.  Default is to look for this in annegGa
#' @param file.base character string used to name graphical output. Default is 'EIC_plots'
#' @param QC.id character identifier for pooled QC samples. Default is 'Pooled_QC'
#' @param mytable character name of table in database to return
#' @param maxEIC numeric How many EICs to plot for each metabolite
#' @param ... parameters to be passed to database functions
#' @return NULL testing
#' @importFrom xcms getEIC
#' @importFrom graphics abline title
plot_ionduplicate = function(anposGa, xpos.cor, annegGa, xneg.cor, Peak.list, gen.plots, file.base, QC.id, mytable, maxEIC, ...) {
    if (missing(Peak.list))
        Peak.list = NULL
    if (missing(gen.plots))
        gen.plots = FALSE
    if (missing(xpos.cor))
        xpos.cor = anposGa@xcmsSet
    if (missing(xneg.cor))
        xneg.cor = annegGa@xcmsSet
    if (missing(QC.id))
        QC.id = "Pooled_QC"
    if (missing(mytable))
        mytable = NULL
    if (missing(file.base))
        file.base = "EIC_plots"
    if (missing(maxEIC))
        maxEIC = 2

    if (is.null(Peak.list)) {
        if (is.null(Peak.list) && is.null(mytable)) {
            stop("Must specify database parameters `myname` and `peak.db` if not providing Peak.list!", call. = FALSE)
        } else {
          # Read Peak.list from database
          Peak.list <- read_tbl(mytable, peak.db)
        }
    }


    # Read the file names for study samples and pooled QCs
    neg.filenames <- row.names(xneg.cor@phenoData)
    neg.QC.files <- neg.filenames[grep(QC.id, neg.filenames)]

    pos.filenames <- row.names(xpos.cor@phenoData)
    pos.QC.files <- pos.filenames[grep(QC.id, pos.filenames)]

    # Get the unique list of Duplicate IDs
    x <- sapply(Peak.list$Duplicate_ID, function(x) sum(as.numeric(Peak.list$Duplicate_ID == x)))
    x
    drops <- Peak.list$Duplicate_ID[x == 1]
    drops  #Duplicate IDs which only appear once

    List.ID <- Peak.list$Duplicate_ID
    List.ID  #List of all duplicate IDs
    res <- !List.ID %in% drops
    Un.ID <- unique(Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE)])
    Un.ID  #List of unique duplicate IDs to get EICs for

    # List of duplicate IDs for both positive and negative modes
    Dup.ID.Pos <- Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE) & sapply(Peak.list$`Ion Mode`, function(x) x ==
        "Pos")]
    Dup.ID.Neg <- Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE) & sapply(Peak.list$`Ion Mode`, function(x) x ==
        "Neg")]

    ## Code to plot EICs for all Duplicate IDs in a loop.####
    if (gen.plots) {
        total = length(Un.ID)

        pdf(file = paste(file.base, ".pdf", sep = ""))
        for (i in Un.ID) {
            EIC.table <- Peak.list[Peak.list$Duplicate_ID %in% i, ]  #Data frame containing features with adduct info
            res <- lapply(colnames(EIC.table), function(ch) grep(QC.id, ch))
            QC.list <- EIC.table[sapply(res, function(x) length(x) > 0)]

            EIC.list <- split(EIC.table, as.factor(EIC.table$`Ion Mode`))
            EIC.pos <- convert_EIC(EIC.list$Pos$EIC_ID)
            EIC.neg <- convert_EIC(EIC.list$Neg$EIC_ID)

            Name.pos <- as.character(EIC.list$Pos$Name)
            Name.neg <- as.character(EIC.list$Neg$Name)

            Index.Pos <- which(Dup.ID.Pos %in% i)
            Index.Neg <- which(Dup.ID.Neg %in% i)

            if (length(EIC.pos) > maxEIC)
              EIC.pos <- EIC.pos[1:maxEIC]
            if (length(EIC.neg) > maxEIC)
              EIC.neg <- EIC.neg[1:maxEIC]

            par(mar = c(5 + 2, 4, 4, 2) + 0.1)

            if (length(EIC.list) == 1) {

            } else if (length(EIC.pos) * length(EIC.neg) != 0) {
                ## has duplicate in both modes
                pos = list()
                cat("Be patient! EIC group No.", which(Un.ID[] %in% i), "out of a total of", length(Un.ID), "\n")
                for (i in 1:length(EIC.pos)) {
                  pos[i] = getEIC(xpos.cor, rt = "corrected", groupidx = EIC.pos[i], sampleidx = pos.QC.files)
                }
                neg = list()
                for (i in 1:length(EIC.neg)) {
                  neg[i] = getEIC(xneg.cor, rt = "corrected", groupidx = EIC.neg[i], sampleidx = neg.QC.files)
                }
                rt <- neg[[1]]@rtrange
                rt.min <- min(rt[, "rtmin"])
                rt.max <- max(rt[, "rtmax"])
                color.palette <- c("black", "blue", "mediumslateblue", "red", "green", "khaki3", "lavenderblush3",
                  "lawngreen", "lightseagreen", "purple", "skyblue1", "darkgreen", "darkgoldenrod1", "seashell4",
                  "sienna3", "salmon4", "hotpink")
                Adduct.No <- length(EIC.pos) + length(EIC.neg)
                nrow <- max(length(EIC.pos),length(EIC.neg))
                if(length(EIC.pos) == length(EIC.neg)) {
                  mymatrix <- matrix(c(1:Adduct.No), nrow = nrow, ncol = 2, byrow = FALSE)
                } else {
                  p = length(EIC.pos)
                  n = length(EIC.neg)
                  if(p<n) {
                    new.p <- c(1:p,rep(0,times = n-p))
                    new.n <- c(max(p)+1:n)
                    mydim <- c(new.p,new.n)
                  } else {
                    new.p <- c(1:p)
                    new.n <- c(max(p)+1:n,rep(0,times = p-n))
                    mydim <- c(new.p,new.n)
                  }
                  mymatrix <- matrix(mydim, nrow = nrow, ncol = 2, byrow = FALSE)
                }
                layout(mymatrix)

                for (j in 1:length(EIC.pos)) {
                  rt <- EIC.table[which(EIC.table$`Ion Mode` %in% "Pos"), "meanRT"] * 60
                  QCmax <- max(QC.list[which(EIC.table$`Ion Mode` %in% "Pos"), ])
                  plot(pos[[j]], col = color.palette, rtrange = cbind(rt.min, rt.max))
                  title(sub = paste(paste("Positive #", Index.Pos[j], ", EIC_ID:", EIC.pos[j]), paste(strwrap(Name.pos[j],
                    width = 0.9 * getOption("width")), collapse = "\n"), sep = "\n"), cex.sub = 1, col.sub = "limegreen",
                    line = 6)
                  abline(v = rt, col = "blue", lty = 2)
                }
                for (j in 1:length(EIC.neg)) {
                  rt <- EIC.table[which(EIC.table$`Ion Mode` %in% "Neg"), "meanRT"] * 60
                  QCmax <- max(QC.list[which(EIC.table$`Ion Mode` %in% "Neg"), ])
                  plot(neg[[j]], col = color.palette, rtrange = cbind(rt.min, rt.max))
                  title(sub = paste(paste("Negative #", Index.Neg[j], ", EIC_ID:", EIC.neg[j]), paste(strwrap(Name.neg[j],
                    width = 0.9 * getOption("width")), collapse = "\n"), sep = "\n"), cex.sub = 1, col.sub = "red",
                    line = 6)
                  abline(v = rt, col = "blue", lty = 2)
                }
            } else
              if (length(EIC.pos) == 0) {
                ## has duplicate in negative mode ONLY
                neg = list()
                for (j in 1:length(EIC.neg)) {
                  neg[j] = getEIC(xneg.cor, rt = "corrected", groupidx = EIC.neg[j], sampleidx = neg.QC.files)
                }
                rt <- neg[[1]]@rtrange
                rt.min <- min(rt[, "rtmin"])
                rt.max <- max(rt[, "rtmax"])
                color.palette <- c("black", "blue", "mediumslateblue", "red", "green", "khaki3", "lavenderblush3",
                  "lawngreen", "lightseagreen", "purple", "skyblue1", "darkgreen", "darkgoldenrod1", "seashell4",
                  "sienna3", "salmon4", "hotpink")
                Adduct.No <- length(EIC.pos) + length(EIC.neg)
                layout(matrix(c(1:Adduct.No), nrow = Adduct.No, ncol = 1, byrow = TRUE))
                for (j in 1:length(EIC.neg)) {
                  rt <- EIC.table[which(EIC.table$`Ion Mode` %in% "Neg"), "meanRT"] * 60
                  QCmax <- max(QC.list[which(EIC.table$`Ion Mode` %in% "Neg"), ])
                  plot(neg[[j]], col = color.palette, rtrange = cbind(rt.min, rt.max))
                  title(sub = paste(paste("Negative #", Index.Neg[j], ", EIC_ID:", EIC.neg[j]), paste(strwrap(Name.neg[j],
                    width = 0.9 * getOption("width")), collapse = "\n"), sep = "\n"), cex.sub = 1, col.sub = "red",
                    line = 6)
                  abline(v = rt, col = "blue", lty = 2)
                }
            } else {
                ## has duplicate in positive mode ONLY
                pos = list()
                for (j in 1:length(EIC.pos)) {
                  pos[j] = getEIC(xpos.cor, rt = "corrected", groupidx = EIC.pos[j], sampleidx = pos.QC.files)
                }
                rt <- pos[[1]]@rtrange
                rt.min <- min(rt[, "rtmin"])
                rt.max <- max(rt[, "rtmax"])
                color.palette <- c("black", "blue", "mediumslateblue", "red", "green", "khaki3", "lavenderblush3",
                  "lawngreen", "lightseagreen", "purple", "skyblue1", "darkgreen", "darkgoldenrod1", "seashell4",
                  "sienna3", "salmon4", "hotpink")
                Adduct.No <- length(EIC.pos) + length(EIC.neg)
                layout(matrix(c(1:Adduct.No), nrow = Adduct.No, ncol = 1, byrow = TRUE))
                for (j in 1:length(EIC.pos)) {
                  rt <- EIC.table[which(EIC.table$`Ion Mode` %in% "Pos"), "meanRT"] * 60
                  QCmax <- max(QC.list[which(EIC.table$`Ion Mode` %in% "Pos"), ])
                  plot(pos[[j]], col = color.palette, rtrange = cbind(rt.min, rt.max))
                  title(sub = paste(paste("Positive #", Index.Pos[j], ", EIC_ID:", EIC.pos[j]), paste(strwrap(Name.pos[j],
                    width = 0.9 * getOption("width")), collapse = "\n"), sep = "\n"), cex.sub = 1, col.sub = "limegreen",
                    line = 6)
                  abline(v = rt, col = "blue", lty = 2)
                }
            }
        }
        dev.off()
        x.pos <- rep(1, length.out = length(Dup.ID.Pos))  #needs to be as long as the number of positive plots in the PDF file
        x.neg <- rep(1, length.out = length(Dup.ID.Neg))  #needs to be as long as the number of negative plots in the PDF file
        return(list(x.pos, x.neg))
        ## End Plotting code####
    }

}

LUMA_order = function(object){
return(object)
}

convert_EIC = function(EIC) {
  EIC_split <- strsplit(EIC, split = ";")
  EIC_num <- lapply(EIC_split, function(x) as.numeric(x))
  EICs <- unlist(EIC_num)
  return(EICs)
}

validate_metgroup <- function(Peak.list.pspec) {
  validate.df <- Peak.list.pspec[order(Peak.list.pspec$metabolite_group), c("MS.ID", "mz", "rt", "Name",
                                                                            "Formula","Annotated.adduct",
                                                                            "isotopes", "adduct",
                                                                            "mono_mass", "metabolite_group",
                                                                            "Correlation.stat")]
  mytbl <- validate.df %>%
    group_by(metabolite_group)

  mytbl <- mytbl %>%
    summarise(max = max(Correlation.stat[Correlation.stat!=max(Correlation.stat)]),
              min = min(Correlation.stat))
  f1 <- which(mytbl$min >= 0.7)
  f2 <- which(mytbl$max <= 0.2)
  clear.groups <- mytbl[(union(f1,f2)),"metabolite_group"] %>%
    data.frame()
  TP.groups <- mytbl[f1,"metabolite_group"] %>%
    data.frame()
  TN.groups <- mytbl[f2,"metabolite_group"] %>%
    data.frame()
  muddy.groups <- mytbl[-(union(f1,f2)),"metabolite_group"] %>%
    data.frame()
  ##Error Check
  if(nrow(clear.groups) != nrow(TP.groups) + nrow(TN.groups))
    stop("The number of true positives and true negatives does not equal the number of clear cut metabolite groups!", call. = FALSE)
  n.col <- ncol(validate.df)
  col.names <- colnames(validate.df)
  validate.df[,(n.col+1:(n.col+3))] <- NA
  colnames(validate.df) <- c(col.names,"Category Score1","Category Score2","Category Score3")
  validate.df$`Category Score2` <- validate.df$Correlation.stat
  validate.df[which(validate.df$metabolite_group %in% TP.groups$metabolite_group),
              c("Category Score1","Category Score3")] <- 1
  validate.df[which(validate.df$metabolite_group %in% TN.groups$metabolite_group),
              c("Category Score1","Category Score3")] <- 0
  x <- which(validate.df$metabolite_group %in% clear.groups$metabolite_group)
  y <- which(validate.df$metabolite_group %in% muddy.groups$metabolite_group)
  myfactor <- factor(NA, levels = c("clear","muddy"))
  myfactor[x] <- "clear"
  myfactor[y] <- "muddy"
  validate.list <- split(validate.df,myfactor)
  return(validate.list)
}
