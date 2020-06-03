#' @title Loop-based EIC Plotter for metabolite groups
#'
#' @export
#' @description Plots EICs and Pseudospectra from parsed metabolite groups using
#'   \code{plotEICs-methods, plotPsSpectrum-methods} from \code{CAMERA}. Also plots
#'   correlation matrices and clustered dendrograms for each parsed metabolite
#'   group.
#' @param Sample.df a data frame with class info as columns.  Must contain a
#'   separate row entry for each unique sex/class combination. Must contain the
#'   columns \code{"Sex","Class","n","Endogenous"}.
#' @param Peak.list a data frame from CAMERA that has been parsed.  Must contain
#'   all output columns from \code{XCMS, CAMERA, ParseCAMERA}.Can contain
#'   additional columns from \code{match_Annotation, calc_minfrac}.
#' @param center numeric value indicating which sample to pick for plotting
#'   purposes
#' @param BLANK a logical indicating whether blanks are being evaluated. Default
#'   is \code{FALSE}.
#' @param gen.plots a logical indicating whether to create plots for metabolite
#'   groups. Default is \code{FALSE}.
#' @param IonMode a character string defining the ionization mode.  Must be one
#'   of \code{c("Positive","Negative")}.
#' @param CAMERA.obj \code{xsannotate} object in parent environment with
#'   isotopes, ion adducts and fragments for positive mode.
#' @param file.base character string used to name graphical output.  Will be
#'   appended with \code{"_CorrPlots.pdf"}.
#' @param QC.id character identifier for pooled QC samples. Default is
#'   \code{"Pooled_QC"}.
#' @param maxlabel numeric How many m/z labels to print
#' @return List of length 2.  1st element is a data frame with all columns as
#'   the original data frame with column \code{"Correlation.stat"}. 2nd element
#'   is a list of objects used to validate CAMERA results.
#' @importFrom CAMERA plotEICs plotPsSpectrum
#' @importFrom grDevices dev.off pdf rainbow graphics.off
#' @importFrom graphics layout par plot text
#' @importFrom Hmisc rcorr
#' @importFrom corrplot corrplot corrMatOrder
#' @importFrom dplyr group_by
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'
#'   file <- system.file("extdata/Sample_Class.txt", package = "lcmsfishdata")
#'   Sample.df <- read.table(file, header = TRUE, sep = "\t")
#'   file2 <- system.file("extdata/CAMERA_objects_Pos.Rdata", package = "lcmsfishdata")
#'   load(file2, envir = environment())
#'   Peak.list <- lcmsfishdata::Peaklist_Pos[["input_parsed"]]
#'   file3 <- system.file("extdata/Sample_Data.csv", package = "lcmsfishdata")
#'   sample_data <- read.table(file3, header = TRUE, sep = ",")
#'   mzdatafiles <- sample_data$CT.ID
#'
#'   file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode
#'   = "Positive", ion.id = c("Pos","Neg"))
#'
#'   test <- plot_metgroup(CAMERA.obj = anposGa, Sample.df = Sample.df,
#'   Peak.list = Peak.list, center = 2, BLANK = FALSE, gen.plots = FALSE,
#'   IonMode = "Positive", file.base = file.base, QC.id = "Pooled_QC", maxlabel
#'   = 10)
#'   class(test) ##is list
#'   length(test) ## with 2 elements
#'   test2 <- test[[1]]
#'   colnames(test2)[(which(!colnames(test2) %in% colnames(Peak.list)))] #Adds new column
#'   test2[["Correlation.stat"]][1:10]
#'
#'
#'   \dontrun{
#'   #Runs with pdf plotting. This requires access to raw datafiles and won't
#'   work with lcmsfishdata. Better to use your own data here.
#'   test <- plot_metgroup(CAMERA.obj = anposGa, Sample.df = Sample.df,
#'   Peak.list = Peak.list, center = 2, BLANK = FALSE, gen.plots = TRUE, IonMode
#'   = "Positive", file.base = file.base, QC.id = "Pooled_QC", maxlabel = 10)
#'
#'   }
#' }
plot_metgroup = function(CAMERA.obj, Sample.df, Peak.list, center, BLANK, gen.plots, IonMode, file.base, QC.id, maxlabel) {


    #Set Default Values
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


    #Bind CAMERA objects from parent.frame
    new_CAMERA.obj <- CAMERA.obj

    new_CAMERA.obj@pspectra <- X

    ## Need to set a numeric value in this slot, else plotEICs will give error
    new_CAMERA.obj@sample <- c(1:nrow(new_CAMERA.obj@xcmsSet@phenoData))

    # Gets the EICs and Psspectra for each metabolite group that has more than one feature
    pspec.length <- sapply(new_CAMERA.obj@pspectra, function(x) length(x))
    get.mg <- which(pspec.length > 1)
    names(get.mg) <- NULL

    Peak.list.pspec <- calc_corrstat(Sample.df, Peak.list, get.mg, BLANK, IonMode)
    validate.list <- .validate_metgroup(Peak.list.pspec)

    Peak.list.new <- list(Peak.list.pspec, validate.list)
    # ! Very important!! Needs a sample index for each new psspectrum; without it, it will only find files ! for
    # the first n spectra, where n is the original number of psspectra. ! Find out how CAMERA performs automatic
    # selection, then replicate here for the new psspectra ! For now, I am using the center sample for every
    # psspectra
    new_CAMERA.obj@psSamples <- rep(center, length(X))  #Tells CAMERA to select all psspectra from the center file
    # new_anposGa@psSamples <- rep(-1, length(X)) #Experimental new_anposGa@sample <- as.numeric(NA) #doesn't work
    # for some reason

    # New code to plot correlation matrix plots for each metabolite group containing more than one feature
    if (gen.plots && !BLANK) {


        graphics.off()
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
                  text(x = 0.5, y = 0.5, paste("Dendrograms can't be plotted for two features.\n", "Use the correlation plot to the left to tell \n",
                    "if two features are from the same metabolite"), cex = 1.6, col = "black")
                  res2 <- rcorr(test.mat)
                  M <- res2$r
                  P <- res2$P
                  # Insignificant correlations are leaved blank
                  order.hc2 <- corrplot::corrMatOrder(M, order = "hclust", hclust.method = "complete")
                  col.new = rainbow(maxlabel)[order.hc2]
                  corrplot(M, type = "upper", order = "hclust", hclust.method = "complete", p.mat = P, sig.level = 0.01, insig = "blank",
                    tl.col = col.new, cl.ratio = 0.15, cl.align.text = "l", cl.cex = 0.6, cl.offset = 0.2)

                  # Plots the EICs and Pseudo-spectra for each metabolite group containing more than one feature in a for loop
                  EIC.plots <- plotEICs(new_CAMERA.obj, pspec = get.mg[i], maxlabel = maxlabel, sleep = 0)  ## Plot the EICs

                  a1 <- new_CAMERA.obj@pspectra[[get.mg[i]]]
                  b1 <- new_CAMERA.obj@groupInfo[a1,"mz"]
                  mz1 <- min(b1) * 0.80
                  mz2 <- max(b1) * 1.20

                  plotPsSpectrum(new_CAMERA.obj, pspec = get.mg[i], maxlabel = maxlabel, sleep = 0, cex.main = 0.8, cexMulti = 1.0, mzrange = c(mz1,mz2))  #Plot the Mass Spectra

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
                    tl.col = col.new, cl.ratio = 0.15, cl.align.text = "l", cl.cex = 0.6, cl.offset = 0.2)
                  # Plots the EICs and Pseudo-spectra for each metabolite group containing more than one feature in a for loop
                  EIC.plots <- plotEICs(new_CAMERA.obj, pspec = get.mg[i], maxlabel = maxlabel, sleep = 0)  ## Plot the EICs

                  a1 <- new_CAMERA.obj@pspectra[[get.mg[i]]]
                  b1 <- new_CAMERA.obj@groupInfo[a1,"mz"]
                  mz1 <- min(b1) * 0.80
                  mz2 <- max(b1) * 1.20

                  plotPsSpectrum(new_CAMERA.obj, pspec = get.mg[i], maxlabel = maxlabel, sleep = 0, cex.main = 0.8, cexMulti = 1.0, mzrange = c(mz1, mz2))  #Plot the Mass Spectra
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
#' @description Plots EICs from ion mode duplicates using
#'   \code{plotEICs-methods} from \code{CAMERA}.
#' @param Peak.list data.frame containing combined ion mode peaklist with column
#'   \code{"Duplicate_ID"}.  Alternatively can be retrieved from databases.
#'   Default is \code{NULL}.
#' @param gen.plots logical indicating whether to create plots for ion mode
#'   duplicates.  Default is \code{FALSE}.
#' @param anposGa \code{xsannotate} object in parent environment with
#'   isotopes, ion adducts and fragments for positive mode.
#' @param xpos \code{xcmsSet} object shoud have grouping, retention time
#'   correction and fillPeaks applied.  Default is to look for this in
#'   \code{anposGa}.
#' @param annegGa \code{xsannotate} object in parent environment with isotopes,
#'   ion adducts and fragments for negative mode.
#' @param xneg \code{xcmsSet} object shoud have grouping, retention time
#'   correction and fillPeaks applied.  Default is to look for this in
#'   \code{annegGa}.
#' @param rt.method Which method to use for EIC. Can be one of
#'   \code{c("corrected","raw")}.
#' @param file.base character string used to name graphical output. Default is
#'   \code{"EIC_plots"}.
#' @param QC.id character identifier for pooled QC samples. Default is
#'   \code{"Pooled_QC"}.
#' @param mytable character name of table in database to return
#' @param maxEIC numeric Max number of features for which to \code{plotEICs} for
#'   each metabolite.
#' @param maxQC numeric Max number of QCs used to \code{plotEICs}.
#' @param ... parameters to be passed to database functions
#' @return list of length 2 EIC indices for the ion duplicate plots
#' @importFrom xcms getEIC
#' @importFrom graphics abline title
#' @examples
#' library(LUMA)
#' if(require(lcmsfishdata, quietly = TRUE)) {
#'
#' file <- system.file("extdata/CAMERA_objects_Pos.Rdata", package = "lcmsfishdata")
#' load(file, envir = environment())
#' file2 <- system.file("extdata/CAMERA_objects_Neg.Rdata", package = "lcmsfishdata")
#' load(file2, envir = environment())
#'
#' Peak.list <- lcmsfishdata::Peaklist_db[["Peaklist_Combined_with_Duplicate_IDs"]]
#'
#' test <- plot_ionduplicate(anposGa = anposGa, annegGa = annegGa, Peak.list =
#' Peak.list, gen.plots = FALSE)
#' class(test) ##is list
#' length(test) ## with 2 elements
#'
#'   \dontrun{
#'   #Runs with pdf plotting. This requires access to raw datafiles and won't
#'   work with lcmsfishdata. Better to use your own data here.
#'   test <- plot_ionduplicate(anposGa = anposGa, annegGa = annegGa, Peak.list =
#'   Peak.list, gen.plots = TRUE)
#'
#'   }
#' }
plot_ionduplicate = function(anposGa, xpos, annegGa, xneg, rt.method, Peak.list, gen.plots,
                             file.base, QC.id, mytable, maxEIC, maxQC, ...) {

    #Bind CAMERA objects from parent.frame
    anposGa <- anposGa
    annegGa <- annegGa

    #Set default values
    if (missing(Peak.list))
        Peak.list = NULL
    if (missing(gen.plots))
        gen.plots = FALSE
    if (missing(xpos))
        xpos = anposGa@xcmsSet
    if (missing(xneg))
        xneg = annegGa@xcmsSet
    if (missing(rt.method))
      rt.method = "corrected"
    if (missing(QC.id))
        QC.id = "Pooled_QC"
    if (missing(mytable))
        mytable = NULL
    if (missing(file.base))
        file.base = "EIC_plots"
    if (missing(maxEIC))
        maxEIC = 2
    if (missing(maxQC))
        maxQC = 1

    if (is.null(Peak.list)) {
        if (is.null(Peak.list) && is.null(mytable)) {
            stop("Must specify database parameters `myname` and `peak.db` if not providing Peak.list!", call. = FALSE)
        } else {
          # Read Peak.list from database
          Peak.list <- read_tbl(mytable, ...)
        }
    }


  ## Read the file names for study samples and pooled QCs

    #Selects QC file(s) to plot EICs for Negative Mode
    neg.filenames <- row.names(xneg@phenoData)
    neg.QC.files <- neg.filenames[grep(QC.id, neg.filenames)]

    if(maxQC == 1) {
      mymat <- matrix(nrow = nrow(xneg@phenoData), ncol = 2)
      mymat[,1] <- neg.filenames
      for(i in seq_along(neg.filenames)) {
        raw <- xneg@rt$raw[[i]]
        corrected <- xneg@rt$corrected[[i]]
        mymat[i,2] <- identical(raw,corrected, num.eq = TRUE)
      }
      neg.QC.files <- mymat[which(mymat[,2] %in% TRUE),1]
    } else {
      if(maxQC > length(neg.QC.files)) {
        stop("maxQC cannot be larger than the number of QC files!")
      } else {
        rand.QC.ind <- sample(seq_along(neg.QC.files), size = maxQC, replace = FALSE)
        neg.QC.files <- neg.QC.files[rand.QC.ind]
      }
    }


    #Selects QC file(s) to plot EICs for Positive Mode
    pos.filenames <- row.names(xpos@phenoData)
    pos.QC.files <- pos.filenames[grep(QC.id, pos.filenames)]

    if(maxQC == 1) {
      mymat <- matrix(nrow = nrow(xpos@phenoData), ncol = 2)
      mymat[,1] <- pos.filenames
      for(i in seq_along(pos.filenames)) {
        raw <- xpos@rt$raw[[i]]
        corrected <- xpos@rt$corrected[[i]]
        mymat[i,2] <- identical(raw,corrected, num.eq = TRUE)
      }
      pos.QC.files <- mymat[which(mymat[,2] %in% TRUE),1]
    } else {
      if(maxQC > length(pos.QC.files)) {
        stop("maxQC cannot be larger than the number of QC files!")
      } else {
        pos.QC.files <- pos.QC.files[rand.QC.ind]
      }
    }

    # Get the unique list of Duplicate IDs
    x <- sapply(Peak.list$Duplicate_ID, function(x) sum(as.numeric(Peak.list$Duplicate_ID == x)))

    drops <- Peak.list$Duplicate_ID[x == 1]


    List.ID <- Peak.list$Duplicate_ID

    res <- !List.ID %in% drops
    Un.ID <- unique(Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE)])


    # List of duplicate IDs for both positive and negative modes
    Dup.ID.Pos <- Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE) &
                                           sapply(Peak.list$Ion.Mode, function(x) x == "Pos")]
    Dup.ID.Neg <- Peak.list$Duplicate_ID[sapply(res, function(x) x == TRUE) &
                                           sapply(Peak.list$Ion.Mode, function(x) x == "Neg")]

  ## Code to plot EICs for all Duplicate IDs in a loop.####
    if (gen.plots) {

        graphics.off()

        total = length(Un.ID)

        pdf(file = paste(file.base, ".pdf", sep = ""))
        for (i in Un.ID) {
            EIC.table <- Peak.list[Peak.list$Duplicate_ID %in% i, ]  #Data frame containing features with adduct info
            res <- lapply(colnames(EIC.table), function(ch) grep(QC.id, ch))
            QC.list <- EIC.table[sapply(res, function(x) length(x) > 0)]

            EIC.list <- split(EIC.table, as.factor(EIC.table$Ion.Mode))
            EIC.pos <- .convert_EIC(EIC.list$Pos$EIC_ID)
            EIC.neg <- .convert_EIC(EIC.list$Neg$EIC_ID)

            Name.pos <- as.character(EIC.list$Pos$Name)
            Name.neg <- as.character(EIC.list$Neg$Name)

            Index.Pos <- which(Dup.ID.Pos %in% i)
            Index.Neg <- which(Dup.ID.Neg %in% i)

            if (length(EIC.pos) > maxEIC)
              EIC.pos <- EIC.pos[1:maxEIC]
            if (length(EIC.neg) > maxEIC)
              EIC.neg <- EIC.neg[1:maxEIC]

            par(mar = c(5 + 2, 4, 4, 2) + 0.1)

            if (length(EIC.list) == 0) {

            } else if (length(EIC.pos) * length(EIC.neg) != 0) {
                ## has duplicate in both modes
                pos = list()
                cat("Be patient! EIC group No.", which(Un.ID[] %in% i), "out of a total of", length(Un.ID), "\n")
                for (i in 1:length(EIC.pos)) {
                  pos[[i]] = getEIC(xpos, rt = rt.method, groupidx = EIC.pos[i], sampleidx = pos.QC.files)
                }
                neg = list()
                for (i in 1:length(EIC.neg)) {
                  neg[[i]] = getEIC(xneg, rt = rt.method, groupidx = EIC.neg[i], sampleidx = neg.QC.files)
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
                  rt <- EIC.table[which(EIC.table$Ion.Mode %in% "Pos"), "meanRT"] * 60
                  rt <- as.numeric(unlist(rt))
                  QCmax <- max(QC.list[which(EIC.table$Ion.Mode %in% "Pos"), ])
                  plot(pos[[j]], col = color.palette, rtrange = cbind(rt.min, rt.max))
                  title(sub = paste(paste("Positive #", Index.Pos[j], ", EIC_ID:", EIC.pos[j]), paste(strwrap(Name.pos[j],
                    width = 0.9 * getOption("width")), collapse = "\n"), sep = "\n"), cex.sub = 1, col.sub = "limegreen",
                    line = 6)
                  abline(v = rt, col = "blue", lty = 2)
                }
                for (j in 1:length(EIC.neg)) {
                  rt <- EIC.table[which(EIC.table$Ion.Mode %in% "Neg"), "meanRT"] * 60
                  rt <- as.numeric(unlist(rt))
                  QCmax <- max(QC.list[which(EIC.table$Ion.Mode %in% "Neg"), ])
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
                cat("Neg mode only! EIC group No.", which(Un.ID[] %in% i), "out of a total of", length(Un.ID), "\n")
                for (j in 1:length(EIC.neg)) {
                  neg[[j]] = getEIC(xneg, rt = rt.method, groupidx = EIC.neg[j], sampleidx = neg.QC.files)
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
                  rt <- EIC.table[which(EIC.table$Ion.Mode %in% "Neg"), "meanRT"] * 60
                  rt <- as.numeric(unlist(rt))
                  QCmax <- max(QC.list[which(EIC.table$Ion.Mode %in% "Neg"), ])
                  plot(neg[[j]], col = color.palette, rtrange = cbind(rt.min, rt.max))
                  title(sub = paste(paste("Negative #", Index.Neg[j], ", EIC_ID:", EIC.neg[j]), paste(strwrap(Name.neg[j],
                                    width = 0.9 * getOption("width")), collapse = "\n"), sep = "\n"), cex.sub = 1, col.sub = "red",
                                    line = 6)
                  abline(v = rt, col = "blue", lty = 2)
                }
              } else {
                ## has duplicate in positive mode ONLY
                pos = list()
                cat("\n\nPos mode only! EIC group No.", which(Un.ID[] %in% i), "out of a total of", length(Un.ID), "\n\n\n")
                for (j in 1:length(EIC.pos)) {
                  pos[[j]] = getEIC(xpos, rt = rt.method, groupidx = EIC.pos[j], sampleidx = pos.QC.files)
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
                  rt <- EIC.table[which(EIC.table$Ion.Mode %in% "Pos"), "meanRT"] * 60
                  rt <- as.numeric(unlist(rt))
                  QCmax <- max(QC.list[which(EIC.table$Ion.Mode %in% "Pos"), ])
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
    } else {
      x.pos <- rep(1, length.out = length(Dup.ID.Pos))  #needs to be as long as the number of positive plots in the PDF file
      x.neg <- rep(1, length.out = length(Dup.ID.Neg))  #needs to be as long as the number of negative plots in the PDF file
      return(list(x.pos, x.neg))
    }

}

