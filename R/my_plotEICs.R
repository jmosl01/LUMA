#
# This function is a modified version of CAMERA's plotEICs.  It implements a bug fix that
# addresses an error with the color labels in the plots. Bug fix was first proposed by
# Jonathan Mosley.
#
my_plotEICs <- function (object, pspec = 1:length(object@pspectra), maxlabel = 0,
                         sleep = 0, method = "bin")
{
    smpls <- unique(object@psSamples[pspec])
    xeic <- new("xcmsEIC")
    xeic@rtrange <- matrix(nrow = length(pspec), ncol = 2)
    xeic@mzrange <- matrix(nrow = length(pspec), ncol = 2)
    pcpos <- 1
    rtmargin <- 1
    for (a in seq(along = smpls)) {
      xraw <- xcmsRaw(object@xcmsSet@filepaths[smpls[a]],
                      profmethod = method)
      pspecS <- pspec[which(object@psSamples[pspec] == smpls[a])]
      peaks <- CAMERA:::getPeaks(object@xcmsSet, smpls[a])
      eic <- lapply(pspecS, function(pc) {
        pidx <- object@pspectra[[pc]]
        pks <- peaks[pidx, , drop = FALSE]
        gks <- object@groupInfo[pidx, , drop = FALSE]
        nap <- which(is.na(pks[, 1]))
        pks[nap, ] <- cbind(gks[nap, c(1:6), drop = FALSE],
                            matrix(nrow = length(nap), ncol = 5, 0))
        bbox <- c(rtmin = min(pks[, "rtmin"]) - rtmargin,
                  rtmax = max(pks[, "rtmax"]) + rtmargin,
                  mzmin = min(pks[,"mzmin"]),
                  mzmax = max(pks[, "mzmax"]))
        eic <- xcms:::getEIC(xraw, rtrange = pks[, c("rtmin","rtmax"), drop = FALSE],
                             mzrange = pks[, c("mzmin","mzmax"), drop = FALSE])
        xeic@rtrange[pcpos, ] <<- bbox[c("rtmin", "rtmax")]
        xeic@mzrange[pcpos, ] <<- bbox[c("mzmin", "mzmax")]
        cat("-->", pcpos, "\n")
        pcpos <<- pcpos + 1
        eic@eic[[1]]
      })
      xeic@eic <- c(xeic@eic, eic)
    }
    names(xeic@eic) <- paste("Pseudospectrum ", pspec, sep = "")
    if (maxlabel > 0) {
      col <- rainbow(maxlabel)
    }
    else {
      col <- c()
    }
    for (ps in seq(along = pspec)) {
      EIC <- xeic@eic[[ps]]
      pidx <- object@pspectra[[pspec[ps]]]
      peaks <- CAMERA:::getPeaks(object@xcmsSet, object@psSamples[pspec[ps]])[pidx,, drop = FALSE]
      grps <- object@groupInfo[pidx, ]
      nap <- which(is.na(peaks[, 1]))
      naps <- rep(FALSE, nrow(peaks))
      if (length(nap) > 0) {
        naps[nap] <- TRUE
        peaks[nap, ] <- cbind(grps[nap, c(1:6), drop = FALSE],
                              matrix(nrow = length(nap), ncol = 5, 0))
      }
      main <- paste("Pseudospectrum ", pspec[ps], sep = "")
      neics <- length(pidx)
      lmaxlabel <- min(maxlabel, neics)
      eicidx <- 1:neics
      maxint <- numeric(length(eicidx))
      for (j in eicidx) {
        maxint[j] <- max(EIC[[j]][, "intensity"])
      }
      o <- order(maxint, decreasing = TRUE)
      rt <- xeic@rtrange[ps, ]
      rt.min <- round(mean(peaks[, "rtmin"]), digits = 3)
      rt.med <- round(peaks[o[1], "rt"], digits = 3)
      rt.max <- round(mean(peaks[, "rtmax"]), digits = 3)
      plot(0, 0, type = "n", xlim = rt, ylim = c(0, max(maxint)),
           xaxs = "i", xlab = "Retention Time (seconds)",
           ylab = "Intensity", main = paste("Extracted Ion Chromatograms for ",
                                            main, "\nTime: From", rt.min, "to", rt.max,
                                            ", mean", rt.med))
      lcol <- rgb(0.6, 0.6, 0.6)
      lcol <- c(col, rep(lcol, max(nrow(peaks) - maxlabel,0)))
      cnt <- 1
      for (j in eicidx[o]) {
        pts <- xeic@eic[[ps]][[j]]
        points(pts, type = "l", col = lcol[cnt])
        peakrange <- peaks[, c("rtmin", "rtmax"), drop = FALSE]
        ptsidx <- pts[, "rt"] >= peakrange[j, 1] & pts[,"rt"] <= peakrange[j, 2]
        if (naps[j]) {
          points(pts[ptsidx, ], type = "l", col = col[cnt],lwd = 1.3, lty = 3)
        }
        else {
          points(pts[ptsidx, ], type = "l", col = col[cnt],lwd = 1.3)
        }
        cnt <- cnt + 1
      }
      pspectrum <- getpspectra(object, grp = pspec[ps])
      mz <- pspectrum[o, "mz"]
      if (lmaxlabel > 0 & "adduct" %in% colnames(pspectrum)) {
        adduct <- sub("^ ", "", pspectrum[o, "adduct"])
        mass <- sapply(strsplit(adduct, " "), function(x) {
          x[2]
        })
        adduct <- sapply(strsplit(adduct, " "), function(x) {
          x[1]
        })
        umass <- unique(na.omit(mass[1:maxlabel]))
        adduct[is.na(adduct)] <- ""
        test <- vector("list", length = length(mz))
        mz <- format(pspectrum[o[1:lmaxlabel], "mz"],
                     digits = 5)
        if (length(umass) > 0) {
          for (i in 1:length(umass)) {
            ini <- which(mass == umass[i])
            for (ii in 1:length(ini)) {
              firstpart <- strsplit(adduct[ini[ii]],"M")[[1]][1]
              secondpart <- strsplit(adduct[ini[ii]],"M")[[1]][2]
              masspart <- mz[ini[ii]]
              test[[ini[ii]]] <- substitute(paste(masspart," ", firstpart, M[i], secondpart),
                                            list(firstpart = firstpart,
                                                 i = i,
                                                 secondpart = secondpart,
                                                 masspart = masspart))
            }
          }
        }
        for (i in seq(along = lmaxlabel)) {
          if (is.null(test[[i]])) {
            test[[i]] <- mz[i]
          }
        }
        leg <- as.expression(test[1:lmaxlabel])
        legend("topright", legend = leg, col = lcol,lty = 1)
      }
      if (sleep > 0) {
        Sys.sleep(sleep)
      }
    }
}
