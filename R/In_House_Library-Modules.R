#' @title Annotates Peak.list with In House Library
#'
#' @export
#' @description Compares isotope and adduct annotations within Peaklist to
#'   user-defined annotation library. For more details, see
#'   \code{match_Annotation}. For working examples, see \code{InitWorkflow,
#'   AnnotatePeaklist, CombineFeatures, CombinePeaklists, CullBackground,
#'   CullMF, CullCV, CullVoidVolume, FormatForMetaboAnalystR, FormatForSIMCA,
#'   NormalizePeaklists, ParseCAMERA, SimplyPeaklists, FinalWorkflow}.
#' @param from.table from which table should LUMA pull the Peak.list
#' @param to.table to which should LUMA save the modified Peak.list
#' @param lib.db character name of In House Library database. Default is
#'   'Annotated Library'
#' @param ... arguments to pass to match_Annotation.
#' @return NULL
AnnotatePeaklist <- function(from.table,to.table,lib.db,...) {

  #Initialize global variables
  lib_db <- NULL

  #Set default values
  if(missing(lib.db))
    lib.db <- "Annotated Library"

  #Initialize the library database
  lib_db <<- lib_db <- connect_libdb(lib.db)

  Peak.list <- read_tbl(mytable = from.table,
                        peak.db = peak_db)

  #Match Peak.list annotations against In House Library
  Peak.list.Annotated <- match_Annotation(Peak.list = Peak.list,
                                Annotated.library = cbind.data.frame(Name,
                                                                     Formula,
                                                                     Molecular.Weight,
                                                                     RT..Min.),
                                Library.phenodata = Library.phenodata,
                                rules = rules,
                                search.par = data.frame(ppm = ppm.cutoff,
                                                        rt = rt.cutoff,
                                                        Voidrt = Voidrt,
                                                        Corr.stat.pos = Corr.stat.pos,
                                                        Corr.stat.neg = Corr.stat.neg,
                                                        CV = cv.cutoff,
                                                        Minfrac = mf.cutoff,
                                                        Endogenous = Endogenous.thresh,
                                                        Solvent = Solvent.ratio,
                                                        gen.plots = gen.plots,
                                                        keep.singletons = keep.singletons),
                                IonMode = IonMode,
                                lib_db = lib_db,
                                ...)

  write_tbl(mydf = Peak.list.Annotated,
            peak.db = peak_db,
            myname = to.table)
}
