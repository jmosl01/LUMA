#' @title Annotates Peak.list with In House Library
#'
#' @export
#' @description Compares isotope and adduct annotations within Peak.list (inherited from CAMERA) to user-defined annotation library.
#' See search_IHL and match_Annotation for more details.
#' @param from.table from which table should LUMA pull the Peak.list
#' @param to.table to which should LUMA save the modified Peak.list
#' @param lib.db character name of In House Library database.
#' Default is 'Annotated Library'
#' @return NULL
AnnotatePeaklist <- function(from.table,to.table,lib.db) {

  #Initialize global variables
  lib_db <- NULL

  #Set default values
  if(missing(lib.db))
    lib.db <- "Annotated Library"

  #Initialize the library database
  lib_db <<- lib_db <- connect_peakdb(lib.db)

  Peak.list <- read_tbl(mytable = from.table,
                        peak.db = peak_db)

  #Match Peak.list annotations against In House Library
  Peak.list.Annotated <- match_Annotation(Peak.list = Peak.list,
                                Annotated.library = data.frame(Name = Name,
                                                               Formula = Formula,
                                                               Molecular.Weight = Molecular.Weight,
                                                               RT..Min. = RT..Min.),
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
                                ion.mode = ion.mode,
                                lib_db = lib_db)

  write_tbl(mydf = Peak.list.Annotated,
            peak.db = peak_db,
            myname = to.table)
}
