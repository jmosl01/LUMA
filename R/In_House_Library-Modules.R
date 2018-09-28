#' @title Searches In House Library
#'
#' @export
#' @description Compares isotope and adduct annotations within Peak.list (inherited from CAMERA) to user-defined annotation library and returns annotation results
#' @param Peak.list a table of class 'tbl_dbi', 'tbl_sql', 'tbl_lazy', or 'tbl' with samples as columns.  Should contain all output columns from XCMS and CAMERA, both metadata and sample data. Retention times must be in min.
#' @param Annotated.library a data frame with annotated metabolite entries. Must contain columns called 'Name', 'Formula', Molecular.Weight' and 'RT..Min.'.  Can contain additional info as separate columns.
#' @param rules a data frame containing the rule list used by CAMERA to annotate ion adducts and fragments.  Must contain the columns 'name','nmol','charge','massdiff','oidscore','quasi','ips'.
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @param lib.db character name of In House Library database.
#' Default is 'Annotated Library'
#' @return data frame containing the original table with added columns 'Name','MS.ID','Formula','Annotated.adduct' and any additional info columns from Annotated.library
search_IHL <- function(Peak.list, Annotated.library, rules, search.par, ion.mode, lib.db) {

  #Set default values
  if(missing(lib.db))
    lib.db <- "Annotated Library"

  #Initialize the library database
  lib_db <- connect_peakdb(lib.db,db.dir)

  #Match Peak.list annotations against In House Library
  Peak.list <- match_Annotation(Peak.list,Annotated.library,rules,search.par,ion.mode,lib_db)

  return(Peak.list)
}
