#' @title Runs _mummichog_
#'
#' @description Runs _mummichog_ within the R environment.  Basically a wrapper script built from the "MS Peaks to Pathways" vignette from MetaboAnalystR.
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Peak.list a table of class 'tbl_df',tbl' or 'data.frame' with variables as columns.  Should contain all output columns from XCMS and CAMERA.
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @param stat.method statistical method to use. Can be "TTest" or "ANOVA"
#' @param statistic what statistic should be used in _mummichog_. Can be "t-score" or "Fold change"
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @importFrom MetaboAnalystR InitDataObjects
#' @return NULL testing
run_mummichog = function(Sample.df,Peak.list,ion.mode,stat.method,search.par) {
  # Create objects for storing processed data from the MS peaks to pathways module
  mSet <- InitDataObjects("conc","stat", FALSE)

  # Generate peak-list data and store in object

  # Set parameters for analysis

  # Sanity check for uploaded data

  # Perform the mummichog algorithm

}


