#' @title Runs \emph{mummichog}
#'
#' @export
#' @description Runs \emph{mummichog} within the R environment.  Basically a wrapper script built from the \bold{"MS Peaks to Pathways"} vignette from MetaboAnalystR.
#' @param Sample.df a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns 'Sex','Class','n','Endogenous'.
#' @param Peak.list a table of class 'tbl_df',tbl' or 'data.frame' with variables as columns.  Should contain all output columns from XCMS and CAMERA.
#' @param ion.mode a character string defining the ionization mode.  Must be either 'Positive' or 'Negative'
#' @param stat.method statistical method to use. Can be "TTest" or "ANOVA"
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns 'ppm','rt','Voidrt','Corr.stat.pos','Corr.stat.neg','CV','Minfrac','Endogenous','Solvent','gen.plots','keep.singletons'.
#' @param instrumentOpt Define the mass-spec instrument used to perform untargeted metabolomics
#' @param msModeOpt Define the mass-spec mode of the instrument used to perform untargeted metabolomics
#' @param pvalCutoff Numeric, specify the p-value cutoff to define significant m/z features from reference m/z features
#' @param lib Input the name of the organism library, default is hsa
#' @param enrichOpt Input the method to perform enrichment analysis
#' @param pvalOpt Input the method to calculate p-values
#' @param permNum Numeric, the number of permutations to perform
#' @importFrom MetaboAnalystR InitDataObjects UpdateMummichogParameters SanityCheckMummichogData PerformMummichog
#' @return NULL testing
run_mummichog = function(Sample.df,Peak.list,ion.mode,stat.method,search.par,
                         instrumentOpt,msModeOpt,pvalCutoff,lib,enrichOpt,pvalOpt,permNum) {
  # Create objects for storing processed data from the MS peaks to pathways module
  mSet <- InitDataObjects("conc","stat", FALSE)

  # Generate peak-list data and store in object
  Peak.list <- format_stats(Peak.list = Peak.list,
                           Sample.df = Sample.df)
  Peak.df <- Peak.list[[1]]
  control <- Peak.list[[2]]
  Peak.stats <- get_pvalue(values = Peak.df[,-1],
                        classes = Peak.df[,1],
                        stat.method = stat.method)

  mSet2 <- set_PeakListData(Peak.stats,mSet,Peak.df,control)
  # Set parameters for analysis
  mSet3 <- UpdateMummichogParameters(mSet2,
                                     instrumentOpt = instrumentOpt,
                                     msModeOpt = msModeOpt,
                                     pvalCutoff = pvalCutoff)

  # Sanity check for uploaded data
  mSet3 <- SanityCheckMummichogData(mSet3)

  # Perform the mummichog algorithm
  mSet4 <- PerformMummichog(mSet3,
                            lib = lib,
                            enrichOpt = enrichOpt,
                            pvalOpt = pvalOpt,
                            permNum = permNum)

  return(mSet4)
}


set_PeakListData = function(input,mSet,Peak.df,control) {
  mSetObj <- mSet
  hits <- c("p.value", "m.z", "t.score") %in% colnames(input)
  if(!all(hits)) {
    print("Missing information, data must contain three columns: 'p.value', 'm.z', 't.score'")
    return(0)
  }
  mSetObj$dataSet$orig <- cbind(input$p.value,input$m.z,input$t.score)
  mSetObj$msgSet$read.msg <- paste("A total of", length(input$p.value),
                                   "m/z features were found in your uploaded data.")
  return(mSetObj)
}
