#' @title Selects retention time column
#'
#' @description Allows the user to select the retention time column from their Annotated Library
#' @param annotated_library a data frame with annotated metabolite entries. Must contain a column for each of `name`, `formula`, `molecular weight`, and `retention time`.  Can contain additional info as separate columns.
#' @importFrom svDialogs dlg_list
#' @return `character` multiple selections from column names in `annotated_library`
select_rt = function(annotated_library){
  choices <- colnames(annotated_library)
  res <- dlg_list(colnames(annotated_library), multiple = TRUE,
                  title = "Select the following four columns, in order, from your Annotated Library:
                  \nname \nformula \nmolecular weight \nretention time\n")$res
  if (!length(res)) {
    cat("You cancelled the choice\n")
  } else {
    cat("You selected:\n")
    print(res)
  }
  return(res)
}
