#' @title run_Fold_Change
#'
#' @export
#' @description  Calculates fold change relative to controls for an untargeted metabolomics dataset
#' @param values a dataframe of numeric values from a metabolomics dataset with samples (rows) and metabolite intensities (cols)
#' @param classes a vector of factor values describing the sample class, one of which must be Control
#' @importFrom utils head
#' @return a dataframe with first column containing descriptors for mean fold change values and standard error of the means. Contains values for each sample class (rows)
#'  and metabolite (cols).
#' @examples
#' run_Fold_Change(metabolites,grouping)
run_Fold_Change <- function(values, classes) {
  bin <- names(values)
  samples <- row.names(values)
  means <- lapply(values, function(x) tapply(as.numeric(x), classes, mean))
  ave.control <- lapply(bin, function(x) means[[x]][["Control"]])
  values[1,]
  sub.change <- lapply(samples, function(x) mapply("-", as.numeric(values[x,]), ave.control))
  fold.change <- as.data.frame(t(sapply(sub.change, function(x) mapply("/", x, ave.control))))
  colnames(fold.change) <- names(values)
  means <- sapply(fold.change, function(x) tapply(as.numeric(x), classes, mean))
  sd <- lapply(fold.change, function(x) tapply(as.numeric(x), classes, sd))
  sq.n <- lapply(fold.change, function(x) sqrt(tapply(as.numeric(x), classes, length)))
  se <- mapply("/", sd, sq.n)
  means
  rownames(means) <- paste("Fold Change",levels(classes), sep = "_")
  se
  rownames(se) <- paste("Standard Error of the Means", levels(classes), sep = "_")
  mydata <- rbind(means,se)
  mynames <- data.frame(x=rownames(mydata))
  myfoldchange <- cbind.data.frame(mynames,mydata)
  return(myfoldchange)
}

#' @title get_pvalue
#'
#' @description  Get p-value and other statistical parameters from univariate statistical analysis of metabolomics dataset
#' @param values a dataframe of numeric values from a metabolomics dataset with samples(rows) and metabolite intensity values (cols)
#' @param classes a vector of factor values describing the sample class
#' @param stat.method which statistical method to use. Can be ANOVA or Ttest
#' @return a matrix containing p-values, F-stats, degrees of freedom for class and
#' degrees of freedom for residuals for each metabolite (row)
#' @importFrom stats t.test
get_pvalue <- function(values, classes, stat.method) {
  class(values) <- c(class(values),stat.method)
  mystats <- run_stats(values, classes)
  return(mystats)
}

run_stats <- function(values, ...) {
  UseMethod("run_stats", values)
}

run_stats.ANOVA <- function(values, classes) {
  test <- lapply(values, function(x) aov(x~classes))
  test.summary <- lapply(test, function(x) summary(x))
  name <- names(test.summary)
  p.value<-sapply(name, function(x) test.summary[[x]][[1]][["Pr(>F)"]][1])
  DF_Class<-sapply(name, function(x) test.summary[[x]][[1]][["Df"]][1])
  DF_Resid<-sapply(name, function(x) test.summary[[x]][[1]][["Df"]][2])
  F_Stat<-sapply(name, function(x) test.summary[[x]][[1]][["F value"]][1])
  mystats <- cbind(name,p.value,DF_Class,DF_Resid,F_Stat)
  return(mystats)
}

run_stats.Ttest <- function(values, classes) {
  if(!is.factor(classes)) {
    classes <- as.factor(classes)
  }
  if(length(unique(classes)) > 2) {
    print("Student's t-test cannot be run on factors with more than two levels!")
  } else {
    # values <- values[,2:length(colnames(values))]
    test <- lapply(values, function(x) t.test(as.numeric(x)~classes, values = values))
    name <- names(test)
    p.value <- sapply(name, function(x) test[[x]]$p.value[1])
    CI.low <- sapply(name, function(x) test[[x]]$conf.int[1])
    CI.high <- sapply(name, function(x) test[[x]]$conf.int[2])
    f.stat <- sapply(name, function(x) test[[x]]$statistic[1])
    mystats <- cbind(name,p.value,CI.low,CI.high,f.stat)
    return(mystats)
  }

}
