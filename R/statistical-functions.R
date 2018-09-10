#' @title run_Fold_Change
#'
#' @export
#' @description  Calculates fold change relative to controls for an untargeted metabolomics dataset
#' @param values a dataframe of numeric values from a metabolomics dataset with samples (rows) and metabolite intensities (cols)
#' @param classes a vector of factor values describing the sample class
#' @param control string name of the relative class against which fold change will be calculated
#' @importFrom utils head
#' @return a dataframe with first column containing descriptors for mean fold change values and standard error of the means. Contains values for each sample class (rows)
#'  and metabolite (cols).
#' @examples
#' run_Fold_Change(metabolites,grouping,"Control")
run_Fold_Change <- function(values, classes, control) {
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
#' @export
#' @description  Get p-value and other statistical parameters from univariate statistical analysis of metabolomics dataset
#' @param values a dataframe of numeric values from a metabolomics dataset with samples(rows) and metabolite intensity values (cols)
#' @param classes a vector of factor values describing the sample class
#' @param stat.method which statistical method to use. Can be ANOVA, ANOVA_FDR, or Ttest
#' @importFrom stats t.test aov
#' @importFrom fdrtool fdrtool
#' @return a matrix containing p-values, F-stats, degrees of freedom for class and
#' degrees of freedom for residuals for each metabolite (row)
#' @examples
#' get_pvalue(metabolites, grouping, "ANOVA")
#' get_pvalue(metabolites, grouping, "Ttest")
#'
#' subset <- which(grouping == "100_Effluent"|grouping == "Control")
#' mydata <- metabolites[subset,]
#' myclass <- grouping[subset]
#' get_pvalue(mydata,myclass,"Ttest")
get_pvalue <- function(values, classes, stat.method) {
  class(values) <- append(class(values),stat.method)
  mystats <- run_stats(values, classes)

  return(mystats)
}

run_stats <- function(values, ...) {
  UseMethod("run_stats", values)
}

run_stats.ANOVA <- function(values, classes) {
  print(paste("Running ANOVA on ", ncol(values),
        "m/z features."))

  test <- lapply(values, function(x) aov(x~classes))
  test.summary <- lapply(test, function(x) summary(x))
  names <- names(test.summary)
  p.value<-sapply(names, function(x) test.summary[[x]][[1]][["Pr(>F)"]][1])
  DF_Class<-sapply(names, function(x) test.summary[[x]][[1]][["Df"]][1])
  DF_Resid<-sapply(names, function(x) test.summary[[x]][[1]][["Df"]][2])
  t.score<-sapply(names, function(x) test.summary[[x]][[1]][["F value"]][1])
  mystats <- cbind.data.frame(m.z = as.double(names),p.value,DF_Class,DF_Resid,t.score)
  return(mystats)
}

run_stats.Ttest <- function(values, classes) {
  if(!is.factor(classes)) {
    classes <- as.factor(classes)
  }
  if(length(unique(classes)) > 2) {
    print("Student's t-test cannot be run on factors with more than two levels!")
    return(NULL)
  } else {
    print(paste("Running T-test on ", ncol(values),
          "m/z features."))

    test <- lapply(values, function(x) t.test(as.numeric(x)~classes, values = values))
    names <- names(test)
    p.value <- sapply(names, function(x) test[[x]]$p.value[1])
    CI.low <- sapply(names, function(x) test[[x]]$conf.int[1])
    CI.high <- sapply(names, function(x) test[[x]]$conf.int[2])
    t.score <- sapply(names, function(x) test[[x]]$statistic[1])
    mystats <- cbind.data.frame(m.z = as.double(names),p.value,CI.low,CI.high,t.score)
    return(mystats)
  }

}

run_stats.ANOVA_FDR <- function(values, classes) {
  print(paste("Running ANOVA with FDR on ", ncol(values),
              "m/z features."))
  test <- lapply(values, function(x) aov(x~classes))
  test.summary <- lapply(test, function(x) summary(x))
  names <- names(test.summary)
  pvalue<-sapply(names, function(x) test.summary[[x]][[1]][["Pr(>F)"]][1])
  DF_Class<-sapply(names, function(x) test.summary[[x]][[1]][["Df"]][1])
  DF_Resid<-sapply(names, function(x) test.summary[[x]][[1]][["Df"]][2])
  t.score<-sapply(names, function(x) test.summary[[x]][[1]][["F value"]][1])
  adj.Strimmer <- fdrtool(pvalue, statistic = "pvalue")
  p.value<-adj.Strimmer$qval
  Strimmer.locFDR<-adj.Strimmer$lfdr
  mystats <- cbind.data.frame(m.z = as.double(names),pvalue,DF_Class,DF_Resid,t.score,p.value,Strimmer.locFDR)
  return(mystats)
}
