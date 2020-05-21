context("test-formatting-functions")

test_that("formats for SIMCA correctly", {
  if(require(lcmsfishdata, quietly = TRUE)) {
    file <- system.file('extdata/Sample_Class.txt', package = "lcmsfishdata")
    Sample.df <- read.table(file, header = TRUE, sep = "\t")
    file2 <- system.file('extdata/Sample_Data.csv', package = "lcmsfishdata")
    Sample.data <- read.table(file2, header = TRUE, sep = ",")
    Peak.list <- Peaklist_db$Peaklist_Normalized
    test <- format_simca(Peak.list = Peak.list, Sample.df = Sample.df, Sample.data = Sample.data)
    expect_true(all(!is.na(colnames(test)))) #Should be no NA's in colnames
    expect_true(all(!is.na((test[[1]])))) #Should be no NA's in first row
  }

})
