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

test_that("formats for MetaboAnalystR correctly", {
  if(require(lcmsfishdata, quietly = TRUE)) {
    file <- system.file("extdata/Sample_Class.txt", package = "LUMA")
    Sample.df <- read.table(file, header = TRUE, sep = "\t")
    file2 <- system.file("extdata/Sample_Data.csv", package = "LUMA")
    Sample.data <- read.table(file2, header = TRUE, sep = ",")
    Peak.list <- Peaklist_db$Peaklist_Normalized
    class(mSetObj) <- "pktable"
    new_mSetObj <- format_MetabolomicData(mSetObj = mSetObj, Peak.list =
                                            Peak.list, Sample.df = Sample.df, Sample.data = Sample.data)
    expect_equal(new_mSetObj$dataSet$orig.cls[1:5],rep("F_100_Effluent",5))

    class(mSetObj) <- "mass_all"
    new_mSetObj <- format_MetabolomicData(mSetObj = mSetObj, Peak.list =
                                            Peak.list, Sample.df = Sample.df, Sample.data = Sample.data)
    expect_null(new_mSetObj$dataSet$orig.cls)
  }
})
