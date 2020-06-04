context("test-calculation-functions")

test_that("Sums features", {
  file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
  search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
  file2 <- system.file('extdata/Sample_Class.txt', package = "LUMA")
  Sample.df <- read.table(file2, sep = "\t", header = TRUE) #Ignore Warning message
  Peak.list <- Peaklist_Pos$output_parsed
  Sum.Peak.list <- sum_features(Peak.list = Peak.list, Sample.df = Sample.df ,
                         search.par = search.par, BLANK = FALSE, IonMode = "Positive")
  my_int <- Peak.list[["metabolite_group"]] %in% Peak.list$metabolite_group[[1]]
  test <- sum(Peak.list[my_int,"Pooled_QC_Pos_1"])
  my_int2 <- Sum.Peak.list[["metabolite_group"]] %in% Sum.Peak.list$metabolite_group[[1]]
  test2 <- Sum.Peak.list[my_int2,"Pooled_QC_Pos_1"][[1]]
  expect_equal(test,test2)
})

test_that("calculates minimum fraction values", {
  if(require(lcmsfishdata, quietly = TRUE)) {
  file <- system.file('extdata/XCMS_objects_Pos.Rdata', package = "lcmsfishdata")
  load(file)
  file2 <- system.file('extdata/Sample_Class.txt', package = "LUMA")
  Sample.df <- read.table(file2, sep = "\t", header = TRUE) #Ignore Warning message
  test <- calc_minfrac(Sample.df = Sample.df, xset4 = xset4, BLANK = FALSE, Peak.list = lcmsfishdata::Peaklist_Pos$From_CAMERA)
  expect_equal(sum(test$MinFrac),4952.333333)
  }
})

test_that("calculates correlation statistic value", {
  if(require(lcmsfishdata, quietly = TRUE)) {

  file <- system.file('extdata/CAMERA_objects_Pos.Rdata', package = "lcmsfishdata")
  load(file)
  pspec.length <- sapply(anposGa@pspectra, function(x) length(x))
  get.mg <- which(pspec.length > 1)
  file2 <- system.file('extdata/Sample_Class.txt', package = "LUMA")
  Sample.df <- read.table(file2, sep = "\t", header = TRUE) #Ignore Warning message
  test <- calc_corrstat(Sample.df = Sample.df, Peak.list = lcmsfishdata::Peaklist_Pos$input_parsed, get.mg = get.mg, BLANK = FALSE, IonMode = "Positive")
  expect_equal(sum(test$Correlation.stat),261.1041613)
  }
})

test_that("calculates coefficient of variation", {
  test <- calc_cv(Peak.list = LUMA::Peaklist_Pos$From_CAMERA)
  expect_equal(sum(test$`X.CV`),9.848621392)
})
