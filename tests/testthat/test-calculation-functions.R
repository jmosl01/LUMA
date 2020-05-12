context("test-calculation-functions")

test_that("Sums features", {
  library(LUMA)
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
