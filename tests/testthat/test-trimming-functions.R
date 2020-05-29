context("test-trimming-functions")


test_that("background removal works as planned", {

  if(require(lcmsfishdata, quietly = TRUE)) {
    file <- system.file('extdata/Sample_Class.txt', package = "LUMA")
    Sample.df <- read.table(file, header = TRUE, sep = "\t")
    file2 <- system.file('extdata/Search_parameters.txt', package = "LUMA")
    search.par <- read.table(file2, header = TRUE, sep = "\t")

    #From combined features
    Peak.list <- list(pos = lcmsfishdata::Peaklist_Pos$Trimmed_by_MinFrac,
                      neg = lcmsfishdata::Peaklist_Neg$Trimmed_by_MinFrac,
                      blanks_pos = lcmsfishdata::Blanks_Pos$Combined_Isotopes_and_Adducts,
                      blanks_neg = lcmsfishdata::Blanks_Neg$Combined_Isotopes_and_Adducts)
    test <- remove_background_peaks(Peak.list = Peak.list, Sample.df = Sample.df,
                                    search.par = search.par, method = "monoMass", mem = TRUE)
    expect_equal(sum(sapply(test, length)),4) #4 Peaklists are returned
    expect_equal(names(test), c("Positive","Negative")) #Two ionization modes are returned
    expect_equal(unname(sapply(test[["Positive"]], nrow)),c(454,36)) #The correct number of metabolites are returned and in the correct order
  }
})

test_that("removes void volume", {
file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
Peak.list <- LUMA::Peaklist_Pos$From_CAMERA
test <- remove_void_volume(Peak.list = Peak.list, search.par = search.par, method = "mz")
expect_equal(nrow(Peak.list) - nrow(test),10)
})

test_that("trims by cv values", {
  file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
  search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
  Peak.list <- LUMA::Peaklist_Pos$From_CAMERA
  test <- trim_cv(Peak.list = Peak.list, search.par = search.par)


  Peak2.list <- LUMA::Peaklist_Pos$Combined_Isotopes_and_Adducts
  test2 <- trim_cv(Peak.list = Peak2.list, search.par = search.par)
  expect_equal(nrow(Peak.list) -  nrow(test), 13)
  expect_equal(nrow(Peak2.list) -  nrow(test2),9)
})

test_that("trims by minimum fraction values", {
  file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
  search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
  method = "mz"
  class(method) = method
  Peak.list <- LUMA::Peaklist_Pos$From_CAMERA_with_MinFrac
  test <- trim_minfrac(Peak.list = Peak.list, search.par = search.par, object = method)
  expect_equal(nrow(Peak.list) -  nrow(test),6)

  method = "monoMass"
  class(method) = method
  Peak2.list <- LUMA::Peaklist_Pos$Combined_Isotopes_and_Adducts
  test <- trim_minfrac(Peak.list = Peak2.list, search.par = search.par, object = method)
  expect_equal(nrow(Peak2.list) - nrow(test),4)
})
