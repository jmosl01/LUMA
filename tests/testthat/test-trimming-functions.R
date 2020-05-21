context("test-trimming-functions")

test_that("removes void volume", {
file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
test <- remove_void_volume(Peak.list = Peaklist_Pos$From_CAMERA, search.par = search.par, method = "mz")
expect_equal(nrow(Peaklist_Pos$From_CAMERA) - nrow(test),10)
})

test_that("trims by cv values", {
  file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
  search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
  test <- trim_cv(Peak.list = Peaklist_Pos$From_CAMERA, search.par = search.par)
  expect_equal(nrow(Peaklist_Pos$From_CAMERA) -  nrow(test), 14)
  test <- trim_cv(Peak.list = Peaklist_Pos$Combined_Isotopes_and_Adducts, search.par = search.par)
  expect_equal(nrow(Peaklist_Pos$Combined_Isotopes_and_Adducts) -  nrow(test),9)
})

test_that("trims by minimum fraction values", {
  file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
  search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
  method = "mz"
  class(method) = method
  test <- trim_minfrac(Peak.list = Peaklist_Pos$From_CAMERA_with_MinFrac, search.par = search.par, object = method)
  expect_equal(nrow(Peaklist_Pos$From_CAMERA_with_MinFrac) -  nrow(test),7)

  method = "monoMass"
  class(method) = method
  test <- trim_minfrac(Peak.list = Peaklist_Pos$Combined_Isotopes_and_Adducts, search.par = search.par, object = method)
  expect_equal(nrow(Peaklist_Pos$Combined_Isotopes_and_Adducts) - nrow(test),4)
})
