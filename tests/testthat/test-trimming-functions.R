context("test-trimming-functions")

test_that("Removes void volume", {
library(LUMA)
file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
test <- remove_void_volume(Peak.list = Peaklist_Pos$From_CAMERA, search.par = search.par, method = "mz")
expect_equal(nrow(Peaklist_Pos$From_CAMERA) - nrow(test),10)
})

