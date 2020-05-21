context("test-utilities")

test_that("Preprocess files names are correct", {
  expect_match(.set_PreProcessFileNames("Positive",TRUE)[[1]],"Blanks_Pos")
  expect_match(.set_PreProcessFileNames("Positive",TRUE)[[2]],"Blanks_Pos")
  expect_match(.set_PreProcessFileNames("Negative",TRUE)[[1]],"Blanks_Neg")
  expect_match(.set_PreProcessFileNames("Negative",TRUE)[[2]],"Blanks_Neg")
  expect_match(.set_PreProcessFileNames("Positive",FALSE)[[1]],"objects_Pos")
  expect_match(.set_PreProcessFileNames("Positive",FALSE)[[2]],"objects_Pos")
  expect_match(.set_PreProcessFileNames("Negative",FALSE)[[1]],"objects_Neg")
  expect_match(.set_PreProcessFileNames("Negative",FALSE)[[2]],"objects_Neg")
})

test_that("Two elements are the same", {
  expect_true(.SameElements(c(1, 2, 3), c(1, 3, 2)))
  expect_false(.SameElements(c(1, 2, 3), c(1, 1, 3, 2)))
})

test_that("EICs are converted to numbers", {
  expect_equal(.convert_EIC("1"),1)
  expect_null(.convert_EIC(NULL))

})

test_that("is a whole number", {
  expect_true(.isWhole(10))
  expect_false(.isWhole(0.1))
})
