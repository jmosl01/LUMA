context("test-trimming-functions")

# test_that("Removes void volume", {
# })

test_that("Two elements are the same", {
  expect_true(.SameElements(c(1, 2, 3), c(1, 3, 2)))
  expect_false(.SameElements(c(1, 2, 3), c(1, 1, 3, 2)))
})

test_that("EICs are converted to numbers", {
  expect_equal(.convert_EIC("1"),1)
  expect_null(.convert_EIC(NULL))

})
