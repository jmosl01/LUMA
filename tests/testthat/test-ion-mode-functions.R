context("test-ion-mode-functions")

test_that("combines ion modes correctly", {

  if(require(lcmsfishdata, quietly = TRUE)) { file <-
  system.file("extdata/Search_parameters.txt", package = "lcmsfishdata")
  search.par <- read.table(file, header = TRUE, sep = "\t")
  class(method) <-  method <- "monoMass"
  Peak.list <- list(Positive = lcmsfishdata::Peaklist_db$Peaklist_Pos_Solvent_Peaks_Removed,
                    Negative = lcmsfishdata::Peaklist_db$Peaklist_Neg_Solvent_Peaks_Removed)

  test <- combine_ion_modes(Peak.list = Peak.list, search.par = search.par,
  method = method)
  test_colnames <- colnames(test)[-which(colnames(test) %in% colnames(Peak.list[["Positive"]]))] #Adds two new columns
  expect_equal(test_colnames, c("Duplicate_ID", "Duplicate_EIC"))
  expect_equal(length(which(duplicated(test[["Duplicate_ID"]]))),56) #number of ion mode duplicates found

  }
})

test_that("removes ion duplicates", {

  if(require(lcmsfishdata, quietly = TRUE)) {

    file <- system.file("extdata/EIC_index_pos.txt", package = "lcmsfishdata")
    EIC_index_pos <- read.table(file, header = TRUE, sep = "\t")
    file2 <- system.file("extdata/EIC_index_neg.txt", package = "lcmsfishdata")
    EIC_index_neg <- read.table(file2, header = TRUE, sep = "\t")

    Peak.list <- lcmsfishdata::Peaklist_db$Peaklist_Combined_with_Duplicate_IDs

    Key.list <- list(Positive = EIC_index_pos, Negative = EIC_index_neg)

    test <- remove_ion_dup(Peak.list = Peak.list,
                           Key.list = Key.list)
    expect_equal(55,nrow(Peak.list) - nrow(test)) ## Number of ion mode duplicates removed

  }


})

test_that("simply combining peaklists doesn't remove any features", {

  if(require(lcmsfishdata, quietly = TRUE)) {

    Peak.list <- list(Positive = lcmsfishdata::Peaklist_Pos[["From_CAMERA"]],
                      Negative = lcmsfishdata::Peaklist_Neg[["From_CAMERA"]])
    test <- pre_combine_ion_modes(Peak.list = Peak.list)
    test2 <- nrow(test) #Combined features from simply combining Peaklists
    expect_equal(test2,sum(sapply(Peak.list, nrow))) #Doesn't remove any ion mode duplicates
  }


})
