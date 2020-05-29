context("test-parser-functions")

test_that("Peaklists are parsed correctly", {

  if(require(lcmsfishdata, quietly = TRUE)) {

    file <- system.file("extdata/primary_adducts_pos.csv", package = "lcmsfishdata")
    rules <- read.csv(file, header = TRUE)
    Peak.list <- as.data.frame(lcmsfishdata::Peaklist_Pos[["Annotated"]])
    test <- parse_pos_results(raw = Peak.list, rule = rules, IonMode = "Positive")
    test_colnames <- colnames(test)[-which(colnames(test) %in% colnames(Peak.list))] ##Adds the following columns
    expect_identical(test_colnames, c("mono_mass","metabolite_group","monoisotopic_flg","adduct_flg",
                                      "isotope_flg","ambiguity_flg","Selection_flg"))
    expect_equal(343,length(unique(Peak.list[["pcgroup"]]))) #Originally were this many metabolite groupings
    expect_equal(3347,length(unique(test[["metabolite_group"]]))) #Now there are many more metabolite groupings

    file <- system.file("extdata/primary_adducts_neg.csv", package = "lcmsfishdata")
    rules <- read.csv(file, header = TRUE)
    Peak.list <- as.data.frame(lcmsfishdata::Peaklist_Neg[["Annotated"]])
    test <- parse_neg_results(raw = Peak.list, rule = rules, IonMode = "Negative")
    test_colnames <- colnames(test)[-which(colnames(test) %in% colnames(Peak.list))] ##Adds the following columns
    expect_identical(test_colnames, c("mono_mass","metabolite_group","monoisotopic_flg","adduct_flg",
                                      "isotope_flg","ambiguity_flg","Selection_flg"))
    expect_equal(264,length(unique(Peak.list[["pcgroup"]]))) #Originally were this many metabolite groupings
    expect_equal(2196,length(unique(test[["metabolite_group"]]))) #Now there are many more metabolite groupings

  }


})
