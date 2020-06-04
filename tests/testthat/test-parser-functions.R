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

test_that("phenotype data is combined properly", {

  file <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
  search.par <- read.table(file, sep = "\t", header = TRUE) #Ignore Warning message
  file2 <- system.file('extdata/Sample_Class.txt', package = "LUMA")
  Sample.df <- read.table(file2, sep = "\t", header = TRUE) #Ignore Warning message
  Peak.list <- LUMA::Peaklist_Pos$output_parsed
  Summed.list <- sum_features(Peak.list = Peak.list, Sample.df = Sample.df ,
                              search.par = search.par, BLANK = FALSE, IonMode = "Positive")
  test <- combine_phenodata(Sample.df = Sample.df, Peak.list = Peak.list,
                            Summed.list = Summed.list, search.par = search.par, BLANK = FALSE, IonMode =
                              "Positive")
  expect_equal(11,nrow(Peak.list) - nrow(test)) ##Combines multiple features into single entries per metabolite
  expect_equal(33727702,test[which(test$metabolite_group %in% 230),"Pooled_QC_Pos_1"][[1]])
  expect_equal(6024959878,test[which(test$metabolite_group %in% 898),"Pooled_QC_Pos_1"][[1]])


})
