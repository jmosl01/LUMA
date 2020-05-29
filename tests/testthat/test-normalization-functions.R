context("test-normalization-functions")

test_that("intensities are unit normalized", {

  if(require(lcmsfishdata, quietly = TRUE)) {

    file <- system.file("extdata/Sample_Class.txt", package = "lcmsfishdata")
    Sample.df <- read.table(file, header = TRUE, sep = "\t")
    Peak.list <- lcmsfishdata::Peaklist_db$Peaklist_Combined_FINAL
    QC.ID <- "Pooled_QC"
    mult <- 1

    test <- replace_zeros_normalize(Peak.list = Peak.list, Sample.df = Sample.df, QC.id = QC.ID, mult = mult)
    expect_equal(11,sum(test[,grep(QC.ID, colnames(test))])) #There are 11 Pooled QC samples, each with unit normalized intensity values

    mult2 <- 1000
    test2 <- replace_zeros_normalize(Peak.list = Peak.list, Sample.df = Sample.df, QC.id = QC.ID, mult = mult2)
    expect_equal(sum(test2[,grep(QC.ID, colnames(test))]),11000) #Now they are each normalized to 1000

    test_pos <- test[grep("Pos",test[["Ion.Mode"]]),grep(QC.ID, colnames(test))]
    test_neg <- test[grep("Neg",test[["Ion.Mode"]]),grep(QC.ID, colnames(test))]

    expect_equal(sum(test_pos), sum(test_neg)) #Ion modes are given equal weighting

  }


})
