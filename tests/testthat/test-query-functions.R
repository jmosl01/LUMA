context("test-query-functions")

test_that("mz features are annotated correctly", {

  if(require(lcmsfishdata, quietly = TRUE)) {
    Peak.list <- Peaklist_Pos$From_CAMERA
    file <- system.file('extdata/primary_adducts_pos.csv', package = "lcmsfishdata")
    rules <- read.table(file, header = TRUE, sep = ",")
    file2 <- system.file('extdata/Search_Parameters.txt', package = "lcmsfishdata")
    search.par <- read.table(file2, header = TRUE, sep = "\t")
    file3 <- system.file('extdata/Annotated_library.csv', package = "lcmsfishdata")
    temp.library <- read.table(file3, header = TRUE, sep = ",")
    Library.phenodata = temp.library[,-which(colnames(temp.library) %in% c("Name","Formula","Molecular.Weight","RT..Min."))]
    Annotated.library = cbind.data.frame(Name = temp.library[["Name"]],
                                         Formula = temp.library[["Formula"]],
                                         Molecular.Weight = temp.library[["Molecular.Weight"]],
                                         RT..Min. = temp.library[["RT..Min."]],
                                         Library.phenodata)

    lib_db <- connect_libdb(lib.db = "Annotated_library", mem = TRUE)
    test <- match_Annotation(Peak.list = Peak.list, Annotated.library = Annotated.library, rules = rules, search.par = search.par, IonMode = "Positive", lib_db = lib_db)
    expect_equal(grep("GUANOSINE",test$Name),19)
    expect_equal(grep("INOSINE",test$Name),c(15,21,22))
  }
})
