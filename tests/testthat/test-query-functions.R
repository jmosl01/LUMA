context("test-query-functions")

test_that("mz features are annotated correctly", {

  if(require(RSQLite, quietly = TRUE)) {
    Peak.list <- Peaklist_Pos$From_CAMERA
    file <- system.file('extdata/primary_adducts_pos.csv', package = "LUMA")
    rules <- read.table(file, header = TRUE, sep = ",")
    file2 <- system.file('extdata/Search_Parameters.txt', package = "LUMA")
    search.par <- read.table(file2, header = TRUE, sep = "\t")
    file3 <- system.file('extdata/Annotated_library.csv', package = "LUMA")
    temp.library <- read.table(file3, header = TRUE, sep = ",")
    Library.phenodata <- temp.library[,-which(colnames(temp.library) %in% c("Name","Formula","Molecular.Weight","RT..Min.","Ion.Mode"))]
    Annotated.library <- cbind.data.frame(Name = temp.library[["Name"]],
                                         Formula = temp.library[["Formula"]],
                                         Molecular.Weight = temp.library[["Molecular.Weight"]],
                                         RT..Min. = temp.library[["RT..Min."]],
                                         Ion.Mode = temp.library[["Ion.Mode"]])

    lib_db <- connect_libdb(lib.db = "Annotated_library", mem = TRUE)
    test <- match_Annotation(Peak.list = Peak.list, Annotated.library = Annotated.library,
                             Library.phenodata = Library.phenodata, rules = rules,
                             search.par = search.par, IonMode = "Positive", lib_db = lib_db)
    expect_equal(grep("GUANOSINE",test$Name),19)
    expect_equal(grep("INOSINE",test$Name),c(15,21,22))
    dbDisconnect(lib_db)
  }
})


test_that("ion duplicates are found", {

  if(require(lcmsfishdata, quietly = TRUE)) {
    file <- system.file("extdata/Search_parameters.txt", package = "lcmsfishdata")
    search.par <- read.table(file, header = TRUE, sep = "\t")
    class(method) <- method <- "monoMass"
    Peak.list.neg <- Peaklist_db$Peaklist_Neg_Solvent_Peaks_Removed
    Peak.list.pos <- Peaklist_db$Peaklist_Pos_Solvent_Peaks_Removed
    test <- search_IonDup(method, Peak.list.pos = Peak.list.pos,
                          Peak.list.neg = Peak.list.neg, search.par = search.par)
    test1 <- colnames(test)[-which(colnames(test) %in% colnames(Peak.list.pos))] #Adds two new columns
    test2 <- length(which(duplicated(test[["Duplicate_ID"]]))) #number of ion mode duplicates found
    expect_equal(length(test1), 2) #Two new columns are added
    expect_equal(test1, c("Duplicate_ID","Duplicate_EIC")) #With these column names
    expect_equal(test2, 56) #Correct number of ion mode duplicates in test data
  }

})
