context("test-plotting-functions")

test_that("generates correlation matrices for direct plotting of metabolite group plots correctly", {

  if(require(lcmsfishdata, quietly = TRUE)) {

    file <- system.file("extdata","Sample_Class.txt", package = "lcmsfishdata")
    Sample.df <- read.table(file, header = TRUE, sep = "\t")
    file2 <- system.file("extdata","CAMERA_objects_Pos.Rdata", package = "lcmsfishdata")
    load(file2, envir = environment())
    Peak.list <- lcmsfishdata::Peaklist_Pos[["input_parsed"]]
    file3 <- system.file("extdata","Sample_Data.csv", package = "lcmsfishdata")
    sample_data <- read.table(file3, header = TRUE, sep = ",")
    mzdatafiles <- sample_data$CT.ID

    file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode
                              = "Positive", ion.id = c("Pos","Neg"))

    test <- plot_metgroup(CAMERA.obj = anposGa, Sample.df = Sample.df,
                          Peak.list = Peak.list, center = 2, BLANK = FALSE, gen.plots = FALSE,
                          IonMode = "Positive", file.base = file.base, QC.id = "Pooled_QC_", maxlabel
                          = 10)
    expect_equal(class(test),"list") ##is list
    expect_equal(2,length(test)) ##with 2 elements
    test2 <- test[[1]]
    expect_equal("Correlation.stat",colnames(test2)[(which(!colnames(test2) %in% colnames(Peak.list)))]) #Adds new column
    expect_equal(8.142202535,sum(test2[["Correlation.stat"]][1:10]))
    expect_error(plot_metgroup(CAMERA.obj = anposGa, Sample.df = Sample.df,
                               Peak.list = Peak.list, center = 2, BLANK = FALSE, gen.plots = TRUE, #Doesn't work with lcmsfishdata; Doesn't have access to raw datafiles
                               IonMode = "Positive", file.base = file.base, QC.id = "Pooled_QC_", maxlabel
                               = 10))
    mydir <- getwd()
    if(grep("testthat", mydir) == 1) {
      myfile <- paste0(mydir,"/Peaklist_Pos_Corrplots.pdf")
    } else {
      myfile <- paste0(mydir, "/tests/testthat/Peaklist_Pos_Corrplots.pdf")
    }
    expect_true(file.exists(myfile)) #Plotting file is created
    dev.off()
    file.remove(myfile) #Cleanup change to working directory

  }

})

test_that("prepares peaklist for direct plotting of ion duplicates correctly", {

  if(require(lcmsfishdata, quietly = TRUE)) {

    file <- system.file("extdata","CAMERA_objects_Pos.Rdata", package = "lcmsfishdata")
    load(file, envir = environment())
    file2 <- system.file("extdata","CAMERA_objects_Neg.Rdata", package = "lcmsfishdata")
    load(file2, envir = environment())

    Peak.list <- lcmsfishdata::Peaklist_db[["Peaklist_Combined_with_Duplicate_IDs"]]

    test <- plot_ionduplicate(anposGa = anposGa, annegGa = annegGa, Peak.list = Peak.list, gen.plots = FALSE)
    expect_equal(class(test),"list") ##is list
    expect_equal(2,length(test)) ##with 2 elements
    test2 <- test[[1]]
    test3 <- test[[2]]
    expect_equal(54, sum(test2))
    expect_equal(58, sum(test3))
    expect_error(plot_ionduplicate(anposGa = anposGa, annegGa = annegGa, Peak.list = Peak.list, gen.plots = TRUE))
    dev.off()

    expect_true(file.exists("EIC_plots.pdf")) #Plotting file is created
    file.remove("EIC_plots.pdf") #Cleanup change to working directory
  }


})
