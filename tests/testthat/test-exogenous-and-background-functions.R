context("test-exogenous-and-background-functions")

test_that("background removal works as planned", {

if(require(lcmsfishdata, quietly = TRUE)) {
file <- system.file('extdata/Sample_Class.txt', package = "LUMA")
Sample.df <- read.table(file, header = TRUE, sep = "\t")
file2 <- system.file('extdata/Search_parameters.txt', package = "LUMA")
search.par <- read.table(file2, header = TRUE, sep = "\t")

#From combined features
Peak.list <- list(pos = Peaklist_Pos_db$Trimmed_by_MinFrac, neg = Peaklist_Neg_db$Trimmed_by_MinFrac, blanks_pos = Blanks_Pos_db$Combined_Isotopes_and_Adducts, blanks_neg = Blanks_Neg_db$Combined_Isotopes_and_Adducts)
test <- remove_background_peaks(Peak.list = Peak.list, Sample.df = Sample.df, search.par = search.par, method = "monoMass", mem = TRUE)
expect_equal(sum(sapply(test, length)),4) #4 Peaklists are returned
expect_equal(names(test), c("Positive","Negative")) #Two ionization modes are returned
expect_equal(unname(sapply(test[["Positive"]], nrow)),c(454,36)) #The correct number of metabolites are returned and in the correct order
 }
})
