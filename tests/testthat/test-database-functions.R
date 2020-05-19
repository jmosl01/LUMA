context("test-database-functions")

test_that("filebases are generated", {
  file <- system.file("extdata/Sample_Data.csv", package =  "LUMA")
  sample_data <- read.table(file, header = T, sep = ",")
  mzdatafiles <- sample_data$CT.ID
  test.spos <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode = "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
  test.sneg <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode = "Negative", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
  test.bpos <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = TRUE, IonMode = "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
  test.bneg <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = TRUE, IonMode = "Negative", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
  expect_equal(test.spos, "Peaklist_Pos")
  expect_equal(test.sneg, "Peaklist_Neg")
  expect_equal(test.bpos, "Blanks_Pos")
  expect_equal(test.bneg, "Blanks_Neg")
})

test_that("valid databases are created", {
  file <- system.file("extdata/Sample_Data.csv", package =  "LUMA")
  sample_data <- read.table(file, header = T, sep = ",")
  mzdatafiles <- sample_data$CT.ID
  file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode = "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
  peak_db <- connect_peakdb(file.base = file.base, mem = T)
  lib_db <- connect_libdb(mem = T)
  samples.pos <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode = "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
  samples.neg <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode = "Negative", ion.id = c("Pos","Neg")) #Returns "Peaklist_Neg"
  blanks.pos <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = TRUE, IonMode = "Positive", ion.id = c("Pos","Neg")) #Returns "Blanks_Pos"
  blanks.neg <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = TRUE, IonMode = "Negative", ion.id = c("Pos","Neg")) #Returns "Blanks_Neg"

  spos_db <- connect_peakdb(file.base = samples.pos, mem = T)
  sneg_db <- connect_peakdb(file.base = samples.neg, mem = T)
  bpos_db <- connect_peakdb(file.base = blanks.pos, mem = T)
  bneg_db <- connect_peakdb(file.base = blanks.neg, mem = T)
  new_db.list <- connect_lumadb(db.list = c("spos_db","sneg_db","bpos_db","bneg_db"), mem = T)

  expect_true(dbIsValid(peak_db))
  expect_true(dbIsValid(lib_db))
  expect_true(all(sapply(new_db.list, function(x) { #All valid databases are created
    dbIsValid(x) })))
})
