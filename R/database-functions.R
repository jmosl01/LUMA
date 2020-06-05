#' @title Generate Filebase
#'
#' @export
#' @description Generates filebase for reading and writing to databases
#' @param mzdatafiles character vector containing data files to process and
#'   store results in databases
#' @param BLANK a logical indicating whether blanks are being evaluated
#' @param IonMode a character string defining the ionization mode.  Must be
#'   either 'Positive' or 'Negative'
#' @param ion.id character vector of length 2 specifying identifier in filename
#'   designating positive or negative ionization mode.  Positive identifier must
#'   come first.
#' @return character
#' @examples
#' library(LUMA)
#' file <- system.file("extdata/Sample_Data.csv", package =  "LUMA")
#' sample_data <- read.table(file, header = TRUE, sep = ",")
#' mzdatafiles <- sample_data$CT.ID
#' test <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode =
#' "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
#' print(test)
gen_filebase = function(mzdatafiles, BLANK, ion.id, IonMode) {
  if (IonMode == "Positive" && BLANK == TRUE) {
    mzdatafiles <- subset(mzdatafiles, subset = grepl(paste(ion.id[1]), mzdatafiles, ignore.case = TRUE))
    file.base = "Blanks_Pos"
    mzdatafiles <- mzdatafiles[c(grep("Blanks", mzdatafiles, ignore.case = TRUE))]
  } else {
    if (IonMode == "Negative" && BLANK == TRUE) {
      mzdatafiles <- subset(mzdatafiles, subset = grepl(paste(ion.id[1]), mzdatafiles, ignore.case = TRUE))
      file.base = "Blanks_Neg"
      mzdatafiles <- mzdatafiles[c(grep("Blanks", mzdatafiles))]
    } else {
      if (IonMode == "Positive" && BLANK == FALSE) {
        mzdatafiles <- subset(mzdatafiles, subset = grepl(paste(ion.id[1]), mzdatafiles, ignore.case = TRUE))
        file.base = "Peaklist_Pos"
        mzdatafiles <- mzdatafiles[-c(grep("Blanks", mzdatafiles))]
      } else {
        if (IonMode == "Negative" && BLANK == FALSE) {
          mzdatafiles <- subset(mzdatafiles, subset = grepl(paste(ion.id[1]), mzdatafiles, ignore.case = TRUE))
          file.base = "Peaklist_Neg"
          mzdatafiles <- mzdatafiles[-c(grep("Blanks", mzdatafiles))]
        } else {
          cat("Ion mode must be Positive or Negative.\nBe sure to specify whether to analyze blanks with logical indicator")
          stop()
        }
      }
    }
  }
  return(file.base)
}

#' @title Connects to peak database
#'
#' @export
#' @description Establishes a connection to an SQLite database for storing
#'   Peak.list; if doesn't exist, creates new database
#' @param file.base character return from gen_filebase function
#' @param db.dir character what should the database directory be called.
#'   Default is 'db'
#' @param mem logical should database be in-memory. Default is FALSE
#' @return Formal class SQLiteConnection
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
#' @examples
#' library(LUMA)
#' if(require(RSQLite, quietly = TRUE)) {
#' file <- system.file("extdata/Sample_Data.csv", package =  "LUMA")
#' sample_data <- read.table(file, header = TRUE, sep = ",")
#' mzdatafiles <- sample_data$CT.ID
#' file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode =
#' "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
#' peak_db <- connect_peakdb(file.base = file.base, mem = TRUE)
#' dbIsValid(peak_db) #valid database is created
#' dbDisconnect(peak_db)
#' }
connect_peakdb = function(file.base, db.dir, mem) {

    #Set default variables
    if (missing(db.dir))
        db.dir = "db"
    if (missing(mem))
        mem = F

    ## Uncomment the following line if you dare
    # peak_db_file <- paste(file.base, db.dir, sep = "_")
    if(mem) {

      peak_db <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

    } else {

      dir.create(db.dir, recursive = FALSE, showWarnings = FALSE)

      peak_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(paste(db.dir, file.base, sep = "/"),".SQLite"))

    }

    return(peak_db)
}

#' @title Connects to library database
#'
#' @export
#' @description Establishes a connection to an SQLite database for searching
#'   Peak.list against library; if doesn't exist, creates new database
#' @param lib.db character name of database
#' @param db.dir character directory containing the database. Default is 'db'
#' @param mem logical should database be in-memory. Default is TRUE
#' @return Formal class SQLiteConnection
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
#' @examples
#' library(LUMA)
#' if(require(RSQLite, quietly = TRUE)) {
#' lib_db <- connect_libdb(mem = TRUE)
#' dbIsValid(lib_db) #Valid database is created
#' dbDisconnect(lib_db)
#' }
connect_libdb = function(lib.db, db.dir, mem) {

  #Set default variables
  if (missing(db.dir))
    db.dir = "db"
  if (missing(mem))
    mem = T


  if(mem) {

    lib_db <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")

  } else {

    dir.create(db.dir, recursive = FALSE, showWarnings = FALSE)

    lib_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(paste(db.dir, lib.db, sep = "/"),".SQLite"))

  }

  return(lib_db)
}

#' @title Connects to LUMA database
#'
#' @export
#' @description Establishes a connection to an RSQLite database for combining
#'   Peak.lists together from two different ionization modes. Must have
#'   previously saved SQLite databases to hard disk.
#' @param db.list list character names of databases containing results from
#'   processing positive mode and negative mode data for samples and blanks
#' @param db.dir character directory containing the databases
#' @param new.db character what should the new database be called.
#' @param mem logical should database be in-memory. Default is FALSE
#' @return list of Formal class SQLiteConnections, starting with new.db entry
#'   followed by one for each db.list entry
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
#' @examples
#' library(LUMA)
#' if(require(RSQLite, quietly = TRUE)) {
#' file <- system.file("extdata/Sample_Data.csv", package =  "LUMA")
#' sample_data <- read.table(file, header = TRUE, sep = ",")
#' mzdatafiles <- sample_data$CT.ID
#' samples.pos <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode
#' = "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos" samples.neg <-
#' gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode = "Negative",
#' ion.id = c("Pos","Neg")) #Returns "Peaklist_Neg" blanks.pos <-
#' gen_filebase(mzdatafiles = mzdatafiles, BLANK = TRUE, IonMode = "Positive",
#' ion.id = c("Pos","Neg")) #Returns "Blanks_Pos" blanks.neg <-
#' gen_filebase(mzdatafiles = mzdatafiles, BLANK = TRUE, IonMode = "Negative",
#' ion.id = c("Pos","Neg")) #Returns "Blanks_Neg"
#'
#' spos_db <- connect_peakdb(file.base = samples.pos, mem = TRUE)
#' sneg_db <- connect_peakdb(file.base = samples.neg, mem = TRUE)
#' bpos_db <- connect_peakdb(file.base = blanks.pos, mem = TRUE)
#' bneg_db <- connect_peakdb(file.base = blanks.neg, mem = TRUE)
#' new_db.list <- connect_lumadb(db.list = c("spos_db","sneg_db","bpos_db","bneg_db"), mem = TRUE)
#' all(sapply(new_db.list, function(x) { #All valid databases are created
#'   dbIsValid(x)
#' }))
#' all(sapply(new_db.list, dbDisconnect))
#' }
connect_lumadb = function(db.list, db.dir, new.db, mem) {

    #Set default variables
    if (missing(db.dir))
      db.dir = "db"
    if (missing(new.db))
        new.db = "Peaklist_db"
    if (missing(mem))
      mem = F

    if(mem) {

      peak_db <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
      pos_db <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
      blanks_pos_db <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
      neg_db <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
      blanks_neg_db <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")


    } else {


      peak_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(paste(db.dir, new.db, sep = "/"),".SQLite"))
      pos_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(paste(db.dir, db.list[[1]], sep = "/"),".SQLite"))
      blanks_pos_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(paste(db.dir, db.list[[3]], sep = "/"),".SQLite"))
      neg_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(paste(db.dir, db.list[[2]], sep = "/"),".SQLite"))
      blanks_neg_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(paste(db.dir, db.list[[4]], sep = "/"),".SQLite"))

    }


    return(list(peak_db = peak_db, pos_db = pos_db, neg_db = neg_db, blanks_pos_db = blanks_pos_db, blanks_neg_db = blanks_neg_db))

}

#' @title Reads Peak.list from database
#'
#' @export
#' @description Extracts Peak.list from an SQLite database as a tibble.
#' Alternatively, load the Peak.list into memory as a data frame
#' @param mytable character name of table in database to return
#' @param peak.db Formal class SQLiteConnection
#' @param asdf logical indicating whether to return a data frame instead of a tibble. Default is FALSE
#' @return tbl alternatively a data frame
#' @examples
#' library(LUMA)
#' if(require(RSQLite, quietly = TRUE)) {
#' file <- system.file("extdata/Sample_Data.csv", package =  "LUMA")
#' sample_data <- read.table(file, header = TRUE, sep = ",")
#' mzdatafiles <- sample_data$CT.ID
#' file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode =
#' "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
#' peak_db <- connect_peakdb(file.base = file.base, mem = TRUE)
#' dbIsValid(peak_db) #Database is valid
#' dbListTables(peak_db) #But no tables yet
#' mydf <- Peaklist_Pos$From_CAMERA
#' write_tbl(mydf = mydf, myname = "From_CAMERA", peak.db = peak_db)
#' test <- read_tbl(mytable = "From_CAMERA", peak.db = peak_db)
#' identical(test,mydf) ## Objects are identical of class "tbl"
#' test2 <- read_tbl(mytable = "From_CAMERA", peak.db = peak_db, asdf = TRUE)
#' !identical(test2,mydf) ## test2 is a data.frame
#' }
read_tbl = function(mytable, peak.db, asdf) {
    if (missing(asdf))
        asdf = FALSE
    if (asdf) {
        mydf <- dplyr::tbl(peak.db, mytable) %>% dplyr::collect() %>% data.frame
        return(mydf)
    } else {
        mytibble <- dplyr::tbl(peak.db, mytable) %>% dplyr::collect()
        return(mytibble)
    }
}

#' @title Writes Peak.list to database
#'
#' @export
#' @description Writes Peak.list to an RSQLite database. Peak.list can be a tibble or data frame.
#' @param mydf tbl or dataframe to write
#' @param peak.db Formal class SQLiteConnection
#' @param myname character what should the table be called
#' @return a tbl object in the remote source
#' @examples
#' library(LUMA)
#' if(require(RSQLite, quietly = TRUE)) {
#' file <- system.file("extdata/Sample_Data.csv", package =  "LUMA")
#' sample_data <- read.table(file, header = TRUE, sep = ",")
#' mzdatafiles <- sample_data$CT.ID
#' file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode =
#' "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
#' peak_db <- connect_peakdb(file.base = file.base, mem = TRUE)
#' dbIsValid(peak_db) #Database is valid
#' dbListTables(peak_db) #But no tables yet
#' mydf <- Peaklist_Pos$From_CAMERA
#' write_tbl(mydf = mydf, myname = "From_CAMERA", peak.db = peak_db)
#' dbListTables(peak_db) #Now has tables
#' dbDisconnect(peak_db)
#' }
write_tbl = function(mydf, peak.db, myname) {
    copy_to(peak.db, mydf, name = myname, temporary = FALSE, overwrite = TRUE)
}

#' @title Retrieves features from Peak.list in a database
#'
#' @description Returns mz/rt features from Peak.list stored in SQLite database
#' @param mytbl character name of table in database to return
#' @param peak.db Formal class SQLiteConnection
#' @param asdf logical indicating whether to return a data frame instead of a tibble. Default is FALSE
#' @return tbl alternatively a data frame
#' @examples
#' library(LUMA)
#' if(require(RSQLite, quietly = TRUE)) {
#' file <- system.file("extdata/Sample_Data.csv", package =  "LUMA")
#' sample_data <- read.table(file, header = TRUE, sep = ",")
#' mzdatafiles <- sample_data$CT.ID
#' file.base <- gen_filebase(mzdatafiles = mzdatafiles, BLANK = FALSE, IonMode =
#' "Positive", ion.id = c("Pos","Neg")) #Returns "Peaklist_Pos"
#' peak_db <- connect_peakdb(file.base = file.base, mem = TRUE)
#' dbIsValid(peak_db) #Database is valid
#' dbListTables(peak_db) #But no tables yet
#' mydf <- Peaklist_Pos$From_CAMERA
#' write_tbl(mydf = mydf, myname = "From_CAMERA", peak.db = peak_db)
#' test_tbl <- LUMA:::get_features(mytbl = "From_CAMERA", peak.db = peak_db)
#' test_tbl #Returns a tibble
#'
#' test_df <- LUMA:::get_features(mytbl = "From_CAMERA", peak.db = peak_db, asdf = TRUE)
#' test_df #Returns a data frame
#'
#' dbDisconnect(peak_db)
#' }
get_features = function(mytbl, peak.db, asdf) {
    if (missing(asdf))
        asdf = FALSE
    if (asdf) {
        mydf <- dplyr::tbl(peak.db, mytbl) %>%
          select(EIC_ID, mz, rt) %>%
          dplyr::collect() %>%
          data.frame
        return(mydf)
    } else {
        mytibble <- dplyr::tbl(peak.db, mytbl) %>%
          select(EIC_ID, mz, rt) %>%
          dplyr::collect()
        return(mytibble)
    }
}
