#' @title Filebase.gen
#'
#' @export
#' @description Generates filebase for reading and writing to databases
#' @param mzdatafiles character vector containing data files to process and store results in databases
#' @param BLANK a logical indicating whether blanks are being evaluated
#' @param ion.mode a character string defining the ionization mode.  Must be either "Positive" or "Negative"
#' @param ion.id character vector of length 2 specifying identifier in filename designating positive or negative ionization mode.  Positive identifier must come first.
#' @return character
filebase.gen=function(mzdatafiles,BLANK,ion.id,ion.mode) {
  if(ion.mode == "Positive" && BLANK == TRUE){
    mzdatafiles <- subset(mzdatafiles, subset = grepl(paste(ion.id[1]), mzdatafiles, ignore.case = TRUE))
    file.base = "Blanks_Pos"
    mzdatafiles <- mzdatafiles[c(grep("Blanks",mzdatafiles, ignore.case = TRUE))]
  } else {
    if(ion.mode == "Negative" && BLANK == TRUE){
      mzdatafiles <- subset(mzdatafiles, subset = grepl(paste(ion.id[1]), mzdatafiles, ignore.case = TRUE))
      file.base = "Blanks_Neg"
      mzdatafiles <- mzdatafiles[c(grep("Blanks",mzdatafiles))]
    } else {
      if(ion.mode == "Positive" && BLANK == FALSE){
        mzdatafiles <- subset(mzdatafiles, subset = grepl(paste(ion.id[1]), mzdatafiles, ignore.case = TRUE))
        file.base = "Peaklist_Pos"
        mzdatafiles <- mzdatafiles[-c(grep("Blanks",mzdatafiles))]
      } else {
        if(ion.mode == "Negative" && BLANK == FALSE){
          mzdatafiles <- subset(mzdatafiles, subset = grepl(paste(ion.id[1]), mzdatafiles, ignore.case = TRUE))
          file.base = "Peaklist_Neg"
          mzdatafiles <- mzdatafiles[-c(grep("Blanks",mzdatafiles))]
        } else {
          cat("Ion mode must be Positive or Negative.\nBe sure to specify whether to analyze blanks with logical indicator")
          stop()
        }
      }

    }

  }
  return(file.base)
}


#' @title peakdbConnect
#'
#' @export
#' @description Establishes a connection to an RSQLite database for storing data; if doesn't exist, creates new database
#' @param file.base character return from filebase.gen function
#' @param db.dir character what should the database directory be called.  Default is "db"
#' @return Formal class SQLiteConnection
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
peakdbConnect=function(file.base,db.dir) {
  if(missing(db.dir)) db.dir = "db"
  peak_db_file <- paste(file.base,db.dir,sep = "_")
  dir.create(db.dir, recursive = FALSE, showWarnings = FALSE)
  peak_db <- DBI::dbConnect(RSQLite::SQLite(),  paste(db.dir,peak_db_file, sep = "/"))
  return(peak_db)
}

#' @title libdbConnect
#'
#' @export
#' @description Establishes a connection to an RSQLite database for library searching; if doesn't exist, creates new database
#' @param lib.db character name of database
#' @param db.dir character directory containing the database
#' @return Formal class SQLiteConnection
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
libdbConnect=function(lib.db,db.dir) {
  lib_db <- DBI::dbConnect(RSQLite::SQLite(),  paste(db.dir,lib.db, sep = "/"))
  return(lib_db)
}

#' @title LUMA_dbConnect
#'
#' @export
#' @description Establishes a connection to an RSQLite database for combining two datasets together from two different ionization modes
#' @param db.list list chracter names of databases containing results from processing positive mode (1,3) and negative mode (2,4) data for samples (1,2) and blanks (3,4)
#' @param db.dir character directory containing the databases
#' @param new.db character what should the new database be called
#' @return list of Formal class SQLiteConnections, starting with new.db entry followed by one for each db.list entry and
#' @importFrom DBI dbConnect
#' @importFrom RSQLite SQLite
LUMA_dbConnect=function(db.list,db.dir,new.db) {
  if(missing(new.db)) new.db = "Peaklist_db"
  peak_db <- DBI::dbConnect(RSQLite::SQLite(), paste(db.dir,new.db, sep = "/"))
  pos_db <- DBI::dbConnect(RSQLite::SQLite(), paste(db.dir,db.list[[1]],sep = "/"))
  blanks_pos_db <- DBI::dbConnect(RSQLite::SQLite(), paste(db.dir,db.list[[3]],sep = "/"))
  neg_db <- DBI::dbConnect(RSQLite::SQLite(), paste(db.dir,db.list[[2]],sep = "/"))
  blanks_neg_db <- DBI::dbConnect(RSQLite::SQLite(), paste(db.dir,db.list[[4]],sep = "/"))
  return(list(peak_db = peak_db,pos_db=pos_db,neg_db=neg_db,blanks_pos_db=blanks_pos_db,blanks_neg_db=blanks_neg_db))
}

#' @title Readtbl
#'
#' @export
#' @description Extract table from an RSQLite database as a tibble.  Alternatively load into memory as a data frame
#' @param myname character name of table in database to return
#' @param peak.db Formal class SQLiteConnection
#' @param asdf logical indicating whether to return a data frame instead of a tibble. Default is FALSE
#' @return tbl alternatively a data frame
Readtbl=function(myname,peak.db,asdf) {
  if(missing(asdf)) asdf = FALSE
  if(asdf) {
    mydf <- dplyr::tbl(peak.db, myname) %>%
      dplyr::collect() %>%
      data.frame
    return(mydf)
  } else {
    mytibble <- dplyr::tbl(peak.db, myname) %>%
      dplyr::collect()
    return(mytibble)
  }
}

#' @title Writetbl
#'
#' @export
#' @description Writes tbl or dataframe to an RSQLite database
#' @param mytbl tbl or dataframe to write
#' @param peak.db Formal class SQLiteConnection
#' @param myname character what should the table be called
#' @return a tbl object in the remote source
Writetbl=function(mytbl,peak.db,myname) {
  copy_to(peak.db, mytbl, name = myname, temporary = FALSE, overwrite = TRUE)
}

#' @title GetFeatures
#'
#' @description Returns mz/rt features from a RSQLite database
#' @param myname character name of table in database to return
#' @param peak.db Formal class SQLiteConnection
#' @param asdf logical indicating whether to return a data frame instead of a tibble. Default is FALSE
#' @return tbl alternatively a data frame
GetFeatures=function(myname,peak.db,asdf) {
  if(missing(asdf)) asdf = FALSE
  if(asdf) {
    mydf <- dplyr::tbl(peak.db, myname) %>%
      select(EIC_ID,mz,rt) %>%
      dplyr::collect() %>%
      data.frame
    return(mydf)
  } else {
    mytibble <- dplyr::tbl(peak.db, myname) %>%
      select(EIC_ID,mz,rt) %>%
      dplyr::collect()
    return(mytibble)
  }
}
