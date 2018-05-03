#' XCMS Wrapper
#'
#' Run XCMS with user defined input parameters and return xcms objects
#' @param mzdatafiles a character vector of data files with full path names
#' @param XCMS.par a single-row data frame with 13 variables containing XCMS parameters. The column names must be c("Peakwidth1","Peakwidth2","ppm","noise","snthresh","mzdiff","prefilter1","prefilter2","center","gapInit","bw","mzwid","minfrac")
#' @return two XCMS objects xset and xset4 without and with retention time alignment, peak grouping, and imputing missing values
#' @export
xcmsWrap=function(mzdatafiles,XCMS.par) {
  xset <- xcmsSet(files = mzdatafiles ,
                  method="centWave",
                  peakwidth=c(XCMS.par$Peakwidth1, XCMS.par$Peakwidth2),
                  ppm=XCMS.par$ppm,
                  noise=XCMS.par$noise,
                  snthresh=XCMS.par$snthresh,
                  mzdiff=XCMS.par$mzdiff,
                  prefilter=c(XCMS.par$prefilter1, XCMS.par$prefilter2),
                  mzCenterFun="wMean",
                  integrate=1,
                  fitgauss=FALSE,
                  verbose.columns=FALSE,
                  BPPARAM = SnowParam(workers = snowWorkers(), type = "SOCK", stop.on.error = TRUE,
                                      progressbar = TRUE))
  pdf(file = paste(file.base, "RTDev Plot.pdf", sep = "_"))
  xset2 <- retcor(xset, method="obiwarp",
                  plottype="deviation",
                  distFunc="cor_opt",
                  profStep=1,
                  center=XCMS.par$center,
                  response=1,
                  gapInit=XCMS.par$gapInit,
                  gapExtend=2.7,
                  factorDiag=2,
                  factorGap=1,
                  localAlignment=0)
  dev.off()
  xset3 <- group(xset2, method="density",
                 bw=XCMS.par$bw,
                 mzwid=XCMS.par$mzwid,
                 minfrac=XCMS.par$minfrac,
                 minsamp=1, max=50)

  xset4 <- fillPeaks(xset3, BPPARAM = SnowParam(workers = snowWorkers(), type = "SOCK", stop.on.error = TRUE,
                                                progressbar = TRUE))
  return(xset,xset4)
}

#' CAMERA Wrapper
#'
#' Run CAMERA with user defined input parameters and return xsAnnotate objects
#' @param xset4 an xcms object that has had peak picking, retention time alignment, peak grouping, and imputing missing values performed
#' @param CAMERA.par a single-row data frame with 9 variables containing CAMERA parameters. The column names must be c("perfwhm","sigma","minfrac","mzabs","maxiso","corval_eic","corval_exp","pval","mzabs.1")
#' @param ion.mode a character string defining the ionization mode.  Must be either "positive" or "negative"
#' @return two grouped xsannotate objects mz1setpos and anposGa without and with annotated isotopes and ion adducts and fragments
#' @export
CAMERAWrap = function(xset4,CAMERA.par,ion.mode) {
  best.perfwhm <- CAMERA.par$perfwhm
  best.sigma <- CAMERA.par$sigma
  best.mzabs.iso <- CAMERA.par$mzabs
  best.minfrac <- CAMERA.par$minfrac
  best.maxiso <- CAMERA.par$maxiso
  best.corval_exp <- CAMERA.par$corval_exp
  best.corval_eic <- CAMERA.par$corval_eic
  best.pval <- CAMERA.par$pval
  best.mzabs.add <- CAMERA.par$mzabs.1

    mz1setpos <- xsAnnotate(xs = xset4, sample = NA)
  mz1setpos <- groupFWHM(object = mz1setpos, perfwhm = best.perfwhm, sigma = best.sigma , intval = "into")
  mz1setposIso <- findIsotopes(mz1setpos,maxcharge = 2,maxiso = best.maxiso, ppm = 3,mzabs = best.mzabs.iso, intval = "into",minfrac = best.minfrac,filter = TRUE)
  mz1setposIso.new <- groupCorr(object = mz1setposIso, cor_eic_th = best.corval_eic,
                                cor_exp_th = best.corval_exp, graphMethod = graph_method,
                                pval = best.pval, calcIso = TRUE, calcCiS = TRUE, calcCaS = TRUE,
                                psg_list = NULL, xraw = NULL)
  anposGa   <- findAdducts(mz1setposIso.new, polarity=CAMERA.ion.mode, ppm=3,
                           mzabs = best.mzabs.add, multiplier = 2, rules = rules)
  peakGa <- getPeaklist(object = anposGa)
  EIC_ID<-row.names(peakGa)
  peak_data <- cbind(EIC_ID, peakGa)
  return(mz1setpos,anposGa)
}

#' Compare CAMERA isotope annotation with user-defined annotation library
#' @param Peaklist a table of class "tbl_dbi", "tbl_sql", "tbl_lazy", or "tbl" with samples as columns.  Should contain all output columns from XCMS and CAMERA, both metadata and sample data. Retention times must be in min.
#' @param Annotated.library a data frame with annotated metabolite entries. Must contain columns called "Name", "Formula", Molecular.Weight" and "RT..Min.".  Can contain additional info as separate columns.
#' @param rules a data frame containing the rule list used by CAMERA to annotate ion adducts and fragments.  Must contain the columns "name","nmol","charge","massdiff","oidscore","quasi","ips".
#' @param search.par a single-row data frame with 11 variables containing user-defined search parameters. Must contain the columns "ppm","rt","Voidrt","Corr.stat.pos","Corr.stat.neg","CV","Minfrac","Endogenous","Solvent","gen.plots","keep.singletons".
#' @param ion.mode a character string defining the ionization mode.  Must be either "Positive" or "Negative"
#' @param lib_db RSQLite connection
#' @return data frame containing the original table with added columns "Name","MS.ID","Formula","Annotated.adduct" and any additional info columns from Annotated.Library
#' @export
IHL.search=function(Peaklist,Annotated.library,rules,search.par,ion.mode,lib_db)  {
  search.list <- Peak.list %>%
    select(EIC_ID,mz,rt) %>%
    dplyr::collect()

  ## Creates full adduct list for all compounds in the Annotated library
  IHL <- Annotated.Library[rep(seq_len(nrow(Annotated.Library)), each = nrow(rules)),]
  x <- rules$nmol
  IHL.temp <- sweep(IHL,1,x,"*")
  x <- rules$massdiff
  if(ion.mode == "Positive"){
    IHL.temp <- sweep(IHL.temp,1,x,"+")
    Ion.Mode <- "Pos"
    bin <- paste(Ion.Mode,"_",search.list$EIC_ID,"_", sep = "")

  } else {
    if(ion.mode == "Negative"){
      IHL.temp <- sweep(IHL.temp,1,x,"+")
      IHL.temp <- sweep(IHL.temp,1,-1,"*")
      Ion.Mode <- "Neg"
      bin <- paste(Ion.Mode,"_",search.list$EIC_ID,"_", sep = "")

    } else {
      break("You must include the ionization mode!")
    }
  }
  x <- rules$charge
  IHL.adduct.data <- sweep(IHL.temp,1,x,"/")
  IHL[,'mz'] <- IHL.adduct.data$Molecular.Weight
  IHL[,'adduct'] <- rules$name
  copy_to(lib_db, IHL,name = paste("Annotated Library",Ion.Mode,sep = "_"), temporary = FALSE, overwrite = TRUE)
  rm(IHL,IHL.adduct.data,IHL.temp)
  bin
  IHL <- tbl(lib_db, paste("Annotated Library",Ion.Mode,sep = "_"))
  d.mz <- search.list$mz*as.numeric(search.par[1,"ppm"])/10^6
  d.rt <- as.numeric(search.par[1,"rt"])

  #calculates the min and max range values for searching the in house library
  search.list$mz.min <- search.list$mz - d.mz
  search.list$mz.max <- search.list$mz + d.mz
  search.list$rt.min <- search.list$rt - d.rt
  search.list$rt.max <- search.list$rt + d.rt
  search.list$MS.ID = NA
  search.list$Formula = NA
  search.list$Name = NA
  search.list$Annotated.adduct = NA
  search.list$Conf.Level = NA
  search.list$FISh.Coverage = NA

  ## attempts to match all peaks against the In House Library
  #! Try to build a query using *apply functions to search all of the search list at once; should speed
  #! things up considerably
  # i = 27  ##Used for debugging purposes
  total = nrow(search.list)
  pb = winProgressBar(title = "Annotating Metabolites.", min = 0, max = total, width = 300)
  counter = 1
  for(i in 1:nrow(search.list)) {
    mz.min = search.list$mz.min[i]
    mz.max = search.list$mz.max[i]
    rt.min = search.list$rt.min[i]
    rt.max = search.list$rt.max[i]
    test.list <- IHL %>%
      filter(
        between(mz,mz.min,mz.max)) %>%
      filter(
        between(RT..min.,rt.min,rt.max)) %>%
      dplyr::collect()
    if (nrow(test.list) == 0) {
      search.list$MS.ID[i] = paste(bin[i],"Unidentified",sep = "")
    } else {
      if(nrow(test.list) >= 1 ) {
        cat("\n\n\nMetabolite annotated for match number ", counter, ".\nMatch results below.\n\n\n", sep = "")
        str(test.list)
        search.list$MS.ID[i] = paste(bin[i],"Annotated",sep = "")
        search.list$Name[i] = glue::collapse(test.list$Name, sep = ";", width = Inf, last = " or ")
        search.list$Formula[i] = glue::collapse(unique(gsub(" ", "",test.list$Formula)), sep = ";", width = Inf, last = " or ")
        search.list$Annotated.adduct[i] = glue::collapse(test.list$adduct, sep = ";", width = Inf, last = " or ")
        search.list$Conf.Level[i] = glue::collapse(test.list$Levels, sep = ";", width = Inf, last = " or ")
        search.list$FISh.Coverage[i] = glue::collapse(test.list$Fish.Coverage, sep = ";" , width = Inf, last = " or ")
        counter = counter + 1
      }
    }

    setWinProgressBar(pb,i, title=paste("Annotating Metabolites: ",round(i/total*100,0),"% done", sep = ""))
  }
  close(pb)
  named.peak.list <- search.list[,8:ncol(search.list)]
  colnames(named.peak.list)
  temp <- as.data.frame(Peak.list)
  colnames(Peak.list)
  if(nrow(temp) == nrow(named.peak.list)){
    temp <- cbind(Peak.list,named.peak.list)
  }
  colnames(temp)
  return(temp)
}

#' Calculate Minfrac for sample classes
#'
#' Calculates the fraction of samples within each class that contains a given feature.  Only the minimum fraction across all sample classes is returned.
#' @param temp a data frame with class info as columns.  Must contain a separate row entry for each unique sex/class combination. Must contain the columns "Sex","Class","n","Endogenous".
#' @param xset4 an xcms object that has had peak picking, retention time alignment, peak grouping, and imputing missing values performed
#' @param BLANK a logical indicating whether blanks are being evaluated
#' @param Peak.list a table of class "tbl_df",tbl" or "data.frame" with variables as columns.  Should contain all output columns from XCMS and CAMERA, and additional columns from IHL.search.
#' @return data frame containing the original table with one additional column "Minfrac" at the end, followed by the CAMERA columns "isotopes","adduct","pcgroup"
#' @export
Calc.minfrac=function(temp,xset4,BLANK,Peak.list)  {
  peakSN <- peakTable(xset4, filebase="SN_Livers_All classes", value="sn") #writes the SN peak table to file
  SN.list<-read.table(file = "SN_Livers_All classes.tsv", sep="\t", header = TRUE) #reads the SN peak table
  file.remove("SN_Livers_All classes.tsv")

    if(BLANK == TRUE) {} else {
    sexes <- unique(paste(temp$Sex,"_",sep = ""))
    samples <- vector(mode = "character", length = length(colnames(SN.list)))
    for(i in 1:length(sexes)) {
      rows_loop <- grep(sexes[i],colnames(SN.list))
      samples[rows_loop] <- sexes[i]
    }
    sum.range.list<- SN.list[,samples %in% sexes]
    t.list<-t(sum.range.list) #transposes the SN list
    write.table(t.list,file = "SN_transposed.csv",sep = ",",row.names=FALSE)
    t.list<-read.table(file = "SN_transposed.csv", sep=",", header = TRUE) #reads the transposed SN peak table
    file.remove("SN_transposed.csv")

    ## Grabs the unique sample ID number for each sample, if present in the filename
    bin <- variable.names(t.list, full=TRUE)
    colnames(sum.range.list)
    sample.ID <- sub('\\D*(\\d{6}).*', '\\1',colnames(sum.range.list))
    colnames(sum.range.list)

    ## Creates a new column for grouping by class based on user input
    groups <- paste(temp$Sex,temp$Class,sep = "_")
    group <- vector(mode = "character", length = length(colnames(sum.range.list)))
    for(i in 1:length(groups)) {
      rows_loop <- grep(groups[i],colnames(sum.range.list))
      group[rows_loop] <- groups[i]
    }
    class.n <- unique(temp$n)

    if(length(class.n) > 1) {
      t.list<-cbind(group,t.list) #combines the two dataframes together
      sn.list<- split(t.list,as.factor(group))  #Splits the data frame by a factor
      sn.list.split <- split(sn.list, sapply(sn.list, function(x) nrow(x)))

      #Create matrix to contain all of the minfrac values for each group from each list in sn.list
      min.frac.fulltable <- matrix(nrow = length(bin), ncol = sum(sapply(sn.list, function(x) length(x))))
      # j = 1 This and following line used for debugging purposes
      # k = 1
      for(j in 1:length(sn.list)){
        na_count <- sapply(as.data.frame(sn.list[j]), function(x,y) sum(length(which(is.na(x)))), y=group) #counts the number of na values per feature across each group
        na_count <- data.frame(na_count) #puts the list in a data frame
        for(i in 1:length(groups)) {
          rows_loop <- grep(groups[i],rownames(na_count))
          group[rows_loop] <- groups[i]
        }
        na_count<-cbind(group,na_count) #combines the grouping variable with the na count values
        na.count.list<- split(na_count,as.factor(group))  #Splits the data frame into lists of data frames by the grouping variable
        str(na.count.list)

        all_count <- sapply(as.data.frame(sn.list[j]), function(x,y) length(x), y=group) #counts the number of na values per feature across each group
        all_count <- data.frame(all_count) #puts the list in a data frame
        str(all_count)
        for(i in 1:length(groups)) {
          rows_loop <- grep(groups[i],rownames(all_count))
          group[rows_loop] <- groups[i]
        }
        all_count<-cbind(group,all_count) #combines the grouping variable with the na count values
        all.count.list<- split(all_count,as.factor(group))  #Splits the data frame into lists of data frames by the grouping variable
        str(all.count.list)

        ##The next section creates a new data matrix with the number of rows and columns equal to the number of features and groups, respectively.  The new table is called 'minfrac.fulltable_loop'
        minfrac.fulltable_loop <- matrix(nrow = length(bin), ncol = length(all.count.list))
        colnames(minfrac.fulltable_loop)<-c(names(na.count.list))###Names the columns for the new table as indicated
        rownames(minfrac.fulltable_loop)<-(bin)###Pulls out and names the table rows based on the metabolite label (I do this as a QA to make sure the loop is correctly linking the ANOVA output and metbaolite IDs

        ##this section calculates the minfrac values for each exposure class in a for loop and populates the minfrac.fulltable_loop with these values
        group.names <- names(na.count.list)
        for (i in names(na.count.list)){
          na_loop <- na.count.list[i]
          na_loop <- do.call("rbind", na_loop)
          all_loop <- all.count.list[i]
          all_loop <- do.call("rbind", all_loop)
          min_frac_loop<-rapply(as.data.frame(na_loop$na_count), function(x,y) x/y, y=all_loop$all_count)##calculates the min frac values for each feature by exposure class
          min_frac_loop <- data.frame(min_frac_loop)##changes the format of the min frac results into a data frame
          colnames(min_frac_loop) <- paste(i,"","min_frac")
          min_frac_loop <- min_frac_loop[-1,]##Removes the value corresponding to the grouping variable
          minfrac.fulltable_loop[,i] = min_frac_loop##Pastes the min frac results from an exposure class into the next col of the min frac full table.
        }
        ##Code to calculate min_frac for one of the exposure classes; used to QA the loop output
        i=names(na.count.list[1])
        group.names <- names(na.count.list)
        test.na <- na.count.list[i]
        test.na<-do.call("rbind", test.na)
        test.all <- all.count.list[i]
        test.all<-do.call("rbind", test.all)
        min_frac<-rapply(as.data.frame(test.na$na_count), function(x,y) x/y, y=as.data.frame(test.all$all_count))#Need to find a way to get rid of the subset operator $
        min_frac <- data.frame(min_frac)
        str(min_frac)
        colnames(min_frac) <- paste(i,"","min_frac")
        min_frac <- min_frac[-1,]
        head(min_frac)

        #Bring loop-based minfrac values out of the loop by placing them in the minfrac master table
        min.frac.fulltable[,k:sum(k-1,ncol(minfrac.fulltable_loop))] <- minfrac.fulltable_loop

        #Update the counter to keep track of where to copy the results in the master table
        k <- k + ncol(minfrac.fulltable_loop)
      }
    } else {
      if(length(class.n) == 1) {
        t.list<-cbind(group,t.list) #combines the two dataframes together
        sn.list<- split(t.list,as.factor(group))  #Splits the data frame by a factor

        na_count <- sapply(as.data.frame(sn.list), function(x,y) sum(length(which(is.na(x)))), y=group) #counts the number of na values per feature across each group
        na_count <- data.frame(na_count) #puts the list in a data frame
        # group<-unlist(strsplit(gsub("(([MF])_([RN])S_T([0482]))|.", "\\1", rownames(na_count)), "\\s+")) #generates a character vector of grouping variables for each na count value
        for(i in 1:length(groups)) {
          rows_loop <- grep(groups[i],rownames(na_count))
          group[rows_loop] <- groups[i]
        }
        na_count<-cbind(group,na_count) #combines the grouping variable with the na count values
        na.count.list<- split(na_count,as.factor(group))  #Splits the data frame into lists of data frames by the grouping variable
        all_count <- sapply(as.data.frame(sn.list), function(x,y) length(x), y=group) #counts the number of na values per feature across each group
        all_count <- data.frame(all_count) #puts the list in a data frame
        for(i in 1:length(groups)) {
          rows_loop <- grep(groups[i],rownames(all_count))
          group[rows_loop] <- groups[i]
        }
        all_count<-cbind(group,all_count) #combines the grouping variable with the na count values
        all.count.list<- split(all_count,as.factor(group))  #Splits the data frame into lists of data frames by the grouping variable

        ##The next section creates a new data matrix with the number of rows and columns equal to the number of features and groups, respectively.  The new table is called 'min.frac.fulltable'
        min.frac.fulltable <- matrix(nrow = length(bin), ncol = length(all.count.list))
        colnames(min.frac.fulltable)<-c(names(na.count.list))###Names the columns for the new table as indicated
        rownames(min.frac.fulltable)<-(bin)###Pulls out and names the table rows based on the metabolite label (I do this as a QA to make sure the loop is correctly linking the ANOVA output and metbaolite IDs
        head(min.frac.fulltable)##A double-check to make sure the new table was created properly

        #this section calculates the minfrac values for each exposure class in a for loop and populates the min.frac.fulltable with these values
        group.names <- names(na.count.list)
        for (i in names(na.count.list)){
          na_loop <- na.count.list[i]
          na_loop <- do.call("rbind", na_loop)
          all_loop <- all.count.list[i]
          all_loop <- do.call("rbind", all_loop)
          min_frac_loop<-rapply(as.data.frame(na_loop$na_count), function(x,y) x/y, y=all_loop$all_count)##calculates the min frac values for each feature by exposure class
          min_frac_loop <- data.frame(min_frac_loop)##changes the format of the min frac results into a data frame
          colnames(min_frac_loop) <- paste(i,"","min_frac")
          min_frac_loop <- min_frac_loop[-1,]##Removes the value corresponding to the grouping variable
          min.frac.fulltable[,i] = min_frac_loop##Pastes the min frac results from an exposure class into the next row of the min frac full table.
        }
        #Code to calculate min_frac for one of the exposure classes; used to QA the loop output
        i=names(na.count.list[2])
        group.names <- names(na.count.list)
        test.na <- na.count.list[i]
        test.na<-do.call("rbind", test.na)
        test.all <- all.count.list[i]
        test.all<-do.call("rbind", test.all)
        min_frac<-rapply(as.data.frame(test.na$na_count), function(x,y) x/y, y=as.data.frame(test.all$all_count))#Need to find a way to get rid of the subset operator $
        min_frac <- data.frame(min_frac)
        colnames(min_frac) <- paste(i,"","min_frac")
        min_frac <- min_frac[-1,]
        head(min_frac)

      } else print("Please include the number of samples in each class in the Sample Class file!")
    }
    ##Code to calculate MinFrac as the minimum fraction of samples with a feature across all exposure classes and sexes
    MinFrac.table <- matrix(nrow = length(bin), ncol = 1) #This creates a new table to fill with one MinFrac value for each feature
    for (i in 1:length(bin)) {
      minfrac_loop <- 1 - min(min.frac.fulltable[i,], na.rm = TRUE)
      MinFrac.table[i,] <- minfrac_loop
    }
  }

  CAMERA.list <- Peak.list[,c("isotopes","adduct","pcgroup")]  #Extracts all of the CAMERA columns to a separate dataframe
  drops <- c("isotopes","adduct","pcgroup")
  stripped.list <- Peak.list[, !(names(Peak.list) %in% drops)] #Removes the CAMERA columns
  stripped.list[,"MinFrac"] <- MinFrac.table #Attaches the MInFrac column to the stripped down peak.list
  new.peak.list <- cbind(stripped.list,CAMERA.list)
  raw <- new.peak.list
return(raw)
}

#' CAMERA Parser for positive ionization mode
#'
#' Parses the CAMERA results using well-defined rules in positive ionization mode.
#' @param raw
#' @export
CAMERA.parser.pos=function(){
  ##This function is an alternative of CAMERA parser using MATLAB from the Ressom Omics Lab at Georgetown University
  ##*******************************************************
  ##*******************POSITIVE   MODE*********************
  ##*******************************************************
  ##Parser of the CAMERA output results to assist analysis of ion annotation
  ##Note: the monoisotpic mass means the mass of [M+H]/[M-H] type of ion with the mono isotopes.
  # setwd("C:/Users/jmosley/Desktop/Maumee River ddMS/Data Processing/Scripts")
  # dir()
  # getwd()
  # rm(list = ls())
  # ls()
  # setwd("C:/Users/jmosley/Desktop/Maumee River ddMS/Data Processing/Scripts")
  # home<-setwd("C:/Users/jmosley/Desktop/Maumee River ddMS/Data Processing/Scripts")
  library(openxlsx)

  ## Input CAMERA output and positive Rule file
  # raw <- tbl(peak_db, "From CAMERA with MinFrac") %>%
  #   collect()
  # raw <- read.table(paste(file.base,"From CAMERA with MinFrac.csv", sep = "_"), sep=",", header = TRUE)
  raw.header <- colnames(raw)
  # rule <- read.table(file = "CAMERA_POS_Rules.csv", sep=",", header= TRUE)
  # rule <- read.table(file = "CAMERA_POS_Rules_modi.csv", sep=",", header= TRUE)
  rule <- read.table(file = "primary_adducts_pos.csv", sep=",", header= TRUE)
  rule[,"X"] <- as.numeric(row.names(rule))
  ##Move last column to first in rule dataframe
  rule <- rule[,c(8,1,2,3,4,5,6,7)]
  rule.header <- colnames(rule)

  PROTO_MASS <- 1.00727646677

  ##Read isotope, adduct, pcgroup and m/z information
  adducts<-which(raw.header == "adduct")
  isotopes<-which(raw.header == "isotopes")
  mz<-which(raw.header == "mz")
  pcgroup<-which(raw.header == "pcgroup")

  isotopes<-raw[isotopes]
  adducts<-raw[adducts]
  mz<-raw[mz]
  pcgroup<-raw[pcgroup]

  ##Check if there is two lines with same adduct annotation and same pcgroup
  ##Find same annotation in the same pcgroup, both of them are removed
  for (i in 1:nrow(adducts)){
    if(adducts[i,1]!=""){
      idx<-intersect(which(adducts$adduct %in% adducts[i,1]), which(pcgroup$pcgroup %in% pcgroup[i,1]))
      if (length(idx)>1){
        warning ("Find same annotation in the same pcgroup, both of them are removed: ",as.character(adducts[i,1]),", pcgroup ",pcgroup[i,1])
        adducts[idx,1]<-""
      }
      else if (length(idx)==0){
        stop("fatal error:1")
      }
    }
  }

  ##Alternative code for section above
  # x.temp<-cbind(adducts,pcgroup)
  # y.temp<-x.temp[duplicated(x.temp),]
  # z.temp<-unique(y.temp)
  # for (i in 1:nrow(z.temp)){
  #   if (z.temp$adduct[[i]]!="")
  #   {
  #     a<-which(x.temp$adduct == z.temp$adduct[[i]])
  #     for (k in 1:length(a)){
  #       if (x.temp$pcgroup[[a[k]]]==z.temp$pcgroup[[i]]){
  #         x.temp[a[k],1]<-""
  #       }
  #     }
  #   }
  # }

  ##Prepare a file with cleaned annotation information
  raw_cleaned <- raw
  raw_cleaned[,length(raw)-1] <- adducts
  raw<-raw_cleaned
  pcgroup_no <- sort(unique(raw_cleaned$pcgroup))

  ##Find and split those ions with more than one annotations
  # raw_single_annotation <- matrix(nrow =nrow(raw), ncol = ncol(raw))
  ambiguity_flg <- matrix(nrow = 0, ncol = 1)
  adducts$adduct <- as.character(adducts$adduct)
  for (j in 1:length(adducts$adduct)){
    if (length(grep("\\[.+\\].+\\[.+\\]", toString(adducts$adduct[[j]])))>0){
      warning ("Conflicting annotations are splitted in pc_group ",pcgroup[j,],": ",as.character(adducts[j,]),":::",j)
      single_annotations <- strsplit(toString(adducts$adduct[[j]]),"(?<=[^-]\\d{1}) ", perl=TRUE)
      # % If the multiple adduct annotation coincide with isotope
      # % annotation, the isotope annotation is only preserved for the first
      # % adduct annotation. Otherwise it will cause error. One possible
      # % solution in the future is to add new isotope group and split
      # % related isotope annotations.
      #! Try to include some code to check if one of the conflicting annotations is present in the
      #! Matched.adduct column
      single_annotations <- unique(single_annotations[[1]])
      score<-1000
      a<-0
      b<-0
      for (i in 1:length(single_annotations)){
        temp <- sapply(strsplit(toString(single_annotations[i])," "), "[", 1)##isolate the adducts
        idx<-which(rule$name %in% temp)
        if (rule$X[idx]<score & rule$ips[idx]>=b) {
          score<- rule$X[idx]
          a<-i
          b<-rule$ips[idx]
        }
      }
      adducts$adduct[j] <- single_annotations[a]
      #     single_line <- raw_cleaned[j,]
      #     single_line[] <- lapply((single_line), as.character)
      #     single_line$adduct[[1]]<-single_annotations[a]
      #     raw_single_annotation<-rbind(raw_single_annotation,single_line[1,])
      ambiguity_flg<-rbind(ambiguity_flg,1)
    }
    else{
      #     raw_single_annotation<-rbind(raw_single_annotation,raw_cleaned[j,])
      ambiguity_flg<-rbind(ambiguity_flg,0)
    }
  }

  raw_single_annotation<-raw_cleaned

  raw_single_annotation[,length(raw)-1] <- adducts


  # row.names(raw_single_annotation)<-1:nrow(raw_single_annotation)
  # copy_to(peak_db, raw_single_annotation, name = "duplicate adduct removal for positive")
  write.table(raw_single_annotation,file = "duplicate adduct removal for positive.csv",sep = ",",row.names=FALSE)
  # write.table(ambiguity_no,file="ambiguity_no.csv", sep = ",",row.names=FALSE)
  write.table(ambiguity_flg,file="ambiguity_flg.csv", sep = ",",row.names=FALSE)
  raw_single_annotation<-read.table(file = "duplicate adduct removal for positive.csv", sep=",", header = TRUE)
  # ambiguity_no<-read.table(file = "ambiguity_no.csv", sep=",", header = TRUE)
  ambiguity_flg<-read.table(file = "ambiguity_flg.csv", sep=",", header = TRUE)
  file.remove("duplicate adduct removal for positive.csv")
  file.remove("ambiguity_flg.csv")
  ##Re-extract the annotations
  isotopes <- raw_single_annotation[which(colnames(raw_single_annotation)=="isotopes")]
  adducts <- raw_single_annotation[which(colnames(raw_single_annotation)=="adduct")]
  mz <- raw_single_annotation[which(colnames(raw_single_annotation)=="mz")]
  pcgroup <- raw_single_annotation[which(colnames(raw_single_annotation)=="pcgroup")]
  pcgroup <- data.matrix(pcgroup)
  pcgroup_no<-sort(unique(pcgroup))

  ion_info<-matrix(c(0,0,0,0),ncol=4)
  ion_info_copy <- vector(mode = "numeric", length = 0)
  adduct_flg<-matrix(nrow = length(pcgroup), ncol = 1)
  iso_flg<-matrix(nrow = length(pcgroup), ncol = 1)
  mono_flg<-matrix(nrow = length(pcgroup), ncol = 1)
  ion_grp<-matrix(nrow = length(pcgroup), ncol = 1)
  iso_grp<-matrix(nrow = length(pcgroup), ncol = 1)
  mono_mass<-matrix(nrow = length(pcgroup), ncol = 1)

  for (i in 1:length(pcgroup_no)){
    curr_group<-pcgroup_no[i]
    ion_idx<-which(pcgroup %in% curr_group)
    group_iso<-isotopes[ion_idx,]
    group_add<-adducts[ion_idx,]
    group_mz<-as.numeric(mz[ion_idx,])
    ##Parse the adducts for a pcgroup. If it is [M+H], it is a monoisotopic ion. If it is annotated
    ##otherwise, it is an adduct ion. If there is no annotation, it is NOT an adduct ion at least,
    ##but not necessarily a monoisotopic ion. It should be noted if an ion is annotated as an
    ##adduct, only the one of mono isotopes is annotated. Although the adduct ion can have its
    ##set of isotopic ions. That's the reason why adducts are first processed here.
    for (j in 1:length(group_add)){
      if (group_add[j]!=""){
        curr_rule <- sapply(strsplit(toString(group_add[j])," "), "[", 1)##isolate the adducts
        rule_idx<-which(rule$name %in% curr_rule)
        ##mass_tmp is the equivalent m/z of [M+H] or [M-H] type of ion
        if(length(rule_idx)==0){
          stop("NO match found for adduct pattern ", as.character(curr_rule)," in rule list")
        }
        nmol<-rule[rule_idx,3]
        z<-rule[rule_idx,4]
        massdiff<-rule[rule_idx,5]
        mass_tmp<- (group_mz[j]*z - massdiff)/nmol +PROTO_MASS
        mass_signature<-as.numeric(sapply(strsplit(toString(group_add[j])," "), "[", 2))
        ##[M+H]or[M-H]
        if (rule_idx==1){
          mono_flg[ion_idx[j]]<-1
          mono_mass[ion_idx[j]]<-group_mz[j]
          adduct_flg[ion_idx[j]]<-0
        }
        else{
          adduct_flg[ion_idx[j]]<-1
          mono_flg[ion_idx[j]]<-0
        }
        ##If the ion has the same monoisotopic value and pcgroup as a previos ion, they are
        ##assigned to the same ion group with same monoisotopic mass. If not, a new ion_info
        ##entry is created
        if(any((mass_signature==ion_info[,4])&(curr_group==ion_info[,3]))){
          ion_grp_idx<-intersect(which(ion_info[,4] %in% mass_signature),which(ion_info[,3] %in% curr_group))
          ion_grp[ion_idx[j]]<-ion_info[ion_grp_idx,1]
        }
        else{
          ion_grp[ion_idx[j]]<-max(ion_info[,1])+1
          ion_info_temp <-cbind(max(ion_info[,1])+1,mass_tmp,curr_group,mass_signature)
          ion_info<-rbind(ion_info,ion_info_temp)
        }
      }
      else{
        adduct_flg[ion_idx[j]]<-0
      }
    }
    # %% Parse the isotopes of a pcgroup. If an ion is [M]+ and not an adduct of
    # %% the other ion, it is a monoisotopic ion, otherwise it is an adduct ion
    # %% with monoisotopic mass (but not considered as monoisotopic ion here). If
    # %% an ion is multiple-charged [M]n+(it is labeled under "isotopes" rather than
    # %% "adduct"), its monoisotopic mass with single charge is calculated as
    # %% (n*m/z - (n-1)*mass-of-Proton). If an ion is labeled but neither [M]+ nor
    # %% [M]n+, it is an isotopic ion. If an ion is not labeled, it is not an
    # %% isotopic ion.
    for (j in 1:length(group_iso)){
      if(group_iso[j]!=""){
        ## isolate the isotope group number
        iso_grp[ion_idx[j]]<-sapply(strsplit(toString(group_iso[j]),"\\[|\\]"), "[", 2)
        ## extract the isotopic information
        curr_iso<-sapply(strsplit(toString(group_iso[j]),"\\[\\d+\\]"), "[", 2)
        if(curr_iso=='[M]+'){
          if(adduct_flg[ion_idx[j]]==0){
            mono_flg[ion_idx[j]]<-1
            mono_mass[ion_idx[j]]<-group_mz[j]
            iso_flg[ion_idx[j]]<-0
          }
          else{
            mono_flg[ion_idx[j]]<-0
            iso_flg[ion_idx[j]]<-0
          }
        }
        ## if an ion is multiple-charged [M]n+
        else if(length(grep(".+\\[\\M\\]\\d\\+",toString(group_iso[j])))>0){
          charge_state<-as.numeric(sapply(strsplit(toString(group_iso[j]),"\\M\\]|\\+"), "[", 2))
          mono_flg[ion_idx[j]]<-0
          if(adduct_flg[ion_idx[j]]!=1){
            # new mono_mass is assigned only when this ion has no adduct annotation,
            # otherwise the mono_mass is calculated based on adduct annotation to deal
            # with clusterion (such as [2M+2H]2+)
            mono_mass[ion_idx[j]]<-charge_state*group_mz[j]-(charge_state-1)*PROTO_MASS
          }
          iso_flg[ion_idx[j]]<-0
        }
        else {
          iso_flg[ion_idx[j]]<-1
          mono_flg[ion_idx[j]]<-0
        }
      }
      else{
        iso_flg[ion_idx[j]]<-0
      }
    }
  }

  # %% For ions which are neither adducts nor isotopes and has not been assigned
  # %% a monoisotopic m/z yet, they are assumed to be monoisotopic and the
  # %% monoisotopic m/z are assigned.
  for (i in 1:length(pcgroup)){
    if (adduct_flg[i]==0 && iso_flg[i]==0 && is.na(mono_mass[i])){
      mono_flg[i]<-1
      mono_mass[i]<-mz$mz[i]
    }
  }
  ion_info<-ion_info[2:nrow(ion_info),]
  # %% If multiple ions are in the same ion group, they are adducts for each
  # %% other and share the same monoisotopic mass. If one of them is previously
  # %% assigned a mono isotopic mass, use the assigned one. If none of them has
  # %% monoisotopic mass (e.g. they are [M+Na] and [M+K]), use the one
  # %% calculated by CAMERA and stored in ion_info.
  temp_idx <- 1
  for (k in 1:nrow(ion_info)){
    #! if(is.matrix(ion_info_copy)){
    #!   ion_info_copy <- ion_info_loop
    #! }
    grp_ion_idx<-which(ion_grp %in% ion_info[k,1])
    grp_mass<-mono_mass[grp_ion_idx]
    if (sum(!is.na(grp_mass))==1){
      mono_mass[grp_ion_idx]<-grp_mass[!is.na(grp_mass)]
    }
    else if ((sum(!is.na(grp_mass))>1)&&(sum(mono_flg[grp_ion_idx])==1)){
      # % If there are both [M+2H]2+ and [M]2+ type of annotation, that
      # % means there are both adducts and isotopes of this ion. In
      # % this case, there could be two valid monoisotopic mass from
      # % previous calculation. One from [M+2H] and the other from
      # % [M+H]([M+Na] will not give a mono mass). In this case, we use
      # % the monoisotopic m/z from [M+H]+
      mono_mass[grp_ion_idx]<-grp_mass[mono_flg[grp_ion_idx]==1]
    }
    else if (sum(!is.na(grp_mass))>=2){
      warning("More than one mono mass value for adducts. \nRemoving conflicting adduct annotation.")
      print(raw_single_annotation[grp_ion_idx,(ncol(raw_single_annotation)-2):ncol(raw_single_annotation)])
      if(is.vector(ion_info_copy)){
        ion_info_copy <- ion_info[-k,]
        for (l in 1:length(grp_ion_idx)){
          mass_tmp <- raw_single_annotation[grp_ion_idx[l],"mz"]
          mass_signature <- mass_tmp + PROTO_MASS
          curr_group <- raw_single_annotation[grp_ion_idx[l],"pcgroup"]
          ion_info_temp <-cbind(max(ion_info[,1])+temp_idx,mass_tmp,curr_group,mass_signature)
          ion_grp[grp_ion_idx[l]] <- ion_info_temp[,1]
          ion_info_copy <- rbind(ion_info_copy,ion_info_temp)
          temp_idx <- temp_idx + 1
        }
      } else {
        if(nrow(ion_info_copy) >= nrow(ion_info) +1){
          sub.idx <- nrow(ion_info_copy) - nrow(ion_info)
          ion_info_loop <- ion_info_copy[-(k-sub.idx),]
          for (l in 1:length(grp_ion_idx)){
            mass_tmp <- raw_single_annotation[grp_ion_idx[l],"mz"]
            mass_signature <- mass_tmp - PROTO_MASS
            curr_group <- raw_single_annotation[grp_ion_idx[l],"pcgroup"]
            ion_info_temp <-cbind(max(ion_info[,1])+temp_idx,mass_tmp,curr_group,mass_signature)
            ion_grp[grp_ion_idx[l]] <- ion_info_temp[,1]
            ion_info_loop <- rbind(ion_info_loop,ion_info_temp)
            temp_idx <- temp_idx + 1
          }
          ion_info_copy <- ion_info_loop
        }
      }
      # stop("More than one mono mass value for adducts")
    }
    else{
      mono_mass[grp_ion_idx] = ion_info[k, 2]
    }
  }
  # %% If multiple ions are in the same isotope group, they are isotopes for
  # %% each other and have the same monoisotopic mass. For one isotopic ions to
  # %% appear, its monoisotopic ion must present.

  iso_grp_no<- unique(iso_grp[!is.na(iso_grp)])

  for (k in 1:length(iso_grp_no)){
    iso_grp_idx<-which(iso_grp %in% iso_grp_no[k])
    iso_grp_mass<-mono_mass[iso_grp_idx]
    if(sum(!is.na(iso_grp_mass))==1){
      mono_mass[iso_grp_idx]<-iso_grp_mass[!is.na(iso_grp_mass)]
    }else if(sum(!is.na(iso_grp_mass))>1) {
      stop("more than one mono mass value for isotopes for isotope group ",iso_grp_no[k])
    }else {
      stop("there is no monoisotopic ion found for isotope group ", iso_grp_no[k])
    }
  }

  ##save(pcgroup,iso_grp_no,iso_grp_no,pcgroup_no,group_add,group_iso,adduct_flg,ambiguity_flg,iso_flg,isotopes,mono_flg,mono_mass,ion_info,mz,pcgroup,raw_single_annotation,ion_grp,iso_grp,file="Temporary")
  ##load(file="Temporary")

  ## Group the ions into metabolites
  label_flg<-matrix(0L,nrow = length(pcgroup), ncol = 1)
  metabolite_grp<-matrix(0L,nrow = length(pcgroup), ncol = 1)
  metabolite_idx<-1

  ## The adducts are first grouped into metabolites
  if(is.matrix(ion_info_copy)) {
    for (k in 1:nrow(ion_info_copy)){
      ion_grp_idx<-which(ion_grp %in% ion_info_copy[k,1])
      if(length(ion_grp_idx)>0){
        label_flg[ion_grp_idx]<-1
        metabolite_grp[ion_grp_idx]<-metabolite_idx
        metabolite_idx<-metabolite_idx+1
      }
      else{
        stop("Error: ion group does not exist for ion group ", ion_info_copy[k, 1])
      }
    }
  } else {
    for (k in 1:nrow(ion_info)){
      ion_grp_idx<-which(ion_grp %in% ion_info[k,1])
      if(length(ion_grp_idx)>0){
        label_flg[ion_grp_idx]<-1
        metabolite_grp[ion_grp_idx]<-metabolite_idx
        metabolite_idx<-metabolite_idx+1
      }
      else{
        stop("Error: ion group does not exist for ion group ", ion_info[k, 1])
      }
    }}

  ## The isotopes are then grouped into metabolites. If a metabolite number
  ## already exists from the adducts of the isotopics ions, use the existing
  ## one. Otherwise, use a new one.

  for (k in 1:length(iso_grp_no)) {
    iso_grp_idx<-which(iso_grp %in% iso_grp_no[k])
    curr_metabolite_grp_no<-metabolite_grp[iso_grp_idx]
    if(any(curr_metabolite_grp_no>0)){
      if (sum(curr_metabolite_grp_no>0)==1){
        metabolite_grp[iso_grp_idx]<-curr_metabolite_grp_no[curr_metabolite_grp_no>0]
        label_flg[iso_grp_idx]<-1
      }
      else {
        stop("more than one metabolite group number for isotopes")
      }
    }
    else {
      metabolite_grp[iso_grp_idx]<-metabolite_idx
      label_flg[iso_grp_idx]<-1
      metabolite_idx<-metabolite_idx+1
    }
  }

  ## For the remaining ions, assign a metabolite number
  metabolite_grp[!label_flg] <- metabolite_idx:(metabolite_idx + sum(!label_flg)-1)
  uni_metabolite <- unique(metabolite_grp)
  mono_selector = matrix(0L,nrow = length(mono_flg), ncol = 1)
  for (i in 1:length(uni_metabolite)){
    idx <- which(metabolite_grp %in% uni_metabolite[i])
    idx_mono <- which(mono_flg[idx]>0)
    if (length(idx_mono)==1){
      mono_selector[idx[idx_mono]]<-1
    }
    else if(length(idx_mono)>=2){
      stop("more than one monoisotopic peak")
    }
    else{
      mono_selector[idx[1]]<-1;
    }
  }

  ## Change mono mz value to neutral monoisotopic mass
  mono_mass <- mono_mass - PROTO_MASS
  ## Output the parsed results into a csv file
  Peak.list<-raw_single_annotation
  attributes(Peak.list)
  Peak.list <- cbind(Peak.list,mono_mass,metabolite_grp,mono_flg,adduct_flg,iso_flg,ambiguity_flg,mono_selector)
  colnames(Peak.list)[(ncol(raw_single_annotation)+1):(ncol(raw_single_annotation)+7)]<-cbind("mono_mass","metabolite_group","monoisotopic_flg","adduct_flg","isotope_flg","ambiguity_flg","Selection_flg")
  copy_to(peak_db,Peak.list, name = "input_parsed", temporary = FALSE, overwrite = TRUE)
  rm(Peak.list, raw, raw_cleaned, raw_single_annotation)
  # write.table(Peak.list,file = paste(file.base,"input","parsed.csv",sep = "_"),sep = ",",row.names=FALSE)

}
