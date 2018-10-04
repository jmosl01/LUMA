# CRAN Note avoidance
if(getRversion() >= "2.15.1")
  utils::globalVariables(
    #Constructor function global variables
    c("peak_db",

    #general data utilities global variables
    "ion.mode","opt.dir","BLANK","mzdatafiles","CAMERA.file","XCMS.file",

    #other module global variables
    "Sexes","Classes","no.Samples","Endogenous","Sample.df","XCMS.par","search.par","DataFiles",
    "Corr.stat.pos","Corr.stat.neg","Endogenous.thresh","Peak.list.trimmed", "Solvent.ratio","Voidrt",
    "XCMS.par","cv.cutoff","gen.plots","keep.singletons","mf.cutoff","ppm.cutoff","rt.cutoff",


    #metadata column names from XCMS and CAMERA outputs
    "EIC_ID","mz","rt", ".SD", "metabolite_group",

    #column names from Annotated Library
      "RT..min.",

    #functions used within apply family of functions
    "sd", "as.numeric",

    #xcms and CAMERA variables
    "graph_method","CAMERA.ion.mode","rules","file.base",

    #database variables and functions
    "db.list","db.dir","lib.db","new.db","peak.db",

    #Possible Exogenous_and_Solvent_Peak-Module variables
    "ion.modes","pos_db","neg_db","blanks_pos_db","blanks_neg_db","lib_db","peak_db","mono_mass","meanRT",

    #Plotting function variables
    "Correlation.stat")
  )
