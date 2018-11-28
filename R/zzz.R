# CRAN Note avoidance
if(getRversion() >= "2.15.1")
  utils::globalVariables(
    #Constructor function global variables
    c("peak_db",

    #general data utilities global variables
    "ion.mode","opt.dir","BLANK","mzdatafiles","CAMERA.file","XCMS.file","CAMERA.par",
    "EIC_ID","mz","rt", ".SD", "metabolite_group","anposGa","annegGa","xset4","xset",
    "mz1setpos","mz1setneg",

    #column names from in House Library modules
      "Name","Formula","Molecular.Weight","RT..Min.",

    #functions used within apply family of functions
    "sd", "as.numeric",

    #xcms and CAMERA variables
    "graph_method","CAMERA.ion.mode","rules","file.base",

    #database variables and functions
    "db.list","db.dir","lib.db","new.db","peak.db",

    #Exogenous_and_Solvent_Peak-Module global variables
    "ion.modes","pos_db","neg_db","blanks_pos_db","blanks_neg_db","lib_db","peak_db","mono_mass","meanRT",

    #other module global variables
    "Sexes","Classes","no.Samples","Endogenous","Sample.df","XCMS.par","search.par","DataFiles",
    "Corr.stat.pos","Corr.stat.neg","Endogenous.thresh","Peak.list.trimmed", "Solvent.ratio","Voidrt",
    "XCMS.par","cv.cutoff","gen.plots","keep.singletons","mf.cutoff","ppm.cutoff","rt.cutoff","ion.id",
    "Library.phenodata","Sample.phenodata","CT.ID","Plate.Number","Plate.Position",

    #Plotting function global variables
    "Correlation.stat")
  )
