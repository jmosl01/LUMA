# CRAN Note avoidance
if(getRversion() >= "2.15.1")
  utils::globalVariables(
    #metadata column names from XCMS and CAMERA outputs
    c("EIC_ID","mz","rt", ".SD", "metabolite_group",

    #column names from Annotated Library
      "RT..min.",

    #functions used within apply family of functions
    "sd", "as.numeric",

    #xcms and CAMERA variables
    "graph_method","CAMERA.ion.mode","rules","file.base",

    #database variables
    "db.list","db.dir","lib.db","new.db","peak.db",

    #Possible Exogenous_and_Solvent_Peak-Algorithm variables
    "ion.modes","pos_db","neg_db","blanks_pos_db","blanks_neg_db","lib_db","peak_db","mono_mass","meanRT",

    #Plotting function variables
    "Correlation.stat")
  )
