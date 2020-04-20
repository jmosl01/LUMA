# We are currently receiving an error when running ParseCAMERA with gen.plots = TRUE
#Munge Data Modules
#if positive, run:
# ParseCAMERA(from.table = "Annotated", to.table = "output_parsed", CAMERA.obj = "anposGa")
#
# Error in .local(object, ...) : 
#   xcmsSource: file not found: E:\ORD Backup\Users\jmosley\Desktop\WLSSD 2014 Livers\mzML/Pooled QCs/Pooled_QC_Pos_10.mzML

#To get around this, set the variables in the script below and run *AFTER* running the InitWorkflow() module. Then execute the remaining script code to bypass the bug.

##Set these variables to be equivalent on your local system
data_dir <- "O:/Public/jmosl01/LUMA-Data"


##Execute this code only after setting the above variables for your local system
new_filepaths <- list.files(path = data_dir, recursive = TRUE, full.names = TRUE)

if(BLANK) new_filepaths <- new_filepaths[grepl("Blanks", new_filepaths)]
if(!BLANK) new_filepaths <- new_filepaths[!grepl("Blanks", new_filepaths)]

#If Positive mode, run this:
new_filepaths <- new_filepaths[grepl(ion.id[1],new_filepaths)]

anposGa@xcmsSet@filepaths <- new_filepaths

#If Negative mode, run this:
new_filepaths <- new_filepaths[grepl(ion.id[2],new_filepaths)]

annegGa@xcmsSet@filepaths <- new_filepaths
