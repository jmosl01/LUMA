################################################################
##### Original code to create mSetObj for LUMA examples   ######
################################################################

library(LUMA)

#Fresh start
rm(list = ls())

#Set parameters for example mSetObj
data.type = "pktable"
anal.type = "stat"
paired = FALSE

##-----------------------------------------------------------------------------------------
## Initialize the metabolite object.
##-----------------------------------------------------------------------------------------
dataSet <- list();
dataSet$type <- data.type;
dataSet$design.type <- "regular";    # one factor to two factor
dataSet$cls.type <- "disc";
dataSet$format <- "rowu";
dataSet$paired <- paired;
analSet <- list();
analSet$type <- anal.type;

mSetObj <- list();
mSetObj$dataSet <- dataSet;
mSetObj$analSet <- analSet;
mSetObj$imgSet <- list();
mSetObj$msgSet <- list();                                 # store various message during data processing
mSetObj$msgSet$msg.vec <- vector(mode = "character");     # store error messages
mSetObj$cmdSet <- vector(mode = "character");             # store R command

# Set the class type for the mSetObj.
class(mSetObj) <- data.type

usethis::use_data(mSetObj, compress = "xz", overwrite = T)
