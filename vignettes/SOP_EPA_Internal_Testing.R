## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy.opts = list(width.cutoff=35),tidy = TRUE)

## ----echo=FALSE---------------------------------------------------------------
print("Error in Formula <<- Formula <- Annotated.Library[, \"Formula\"] :")
print("cannot change value of locked binding for 'Formula'")

## ----eval=FALSE---------------------------------------------------------------
#  detach(package:Formula, unload=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  input.dir <- system.file("extdata","Annotated_library.csv", package = "lcmsfishdata", mustWork = T)
#  x <- unlist(gregexpr("/",input.dir))
#  keep <- substr(input.dir, 1, x[length(x)] - 1)
#  keep

## ----warning=FALSE------------------------------------------------------------
library(LUMA)

## ----eval=FALSE---------------------------------------------------------------
#  InitWorkflow()

## ----eval=FALSE---------------------------------------------------------------
#  CullVoidVolume(from.table = "From CAMERA_with Minfrac", to.table = "Trimmed by RT", method = "mz")

## ----eval=FALSE---------------------------------------------------------------
#  AnnotatePeaklist(from.table = "Trimmed by RT", to.table = "Annotated")

## ----eval= FALSE--------------------------------------------------------------
#  ParseCAMERA(from.table = "Annotated", to.table = "output_parsed", CAMERA.obj = "anposGa")

## ----eval=FALSE---------------------------------------------------------------
#  ParseCAMERA(from.table = "Annotated", to.table = "output_parsed", CAMERA.obj = "annegGa")

## ----eval=FALSE---------------------------------------------------------------
#  CombineFeatures(from.table = "output_parsed", to.table = "Combined Isotopes and Adducts")

## ----eval=FALSE---------------------------------------------------------------
#  file.copy(from = paste0(keep,"/Peaklist_Neg_CorrPlots-GlutamicAcid_clear.pdf"), to = "./Peaklist_Neg_CorrPlots-GlutamicAcid_clear.pdf", overwrite = TRUE)
#  file.copy(from = paste0(keep,"/Peaklist_Pos_CorrPlots-Serine_muddy.pdf"), to = "./Peaklist_Pos_CorrPlots-Serine_muddy.pdf", overwrite = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  CullCV(from.table = "Combined Isotopes and Adducts", to.table = "Trimmed by CV")

## ----eval=FALSE---------------------------------------------------------------
#  CullMF(from.table = "Trimmed by CV", to.table = "Trimmed by MinFrac")

## ----eval=FALSE---------------------------------------------------------------
#  CullBackground(from.tables = c("Trimmed by MinFrac","Combined Isotopes and Adducts"),
#                 to.tables =   c("Peaklist_Pos_Solvent Peaks Removed", "Peaklist_Neg_Solvent Peaks Removed",    "Peaklist_Pos_Solvent Peaks Only", "Peaklist_Neg_Solvent Peaks Only"), method = "monoMass")

## ----eval=FALSE---------------------------------------------------------------
#  file.copy(from = paste0(keep,"/EIC_index_pos.txt"), to = "./EIC_index_pos.txt", overwrite = TRUE)
#  file.copy(from = paste0(keep,"/EIC_index_neg.txt"), to = "./EIC_index_neg.txt", overwrite = TRUE)
#  SimplyPeaklists(from.tables = c("Peaklist_Pos_Solvent Peaks Removed", "Peaklist_Neg_Solvent Peaks Removed"), to.table = "Peaklist_Combined", peak.db = "Peaklist_db")

## ----eval=FALSE---------------------------------------------------------------
#  NormalizePeaklists(from.table = "Peaklist_Combined",  to.table = "Peaklist_Normalized")

## ----eval=FALSE---------------------------------------------------------------
#  FormatForMetaboAnalystR(from.table = "Peaklist_Normalized", to.csv = "Peaklist_for_MetaboAnalyst", data.type = "pktable", anal.type = "stat")

## ----eval=FALSE---------------------------------------------------------------
#  FinalWorkflow(peak_db = peak_db, lib_db = lib_db)

## ----eval=FALSE---------------------------------------------------------------
#  sessionInfo()

## ----eval=FALSE---------------------------------------------------------------
#  
#  temp_db <- connect_peakdb(file.base = gen_filebase(mzdatafiles = DataFiles, BLANK, ion.id, IonMode))

## ----eval=FALSE---------------------------------------------------------------
#  mydf <- read_tbl(peak.tbls[1], peak.db = temp_db, asdf = T)
#  mydf

## ----eval=FALSE---------------------------------------------------------------
#  write_tbl(mydf = mydf, peak.db = temp_db, myname = "Reproduced_tbl")

## ----eval=FALSE---------------------------------------------------------------
#  mydf2 <- read_tbl("Reproduced_tbl", peak.db = temp_db, asdf = T)
#  identical(mydf,mydf2)
#  
#  ##Quick query of features in a table
#  LUMA:::get_features(mytbl = "Reproduced_tbl", peak.db = temp_db)

