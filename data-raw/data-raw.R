### grouping Variables ---------------------------------

mydata <- read.csv(file = "inst/extdata/Peaklist.csv", sep = "," , header = TRUE, stringsAsFactors = FALSE)
if("Sex" %in% names(mydata) && "Exposure.Class" %in% names(mydata)){
  break
} else {
  Sexes <- c("M_","F_")
  Exposure.classes <- c("100_Effluent","20_Effluent","5_Effluent","Control")
  sample.names <- mydata[2:sum(match("EIC_ID", mydata[,1]),-1),1]
  class.list <- lapply(Exposure.classes, function(x) grep(x,sample.names))
  names(class.list) <- Exposure.classes
  sex.list <- lapply(Sexes, function(x) grep(x,sample.names))
  names(sex.list) <- Sexes
  col.classes <- merge(stack(sex.list),stack(class.list), by = "values")
  names(col.classes) <- c("Sample.No","Sex","Exposure.Class")
  col.classes$Sex <- gsub("_","",col.classes$Sex)
  rownames(mydata)
  rownames(col.classes) <- sapply(as.numeric(rownames(col.classes)), function(x) x+1)
  new.mydata <- merge(col.classes,mydata, all.y = TRUE, by = "row.names")
  new.mydata$Row.names <- as.numeric(new.mydata$Row.names)
  mydata <- new.mydata[do.call(order,new.mydata),]
  row.names(mydata) <- NULL
  mydata <- mydata[,-c(1:2)]
}

## The next bit of code separates the dataset by Sex and Exposure Class
male.mydata <- subset(mydata, mydata$Sex == "M")
# female.mydata <- subset(mydata, mydata$Sex == "F")
# sets <- c("male.mydata","female.mydata")
#Males
grouping <- male.mydata$Exposure.Class
metabolite.names <- names(male.mydata[4:length(male.mydata)])
metabolites <-male.mydata[, 4:length(male.mydata) ]
str(metabolites)
row.names(metabolites) <- NULL
attributes(metabolites)
myANOVA <- run_stats.ANOVA(metabolites,grouping)
reduced.metabolites <- myANOVA[which(myANOVA[,"p.value"] <= 0.001),]
rownames(reduced.metabolites)
metabolites <- metabolites[,rownames(reduced.metabolites)]
write.csv(metabolites, file = "data-raw/metabolites.csv")
devtools::use_data(metabolites, overwrite = TRUE)
write.csv(grouping, file = "data-raw/grouping.csv")
devtools::use_data(grouping, overwrite = TRUE)
