#' @title CAMERA Parser in positive mode
#'
#' @export
#' @description Parses the CAMERA results using well-defined rules to eliminate conflicting annotations.
#' @param raw a data frame with variables as columns.  Should contain all output columns from XCMS and CAMERA, additional columns from IHL.search and a Minfrac column.  The last columns must be the CAMERA columns "isotopes","adduct","pcgroup" in that order.
#' @param rule a data frame containing the rule list used by CAMERA to annotate ion adducts and fragments.  Must contain the columns "name","nmol","charge","massdiff","oidscore","quasi","ips".
#' @param ion.mode a character string defining the ionization mode.  Must be "Positive"
#' @return data frame parsed version of the original data frame with additional columns "mono_mass","metabolite_group","monoisotopic_flg","adduct_flg","isotope_flg","ambiguity_flg","Selection_flg"
CAMERA.parser.pos=function(raw,rule,ion.mode){
  ##This code is a modified version of CAMERA_parser.m from the Ressom Omics Lab at Georgetown University (http://omics.georgetown.edu/) adapted for the R environment
  ##*******************************************************
  ##*******************POSITIVE   MODE*********************
  ##*******************************************************
  ##Note: the monoisotpic mass means the mass of [M+H]/[M-H] type of ion with the mono isotopes.
  if(ion.mode != "Positive") stop("Error: ion.mode must be \"Positive\".")
  raw.header <- colnames(raw)
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

  return(Peak.list)
}

#' @title CAMERA Parser in negative mode
#'
#' @export
#' @description Parses the CAMERA results using well-defined rules to eliminate conflicting annotations.
#' @param raw a data frame with variables as columns.  Should contain all output columns from XCMS and CAMERA, additional columns from IHL.search and a Minfrac column.  The last columns must be the CAMERA columns "isotopes","adduct","pcgroup" in that order.
#' @param rule a data frame containing the rule list used by CAMERA to annotate ion adducts and fragments.  Must contain the columns "name","nmol","charge","massdiff","oidscore","quasi","ips".
#' @param ion.mode a character string defining the ionization mode.  Must be "Negative"
#' @return data frame parsed version of the original data frame with additional columns "mono_mass","metabolite_group","monoisotopic_flg","adduct_flg","isotope_flg","ambiguity_flg","Selection_flg"
CAMERA.parser.neg=function(raw,rule,ion.mode){
  ##This code is a modified version of CAMERA_parser.m from the Ressom Omics Lab at Georgetown University (http://omics.georgetown.edu/) adapted for the R environment
  ##*******************************************************
  ##*******************NEGATIVE   MODE*********************
  ##*******************************************************
  ##Note: the monoisotpic mass means the mass of [M+H]/[M-H] type of ion with the mono isotopes.
  if(ion.mode != "Negative") stop("Error: ion.mode must be \"Negative\".")

  ## Input CAMERA output and Negitive Rule file
  raw.header <- colnames(raw)
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
      # % adduct annotation. Otherwise it will cause error. One Negsible
      # % solution in the future is to add new isotope group and split
      # % related isotope annotations.
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


  row.names(raw_single_annotation)<-1:nrow(raw_single_annotation)
  write.table(raw_single_annotation,file = "duplicate adduct removal for Negitive.csv",sep = ",",row.names=FALSE)
  #write.table(ambiguity_no,file="ambiguity_no.csv", sep = ",",row.names=FALSE)
  write.table(ambiguity_flg,file="ambiguity_flg.csv", sep = ",",row.names=FALSE)
  raw_single_annotation<-read.table(file = "duplicate adduct removal for Negitive.csv", sep=",", header = TRUE)
  #ambiguity_no<-read.table(file = "ambiguity_no.csv", sep=",", header = TRUE)
  ambiguity_flg<-read.table(file = "ambiguity_flg.csv", sep=",", header = TRUE)
  file.remove("duplicate adduct removal for Negitive.csv")
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
        mass_tmp<- (group_mz[j]*-z - massdiff)/nmol -PROTO_MASS
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
        if(curr_iso=='[M]-'){
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
        ## if an ion is multiple-charged [M]n-
        else if(length(grep(".+\\[\\M\\]\\d\\-",toString(group_iso[j])))>0){
          charge_state<-as.numeric(sapply(strsplit(toString(group_iso[j]),"\\M\\]|\\-"), "[", 2))
          mono_flg[ion_idx[j]]<-0
          if(adduct_flg[ion_idx[j]]!=1){
            # new mono_mass is assigned only when this ion has no adduct annotation,
            # otherwise the mono_mass is calculated based on adduct annotation to deal
            # with clusterion (such as [2M+2H]2+)
            mono_mass[ion_idx[j]]<-charge_state*group_mz[j]+(charge_state-1)*PROTO_MASS
          }
          iso_flg[ion_idx[j]]<-0
        }
        else {
          iso_flg[ion_idx[j]]<-1
          mono_flg[ion_idx[j]]<-0
        }
      } else{
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
      # % This seems to arise from CAMERA annotating two different ions with different
      # % m/z values as the same isotopic mass, presumably because they are close in mass.
      # % In this case, we don't use either monoisotopic m/z, but assume both are monoisotopic
      # % ions and recalculate them from the original m/z values.
      warning("More than one mono mass value for adducts. \nRemoving conflicting adduct annotation.")
      print(raw_single_annotation[grp_ion_idx,(ncol(raw_single_annotation)-2):ncol(raw_single_annotation)])
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
      warning("Error: more than one mono mass value for isotopes for isotope group ",iso_grp_no[k])
    }else {
      warning("Error: there is no monoisotopic ion found for isotope group ", iso_grp_no[k])
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
        stop("Error: more than one metabolite group number for isotopes")
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
      print(raw_single_annotation[idx,(ncol(raw_single_annotation)-2):ncol(raw_single_annotation)])
      stop("Error: more than one monoisotopic peak")
    }
    else{
      mono_selector[idx[1]]<-1;
    }
  }

  ## Change mono mz value to neutral monoisotopic mass
  mono_mass <- mono_mass + PROTO_MASS
  ## Output the parsed results into a csv file
  Peak.list<-raw_single_annotation
  attributes(Peak.list)
  Peak.list <- cbind(Peak.list,mono_mass,metabolite_grp,mono_flg,adduct_flg,iso_flg,ambiguity_flg,mono_selector)
  colnames(Peak.list)[(ncol(raw_single_annotation)+1):(ncol(raw_single_annotation)+7)]<-cbind("mono_mass","metabolite_group","monoisotopic_flg","adduct_flg","isotope_flg","ambiguity_flg","Selection_flg")

  return(Peak.list)
  }
