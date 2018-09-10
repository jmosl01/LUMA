#' A Subset of liver metabolomics samples from the 2014 WLSSD 21 d effluent exposure study
#'
#' A dataset containing the normalized integrated intensity values for 444 endogenous metabolites in 48 liver samples.  Note that we define a metabolite as endogenous if it is present in at least half of our control samples.
#'
#' @format A data frame with 48 rows and 444 variables
#' \describe{
#' \item{Pos_629_Creatine.Creatine}{Normalized intensity values for creatine detected in positive ionization mode}
#' \item{Pos_846_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' \item{Pos_1951_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' \item{Pos_2044_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' \item{Pos_2317_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' \item{Pos_2375_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' \item{Pos_3160_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' \item{Neg_2118_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' \item{Neg_2461_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' \item{Neg_3431_Unidentified}{Normalized intensity values for a unique endogenous metabolite whose identity is unknown}
#' }
#' @source Unpublished data from the EPA metabolomics team in Athens, GA
"metabolites"

#' Factors for the samples (rows) in the values dataset
#'
#' A vector containing factor values to use for statistical analysis for the 48 liver samples in the values dataset
#'
#' @format A factor-formatted vector with 48 values
#' \describe{
#' \item{classes}{grouping factor corresponding to the sample class info for liver samples in the values dataset}
#' }
#' @source Unpublished data from the EPA metabolomics team in Athens, GA
"grouping"

