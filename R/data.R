#' Toxicity probability matrix
#'
#' A 6 dimension array providing the toxicity probability.
#' The dimension of the toxicity probability matrix is 4 3 6 6 3 5 which
#' represents scenario, cycle effect, dose level, cycle number, tox type,
#' and tox grade. Scenarios are order in a way that first scenario is MTD = dose 2,
#' and second is MTD = dose 3, and third is MTD = dose 4, and fourth
#' is MTD = dose level 5. Cycle effect is ordered from decreasing toxicity
#' over cycles, flat, and increasing.
#'
#' @references ???

"prob"

#' Efficacy response generation parameters
#'
#' A list of 4 for the parameters in the efficacy response generation
#'
#' @return
#' \item{Dose_Cycle_Meff}{Dose-cycle mean efficacy matrix. A 4 dimension array
#' providing the mean of the multivariate Gaussian distribution in
#' efficacy data generation. The dimension of the Dose-cycle mean efficacy matrix
#' is 5 5 6 6 which represents dose efficacy pattern, cycle efficacy pattern,
#' dose and cycle. Patterns are ordered from increasing, flat, platform decreasing
#' and quadratic efficacy across dose levels and cycles.}
#'
#' \item{Sigma}{A 6 by 6 matrix, the covariance matrix of the multivariate
#' Guassian distribution in efficacy data generation.}
#'
#' \item{sd_trans}{A positive number controls the skewness of the distribution
#' of the efficacy response}
#'
#' \item{eff.M}{An array recording the mean of the generated efficacy data.
#' The dimension of the Dose-cycle mean efficacy matrix
#' is 5 5 6 6 which represents dose efficacy pattern, cycle efficacy pattern,
#' dose and cycle. Patterns are ordered from increasing, flat, platform decreasing
#' and quadratic efficacy across dose levels and cycles.}
#'
#' @references

"eff"


#' A list of patient information
#'
#' A list of the patient treatment records. The simulation with MTD = 4,
#' dose-toxicity trend, efficacy-toxicity trend and efficacy-cycle trend are all
#' flat.
#'
#' @return
#' \item{PatID}{denotes the patient ID where the elements are specified by
#' cohort and subject number. For example, "cohort2subject3" denotes the third
#' subject in the second cohort}
#' \item{dose}{records the dose level assigned for each patient
#' through the whole treatment}
#' \item{cycle}{shows the treatment cycle information of each record}
#' \item{nTTP}{records th corresponding nTTP score.}
#' \item{dlt}{indicates whether a DLT event is observed or not?}
#' \item{efficacy}{provides the continuous efficacy for each cycle. Required when
#' \code{effcy.flag == T}. The range of efficacy is (0, 1) and use -1 for
#' missing efficacy response.}
#' @references

"patlist_sim"
