#' ################################################################################
#' ####################### Function in ICD_thrd tuning  #######################
#' ################################################################################
#'
#' #' Calculate the dlt drop off rate
#' #'
#' #' Calculate the dlt drop off rate
#' #'
#' #' @param list_simul The simulation data
#' #' @param chSize     The size of each cohort
#' #' @export
#'
#' dlt.rt.c1 <- function(list_simul, chSize){
#'
#'   sum(sapply(list_simul, function(a){
#'     sum(a$patlist$cycle[which(a$patlist$dlt == 1)] == 1)})) /
#'     (sum(sapply(list_simul, function(a){a$n.cohort})) * chSize)
#' }
#'
#' #' @describeIn Calculate the dlt rate for subsequent cycles
#' dlt.rt.subseq <- function(list_simul, chSize){
#'
#'   sum(sapply(list_simul, function(a){
#'     sum(a$patlist$cycle[which(a$patlist$dlt == 1)] > 1)})) /
#'     (sum(sapply(list_simul, function(a){a$n.cohort})) * chSize)
#' }
