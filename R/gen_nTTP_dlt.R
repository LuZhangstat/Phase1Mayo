#' gen_nTTP_dlt <- function(dose_cycle_PatID, tox.prob.M, wm, nTTP.all,
#'                          eff.structure, eff.Sigma, eff.sd_trans,
#'                          patlist){
#'
#'   #' Generate nTTP and DLT and efficacy
#'   #'
#'   #' \code{gen_nTTP_dlt} Generate nTTP and DLT and efficacy
#'   #'
#'   #' @param dose_cycle_PatID dose, cycle, PatID
#'   #' @param tox.prob.M       toxicity probability matrix with dimension: dose cycle type grade
#'   #' Tox.prob.M can be the output of the build-in matrix of function
#'   #' GenToxProb in package phase1RMD
#'   #' @param wm      (numerical matric) Toxicity weighted matrix, with row be
#'   #' the type of the toxicity and column be the toxicity grade
#'   #' @param nTTP.all The output of \code{\link{nTTP.array}}
#'   #' @param eff.structure Parameter setting for generating efficacy data (a 6 by 6 matrix)
#'   #' @param eff.Sigma     Correlation matrix for efficacy across cycles
#'   #' @param eff.sd_trans  The parameter controls the variance and skewness of the efficacy data
#'   #' @param patlist       the simulation data, for searching the history efficacy data
#'   #'
#'   #' @return
#'   #'
#'   #' @examples
#'   #'
#'   #'
#'   #' \dontrun{
#'   #'
#'   #' }
#'   #' @export
#'
#'   # sample the grade level for each toxiciy type
#'   n.tox.grade <- dim(tox.prob.M)[4]
#'   dose <- as.numeric(dose_cycle_PatID[1])
#'   cycle <- as.numeric(dose_cycle_PatID[2])
#'   nttp.indices <- apply(tox.prob.M[dose, cycle, , ], 1,
#'                         function(o){sample(1:n.tox.grade, 1, prob = o)})
#'
#'   # generate the corresponding dlt and the nTTP
#'   y.dlt <- as.numeric(max(wm[cbind(1:nrow(wm), nttp.indices)]) >= 1)
#'   y.nTTP <- nTTP.all[matrix(nttp.indices, nrow = 1)]
#'
#'   # generate efficacy data
#'   if(cycle == 1){
#'     # y.effz record the transfer variable
#'     y.effz <- rnorm(1, mean = eff.structure[dose, 1],  sd = eff.Sigma[1, 1])
#'     y.eff <- pnorm(y.effz, mean = 0, sd = eff.sd_trans)
#'   } else{
#'     ChoLD <- t(chol(eff.Sigma[1:(cycle - 1), 1:(cycle - 1)]))
#'     CinvC <- forwardsolve(ChoLD, eff.Sigma[cycle, 1:(cycle - 1)])
#'     y.effz.old <- patlist$effz[which(patlist$PatID == dose_cycle_PatID[3])]
#'     y.dose.old <- patlist$dose[which(patlist$PatID == dose_cycle_PatID[3])]
#'
#'     CinvM <- forwardsolve(ChoLD, y.effz.old -
#'                             eff.structure[cbind(y.dose.old, 1:(cycle - 1))])
#'     CondMu <- eff.structure[dose, cycle] + sum(CinvC * CinvM)
#'     CondV <- eff.Sigma[cycle, cycle] - sum(CinvC^2)
#'
#'     y.effz <- rnorm(1, mean = CondMu, sd = sqrt(CondV))
#'     y.eff <- pnorm(y.effz, mean = 0, sd = eff.sd_trans)
#'   }
#'
#'   #    if(cycle == 3){y.effz = runif(1, min = 0.2, max = 0.8); y.eff = abs(y.effz)} else{
#'   #    y.effz = -1; y.eff = -1
#'   #   }
#'
#'   return(c(y.nTTP = y.nTTP, y.dlt = y.dlt, y.eff = y.eff, y.effz = y.effz))
#' }
