## summary functions ##

print.sim_phase1 <- function(object, ...){
  thislist <- list(...)
  cat("\n The therotical mean nTTP matrix is:\n")
  print(object$senerio_sum$mnTTP.M)
  cat("\n The therotical probability of observing DLT is: \n")
  print(object$senerio_sum$pDLT.M)
  cat("\n There are in total", length(object$list_simul), "simulations \n ")
  cat("\n The cycle 1 dose allocation table of this simulation is: \n")
  print(paste0(format(round(object$alloc.perc*100, 2), 2), "%"))
}

summary.sim_phase1 <- function(object, ...){
  #' Displays a useful description of a sim_phase1 object
  #'
  #' Displays a useful description of a sim_phase1 object
  #'
  #' @param object sim_phase1 object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #'
  #' @examples
  #' ## Check ?sim_phase1 for example

  thislist <- list(...)
  dose.sbsq <- unlist(sapply(object$list_simul, function(a){
    a$patlist$dose[which(a$patlist$cycle!= 1)]}))
  sbsq.alloc <- table(factor(dose.sbsq[which(dose.sbsq != "early break")],
                             levels = 1:dim(object$senerio_sum$mnTTP.M)[1])) /
    length(dose.sbsq)
  cat("\n The therotical mean nTTP matrix is:\n")
  print(object$senerio_sum$mnTTP.M)
  cat("\n The therotical probability of observing DLT is: \n")
  print(object$senerio_sum$pDLT.M)
  cat("\n There are in total", length(object$list_simul), "simulations \n ")
  cat("\n The dose allocation table for cycle 1 of this simulation is: \n")
  print(paste0(format(round(object$alloc.perc*100, 2), 2), "%"))
  cat("\n The dose allocation table for cycle > 1 of this simulation is: \n")
  print(paste0(format(round(sbsq.alloc*100, 2), 2), "%"))

}
