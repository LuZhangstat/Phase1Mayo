## summary and plot functions for RunPRMD ##

summary.RunPRMD <- function(object, ...){

  #' Summary an RunPRMD object
  #'
  #' Summary an RunPRMD object
  #'
  #' @param object RunPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #'
  #' @return
  #' \item{object}{The output of function \code{RunPRMD}}
  #' \item{mnttp.M}{The mean nTTP for all doses and cycles}
  #' \item{dlt.count.M}{The number of DLT for all doses and cycles}
  #' \item{eff.M}{The mean efficacy for all doses and cycles. Return \code{NULL}
  #' when \code{object$effcy.flag == TRUE}}
  #'
  #' @export
  #'
  #' @examples
  #' ## Check ?RunPRMD for example
  #'

  mnttp.M <- sapply(object$cycles, function(c){
    sapply(object$doses, function(d){
      mean(object$patlist$nTTP[which(object$patlist$cycle == c &
                                       object$patlist$dose == d)])
    })
  })
  colnames(mnttp.M) <- paste0("C", object$cycles)
  rownames(mnttp.M) <- paste0("D", object$doses)

  dlt.count.M <- sapply(object$cycles, function(c){
    sapply(object$doses, function(d){
      sum(object$patlist$dlt[which(object$patlist$cycle == c &
                                       object$patlist$dose == d)])
    })
  })
  colnames(dlt.count.M) <- paste0("C", object$cycles)
  rownames(dlt.count.M) <- paste0("D", object$doses)

  if(object$effcy.flag == T){
    eff.M <- sapply(object$cycles, function(c){
      sapply(object$doses, function(d){
        sum(object$patlist$efficacy[which(object$patlist$cycle == c &
                                       object$patlist$dose == d)])
      })
    })
  }else{eff.M = NULL}
  ans <- list(object = object, mnttp.M = mnttp.M,
              dlt.count.M = dlt.count.M, eff.M = eff.M)
  class(ans) <- "summary.RunPRMD"
  ans
}

print.summary.RunPRMD <- function(x, ...){

  #' Displays a useful description of a summary.RunPRMD object
  #'
  #' Displays a useful description of a summary.RunPRMD object
  #'
  #' @param x summary.RunPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #'
  #' @import RColorBrewer
  #' @importFrom dplyr mutate rowwise mutate_at starts_with select funs vars
  #' @import kableExtra
  #' @import knitr
  #' @importFrom utils capture.output
  #' @export

  if (!inherits(x, "summary.RunPRMD"))
    stop(gettextf("'x' must inherit from class %s",
                  dQuote("summary.RunPRMD")),
         domain = NA)

  n.cycle <- max(x$object$cycles); n.dose <- max(x$object$doses)
  color.pal <- brewer.pal(n = n.dose, name = "RdYlGn")[n.dose:1]
  l <- length(x$object$patlist$PatID)
  uniq_ID <- unique(x$object$patlist$PatID)

  n.patient <- length(uniq_ID)
  report <- matrix(NA, nrow = n.patient, ncol = n.cycle)
  colnames(report) <- c(paste0("cycle", 1:n.cycle))
  rownames(report) <- uniq_ID

  for(i in 1:l){
    report[which(uniq_ID == x$object$patlist$PatID[i]),
           x$object$patlist$cycle[i]] <-
      paste(format(round(x$object$patlist$nTTP[i],3),3), x$object$patlist$dlt[i],
            x$object$patlist$dose[i], sep = ",")
  }

  report.data <- data.frame(report)

  cat("\n The mean nTTP for all doses and cycles: \n")
  print(format(round(x$mnttp.M, 3), 3))
  cat("\n The number of DLT for all doses and cycles: \n")
  print(format(round(x$dlt.count.M, 3), 3))
  if(x$object$effcy.flag == T){
    cat("\n The mean efficacy for all doses and cycles: \n")
    print(format(round(x$eff.M, 3), 3))
  }
  cat("Recommend dose for cycle 1: ", x$object$doseA, "\n")
  cat("\nFor patients: \n", x$object$pat_rec$patID,
      "\non cycle: \n", x$object$pat_rec$cycle,
      "\nWe suggest dose levels: \n", x$object$pat_rec$dose, "\n")

  invisible(x)
  invisible(capture.output(print(
  report.data %>%
  mutate(PatID = row.names(.))%>%
  rowwise()%>%
  mutate_at(vars(starts_with("cycle")), funs(
    ifelse(!is.na(.),
           ifelse(as.numeric(unlist(strsplit(as.character(.), ","))[2]) == 0,
                  cell_spec(., color = "white", bold = T,
                            background = color.pal[as.numeric(strsplit(as.character(.), ",")[[1]])[3]]),
                  cell_spec(., color = "black", bold = T,
                            background = color.pal[as.numeric(strsplit(as.character(.), ",")[[1]])[3]])),
           cell_spec(., color = "white", background = "white"))
  ))%>%
  # rowwise()%>%
  #   mutate_at(vars(starts_with("cycle")), funs(format(round(as.numeric(unlist(
  #     strsplit(as.character(.), ","))[1]), 3), 3)))%>%
  select(PatID, starts_with("cycle"))%>%
  kable(escape = F, format = "html") %>%
  kable_styling() %>%
  footnote((general = "the entries in each cell are formatted in: nTTP, dlt, dose"))
  )))
}

plot.RunPRMD <- function(x, ..., select_cycle = x$cycles){

  #' nTTP and efficacy boxplots of a RunPRMD object
  #'
  #' nTTP and efficacy boxplots of a RunPRMD object
  #'
  #' @param x RunPRMD object to summarise
  #' @param ... other arguments ignored (for compatibility with generic)
  #' @param select_cycle A vector indication the cycle in the boxplot.
  #' Default is \code{cycle} of \code{x}.
  #'
  #' @examples
  #' ## Check ?RunPRMD for example
  #'
  #' @import ggplot2
  #'
  #' @export

  thislist <- list(...)
  Max.nTTP <- max(x$patlist$nTTP)
  select.index <- sapply(x$patlist$cycle, function(a){any(a == select_cycle)})
  nttp.dt <- data.frame(
    nTTP = as.vector(x$patlist$nTTP[select.index]),
    cycle = factor(x$patlist$cycle[select.index], levels = select_cycle),
    dose = factor(x$patlist$dose[select.index], levels = x$doses))
  pnTTP <- ggplot(data = nttp.dt, aes(x = dose, y = nTTP, fill = cycle)) +
    stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() +
    facet_wrap(~cycle) + xlab("Doses") + ylab("nTTP")
    scale_x_discrete(breaks = x$doses, labels = x$doses)
  print(pnTTP)
  if(x$effcy.flag == T){
    readline(prompt="Press [enter] to check the efficacy box-plot")
    effcy.dt <- data.frame(
      effcy = as.vector(x$patlist$efficacy[select.index]),
      cycle = factor(x$patlist$cycle[select.index], levels = select_cycle),
      dose = factor(x$patlist$dose[select.index], levels = x$doses))
    peff <- ggplot(data = effcy.dt, aes(x = dose, y = effcy, fill = cycle)) +
      stat_boxplot(geom ='errorbar', width = 0.5) + geom_boxplot() +
      facet_wrap(~cycle) + xlab("Doses") + ylab("Efficacy")
      scale_x_discrete(breaks = x$doses, labels = x$doses)
    print(peff)
  }
}

