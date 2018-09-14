
RunPRMD <-function(seed = 1234, patlist, patID_act = NULL,
                   cycle_act = NULL, dose_act = NULL, dlt_act= NULL,
                   doses = 1:6, cycles = 1:6,
                   tox.target = 0.28, p_tox1 = 0.2, p_tox2 = 0.2,
                   trialSize = 36, chSize = 3,
                   thrd1 = 0.28, thrd2 = 0.28, proxy.thrd = 0.1,
                   param.ctrl = list(), inits.list.set = list(),
                   n.chains = 1, effcy.flag = T, ICD.flag = T,
                   DLT.drop.flag = T, testedD = T,  IED.flag = T,
                   ICD_thrd = 0.3) {

  #' Simulation for an Adaptive, Multi-Stage Phase I Dose-Finding Design
  #'
  #' A function to implement simulation for an adaptive, multi-stage phase 1
  #' dose-finding design incorporating a longitudinal continuous efficacy
  #' outcome and toxicity data from multiple treatment cycles, proposed by Du et
  #' al(2017) and grade
  #'
  #' @param seed The seed of R's random number generator. Default is 1234
  #' @param patlist A list of the patient treatment records, which must contains
  #'   the following variables: \describe{ \item{PatID}{denotes the patient ID
  #'   where the elements are specified by cohort and subject number. For
  #'   example, "cohort2subject3" denotes the third subject in the second
  #'   cohort} \item{dose}{records the dose level assigned for each patient
  #'   through the whole treatment} \item{cycle}{shows the treatment cycle
  #'   information of each record} \item{nTTP}{records the corresponding nTTP
  #'   score.} \item{dlt}{indicates whether a DLT event is observed or not?}
  #'   \item{efficacy}{provides the continuous efficacy for each cycle. Required
  #'   when \code{effcy.flag == T}. The range of efficacy is (0, 1) and use -1
  #'   for missing efficacy response.} }
  #' @param patID_act A vector records the patients' ID who need dose
  #'   recommendation for next cycle. Default is \code{NULL}
  #' @param cycle_act A vector records the current cycle of patID_act. Default
  #'   is \code{NULL}
  #' @param dlt_act A vector indicates whether a dlt is observed in current
  #'   cycle for current patients. Default is \code{NULL}
  #' @param dose_act A vector records the current dose level of patID_act.
  #'   Default is \code{NULL}
  #' @param doses A vector of doses that users are going to explore. Default is
  #'   1:6, where dose 1 through dose 6 are being tested.
  #' @param cycles A vector of cycles that the treatment plans to go through.
  #'   Default is 1:6, where patients will experience up to 6 cycles of the
  #'   treatment
  #' @param tox.target  The target toxicity of the treatment. Default is 0.28.
  #'   See details below.
  #' @param p_tox1  The probability cutoff for cycle 1 toxicity. Default is 0.2.
  #'   See details below.
  #' @param p_tox2  The probability cutoff for later cycles toxicity beyond
  #'   cycle 1. Default is 0.2. See Details below.
  #' @param trialSize The maximum sample size for trial simulation. Default is
  #'   36. Must be the multiple of cohort size, represented by chSize
  #' @param chSize    The cohort size of patients recruited. Default is 3.
  #' @param thrd1 An upper bound of toxicity for cycle 1 of the treatment.
  #'   Default is 0.28. See Details below.
  #' @param thrd2 An upper bound of toxicity for late cycles of the treatment,
  #'   beyond cycle 1. Default is 0.28. See Details below
  #' @param proxy.thrd A distance parameter to define efficacious doses. Any
  #'   dose whose predicted efficacy is within proxy.thrd away from the largest
  #'   one among the safe doses will be declared an efficacious dose.
  #' @param param.ctrl A list specifying the prior distribution for the
  #'   parameters. p1_beta_intercept, the prior mean of intercept of toxicity
  #'   model assuming a normal prior; p2_beta_intercept, the precision (inverse
  #'   of variance) of intercept of toxicity model assuming a normal prior;
  #'   p1_beta_cycle, the prior mean of cycle effect of toxicity model assuming
  #'   a normal prior; p2_beta_cycle, the precision (inverse of variance) of
  #'   cycle effect of toxicity model assuming a normal prior; p1_beta_dose, the
  #'   prior minimum of dose effect of toxicity model assuming a uniform prior;
  #'   p2_beta_dose, the prior maximum of dose effect of toxicity model assuming
  #'   a uniform prior; p1_alpha, the prior mean vector of the parameters from
  #'   efficacy model assuming a multivariate normal prior; p2_alpha, the prior
  #'   precision matrix (inverse of covariance matrix) of the parameters from
  #'   efficacy model assuming a multivariate normal prior; p1_gamma0, the prior
  #'   mean of association parameter \eqn{\gamma} (See Du et al(2017)) of two
  #'   submodels of the joint model assuming a normal prior; p2_gamma0, the
  #'   prior precision (inverse of variance) of association parameter
  #'   \eqn{\gamma} of two submodels of the joint model assuming a normal prior.
  #'   Default is non-informative priors.
  #' @param inits.list.set Default is list()
  #' @param n.chains  No. of MCMC chains in Bayesian model fitting. Default is 1
  #' @param effcy.flag Whether efficacy data is considered in the model fitting
  #'   or not
  #' @param DLT.drop.flag Whether the patients should suspend the treatment when
  #'   observing DLT. Default is TRUE
  #' @param ICD.flag Whether we allow dose changing for cycle > 1 in stage 1
  #'   model or not? Default is TRUE. See details below
  #' @param testedD Default is TRUE. Whether we only allow ICD or IED among
  #'   cycle 1 tested dose level
  #' @param IED.flag Default is TRUE. Whether we allow dose changing for cycle >
  #'   1 in stage 2 model or not?
  #' @param ICD_thrd The cut-off point of the posterior toxicity probability in
  #'   defining ICD. Default is 0.3. See details below.
  #'
  #' @return \item{patlist}{The input data \code{patlist}} \item{doseA}{The
  #' recommended dose level for cycle 1 for new cohorts} \item{pat_rec}{The
  #' recommended dose for current patients for next cycle} \item{effcy.flag}{The
  #' input argument \code{effcy.flag}} \item{doses}{The input argument
  #' \code{doses}} \item{cycles}{The inpur argument \code{cycles}}
  #'
  #' @examples
  #'
  #' data("patlist_sim")
  #' patlist <- patlist_sim[-((length(patlist_sim$PatID) - 2):
  #'                        length(patlist_sim$PatID))]
  #' patID_act <- patlist_sim$PatID[(length(patlist_sim$PatID) - 2):
  #'                                length(patlist_sim$PatID)]
  #' cycle_act <- patlist_sim$cycle[(length(patlist_sim$PatID) - 2):
  #'                                 length(patlist_sim$PatID)] - 1
  #' dose_act <- rep(1, 3)
  #' dlt_act <- rep(0, 3)
  #'
  #' test <- RunPRMD(patlist = patlist, patID_act = patID_act,
  #'                 cycle_act = cycle_act, dose_act = dose_act,
  #'                 dlt_act = dlt_act, trialSize = 36, chSize = 3)
  #' summary(test)
  #' plot(test)
  #'
  #' test2 <- RunPRMD(patlist = patlist, patID_act = patID_act,
  #'                  cycle_act = cycle_act, dose_act = dose_act,
  #'                  dlt_act = dlt_act, trialSize = 75, chSize = 3)
  #' summary(test2)
  #' plot(test2, select_cycle = 1)
  #'
  #' \dontrun{
  #'
  #' }
  #'
  #' @import rjags
  #' @import coda
  #' @importFrom arrayhelpers vec2array
  #' @importFrom phase1RMD GenToxProb
  #' @importFrom coda gelman.diag
  #' @export


  #######################################################################
  # input argments checking (modified later)
  #######################################################################

  ## test whether trialSize/chSize is integer
  if(trialSize %% chSize != 0){
    stop("trialSize has to be the multiple of cohort size")
  }

  if(ICD.flag == T & IED.flag == F & effcy.flag == T){
    stop("No model with options ICD.flag == T & IED.flag == F & effcy.flag == T \n")
  }

  if( !is.null(patID_act) ){
    if(length(cycle_act)!=length(patID_act) |
       length(dose_act)!=length(patID_act) |
       length(dlt_act)!=length(patID_act))
    stop("cycle_act, dose_act and dlt_act should have the same length as patID_act when patID_act is specified")
  }

  #######################################################################
  # Model fitting report
  #######################################################################

  if(effcy.flag == F){
    if(ICD.flag == T){
      cat("Model: RMD with longitudinal toxicity only, ICD on \n\n")
    } else{
      cat("Model: RMD with longitudinal toxicity only \n\n")
    }
  } else {
    if(IED.flag == F){
      cat("Model: RMD with longitudinal toxicity and efficacy \n\n")
    } else {
      if(ICD.flag == T){
        cat("Model: RMD with longitudinal toxicity and efficacy, ICD on, IED on \n\n")
      } else {
        cat("Model: RMD with longitudinal toxicity and efficacy, ICD off at stage 1, IED (ICD) on at stage 2\n\n")
      }
      }
    }

  n.pats <- length(unique(patlist$PatID))
  Max.cohort <- trialSize / chSize
  if(n.pats %% chSize != 0){
    stop("The number of the patients is not a multiple of the cohort size\n\n")
  } else{
    n.cohort <- n.pats / chSize
  }

  if(is.null(patID_act)){
    act.index <- NULL
  }else{
    if(DLT.drop.flag == T){
      act.index <- which((cycle_act + 1) <= max(cycles) & dlt_act == 0)
    } else {
      act.index <- which((cycle_act + 1) <= max(cycles))
    }
  }

  if(effcy.flag == F){
    # For model without efficacy response. (effcy.flag == F)
    # ##ICD.flag == T, n.cohort < Max.cohort
    # ##1. Recommend dose for next cohort(, and current paitents for next cycle
    # ##(if(patID_act != NULL)))based on longitudinal toxicity
    # ##ICD.flag == F, n.cohort < Max.cohort
    # ##2. Recommend dose for next cohort based on longitudinal toxicity
    # ##ICD.flag == T, n.cohort = Max.cohort, the last cohort not observe dlt and not finish all cycles and not stop the treatment
    # ##3. No more new cohort, recommend dose for current paitents for next cycle based on longitudinal toxicity
    # ##ICD.flag == T, all finish
    # ##4. The trial is ending, and we are declaring MTD for each cycles
    # ##ICD.flag == F, n.cohort == Max.cohort
    # ##5. The trial is ending, and we are declaring MTD for cycle 1
    if(ICD.flag == F){
      if(n.cohort == Max.cohort){
        cat("The trial is ending, and we are declaring MTD for cycle 1\n\n")
      }else{
        cat("Recommend dose for next cohort based on longitudinal toxicity\n\n")
      }
    }else{
      if(n.cohort == Max.cohort){
        if(length(act.index) != 0){
          cat("No more new cohort, recommend dose for current paitents for next")
          cat(" cycle based on longitudinal toxicity\n\n")
        }else{
          cat("The trial is ending, and we are declaring MTD for each cycles\n\n")
        }
      }else{
        if(length(act.index) != 0){
          cat(" Recommend dose for next cohort, and current paitents for next")
          cat(" cycle based on longitudinal toxicity\n\n")
        } else {
          cat(" Recommend dose for next cohort based on longitudinal toxicity\n\n")
        }
      }
    }
  }else{
    # For model with efficacy responses (effcy.flag == T)
    # ##IED.flag == T, ICD.flag == T, n.cohort < Max.cohort/2
    # ##1. Recommend dose for next cohort, and current paitents for next cycle based on longitudinal toxicity
    # ##IED.flag == T, ICD.flag == F, n.cohort < Max.cohort/2
    # ##1. Recommend dose for next cohort based on longitudinal toxicity
    # ##IED.flag == T, n.cohort < Max.cohort & n.cohort >= Max.cohort/2
    # ##2. Recommend dose for next cohort, and current paitents for next cycle
    # ##based on longitudinal toxicity and efficacy
    # ##IED.flag == T, n.cohort = Max.cohort and not all patients finish the treatment
    # ##4. No more new cohort, recommend dose for current paitents for next cycle
    # ##based on longitudinal toxicity and efficacy response
    # ##IED.flag == T, all patients finished their treatment
    # ##6. The trial is ending, and we are declaring efficacious doses,
    # ##as well as recommending the lowest dose that is efficacious for all cycles
    # ##IED.flag == F, ICD.flag == F n.cohort < Max.cohort / 2
    # ##Recommend dose for next cohort based on longitudinal toxicity
    # ##IED.flag == F, n.cohort < Max.cohort & n.cohort >= Max.cohort/2
    # ##3. Recommend dose for next cohort based on longitudinal toxicity and efficacy ()
    # ##IED.flag == F, n.cohort == Max.cohort
    # ##5. The trial is ending, and we are declaring efficacious doses,
    # ##as well as recommending the lowest dose that is efficacious for cycle 1
    if(IED.flag == T){
      if(n.cohort < Max.cohort / 2){
        if(ICD.flag == T & length(act.index) != 0){
          cat("Recommend dose for next cohort, and current patients for next")
          cat(" cycle based on longitudinal toxicity\n\n")
        }else{
          cat("Recommend dose for next cohort based on longitudinal toxicity\n\n")
        }
      }else if(n.cohort < Max.cohort & n.cohort >= Max.cohort / 2){
        if(length(act.index) != 0){
          cat("Recommend dose for next cohort, and current patients for next ")
          cat("cycle based on longitudinal toxicity and efficacy\n\n")
        } else{
          cat("Recommend dose for next cohort based on longitudinal toxicity and efficacy\n\n")
        }
      } else {
        if(length(act.index) != 0){
          cat("No more new cohort, recommend dose for current paitents for next")
          cat(" cycle based on longitudinal toxicity and efficacy response\n\n")
        }else{
          cat("The trial is ending, and we are declaring efficacious doses")
          cat(", as well as recommending the lowest dose that is efficacious for all cycles\n\n")
        }
      }
    }else{
      if(n.cohort < Max.cohort / 2){
        cat("Recommend dose for next cohort based on longitudinal toxicity\n\n\n")
      }else if(n.cohort < Max.cohort & n.cohort >= Max.cohort / 2){
        cat("Recommend dose for next cohort based on longitudinal toxicity and efficacy\n\n\n")
      }else{
        cat("The trial is ending, and we are recommending the lowest dose that")
        cat(" is efficacious for cycle 1\n\n\n")
      }
    }
  }

  ########################################################################
  # precalculate data
  ########################################################################
  MaxCycle <- length(cycles)
  Numdose <- length(doses)

  # default prior settings for all model fitting #
  ctrl_param <- list(p1_beta_intercept = 0, p1_beta_cycle = 0,
                     p2_beta_intercept = 0.001, p2_beta_cycle = 0.001,
                     p1_beta_dose = 0, p2_beta_dose = 1000,
                     p1_alpha = c(0, 0, 0, 0), p2_alpha = diag(rep(0.001, 4)),
                     p1_rho = 0, p2_rho = 0.001)

  ## modifty setting if specific setting is given
  ctrl_param <- modifyList(ctrl_param, param.ctrl)

  #######################################
  # model fitting & dose recommendation #
  #######################################
  doseA <- min(doses)                    # recommend dose for cycle1
  pat_rec <- list(patID = NULL, cycle = NULL, dose = NULL)
  set.seed(seed)
  dose_flag <- 0            # indicating whether there is one dose level or not

  if (n.cohort == 1) {
    ## the first cohort, use 3 + 3 design
    dlt <- patlist$dlt[which(patlist$cycle == 1)]
    if (sum(dlt) == 0) {
      doseA <- doseA + 1          ## no dlt observed, escalate dose to 2
    } else if (sum(dlt <= 2)){
      doseA <- doseA              ## observe 1 or 2 dlts, same dose for cohort 2
    } else{
      cat("early stop in the first cohort, 3 dlt observed\n\n")
      doseA <- NA
      res <- list(doseA = doseA, pat_rec = pat_rec)
      attr(res,'class') <- 'RunPRMD'
      return(res)
    }

    pat_rec$patID <- patID_act[act.index]
    pat_rec$cycle <- cycle_act[act.index] + 1
    pat_rec$dose <- dose_act[act.index]

  } else if(n.cohort == 2 & all(patlist$dose[patlist$cycle == 1] == 1)) {
    ## the second cohort, use 3 + 3 design
    dlt <- patlist$dlt[which(patlist$cycle == 1)]
    ## 6 records of 2 cohort in the first cycle
    if (sum(dlt) == 1) {
      doseA <- doseA + 1          ## 1 dlt out of 6 records, escalate to dose 2
    } else if(sum(dlt) == 2){
      doseA <- doseA              ## 2 dlts out of 6 records, same dose
    } else {
      cat("more than two 2 dlts out of 6 records in the first 2 cohort, ")
      cat("\n stop the study\n\n")
      doseA <- NA
      res <- list(doseA = doseA, pat_rec = pat_rec)
      attr(res,'class') <- 'RunPRMD'
      return(res)
    }

    pat_rec$patID <- patID_act[act.index]
    pat_rec$cycle <- cycle_act[act.index] + 1
    pat_rec$dose <- dose_act[act.index]

  } else if(ifelse(effcy.flag == T,         # If efficacy is true
                   n.cohort < (Max.cohort / 2),
                   n.cohort <= Max.cohort)){
    if(n.cohort == Max.cohort & effcy.flag == F & length(act.index) == 0){
      cat("Stage III\n")
    }else{cat("Stage I\n")}
    ### stage 1 model fitting###
    retrieve_param <- c("beta_dose", "beta_other", "gamma")
    dose_flag <- ifelse(length(unique(patlist$dose)) == 1, 1, 0) # dose_flag = 1 represents one dose level only
    post_samples <- phase1stage1(patlist = patlist,
                                 ctrl_param = ctrl_param,
                                 n.iters = 5000, burn.in = 5000,
                                 retrieve_param = retrieve_param,
                                 dose_flag = dose_flag, n.chains = n.chains,
                                 inits.list.set = inits.list.set)
    if(n.chains > 1){
      ## use the 'potential scale reduction factor' to check convergence mixing, reference
      ## https://blog.stata.com/2016/05/26/gelman-rubin-convergence-diagnostic-using-multiple-chains/
      diag.converge <- c()
      for(k in 1:length(retrieve_param)){
        diag.converge <- c(diag.converge,
                           gelman.diag(as.mcmc.list(post_samples[[k]]),
                                       autoburnin = F)$psrf[ , 1])
      }
      if(any(diag.converge > 1.1)){cat(" MCMC fail to converge\n")}
    }

    ## dose recommendation ##
    if(dose_flag == 1){
      ## if only one dose is assigned in the study
      sim.betas <- as.matrix(rbind(post_samples$beta_other[, , 1],
                                   post_samples$beta_dose[, , 1]))
      mnTTP.dose1 <-
        mean(apply(sim.betas, 2,
                   function(o) { as.numeric(o[1] + o[2] < thrd1)}))

      ## mnTTP.dose1: Pr[(dose = 1) + (cycle = 1) < 0.28]

      if(mnTTP.dose1 < p_tox1) {
        cat("early stop, no recommended dose \n\n") # need to be modified
        doseA <- NA
        res <- list(doseA = doseA, pat_rec = pat_rec)
        attr(res,'class') <- 'RunPRMD'
        return(res)
      } else {
        doseA <- doseA + 1    ## increase the dose for the next cohort
        pat_rec$patID <- patID_act[act.index]
        pat_rec$cycle <- cycle_act[act.index] + 1
        pat_rec$dose <- dose_act[act.index] ## keep the same dose
      }
    } else {
      ## safe dose determination, for early termination##
      if(n.cohort > 2){
        # The termination only works for cohort >= 3 #
        allow.doses <- stg1.safe.dos(post_samples = post_samples,
                                     doses = doses, cycles = cycles,
                                     thrd1 = thrd1, thrd2 = thrd2,
                                     p_tox1 = p_tox1, p_tox2 = p_tox2,
                                     ICD.flag = ICD.flag)
        if (length(allow.doses) == 0) {
          cat("early stop, no allowable dose level \n\n\n")
          doseA <- NA
          res <- list(doseA = doseA, pat_rec = pat_rec)
          attr(res,'class') <- 'RunPRMD'
          return(res)
        }
      }
      Max_tested_doseA = max(patlist$dose[which(patlist$cycle == 1)]) ## calculate the maximum tested dose level for cycle1

      ## recommend dose for new cohort ##
      if(effcy.flag == F & n.cohort == Max.cohort & length(act.index)==0){
        doseA <- stg3.dos.rec(post_samples = post_samples, doses = doses,
                              tox.target = tox.target,
                              Max_tested_doseA = Max_tested_doseA)
      }else{
        doseA <- stg1.dos.rec(post_samples = post_samples, doses = doses,
                              tox.target = tox.target,
                              Max_tested_doseA = Max_tested_doseA)
      }

      ################################################
      # update the active patient for the next cycle #
      ################################################

      if(ICD.flag == T){
        ### allow dose modification ###
        ## if more than one dose is assigned in the study
        uniq_ID <- unique(patlist$PatID)

        ## safe dose for active patients for next cycle ##
        dos.rec.i.result <-
          stg1.dos.rec.i(post_samples = post_samples, uniq_ID = uniq_ID,
                         patID_act = patID_act, cycle_act = cycle_act,
                         rec_dose_act = dose_act,
                         Max_tested_doseA = Max_tested_doseA,
                         doses = doses, cycles = cycles, c1 = tox.target,
                         p1 = ICD_thrd, DLT.drop.flag = DLT.drop.flag,
                         y.dlt = dlt_act, testedD = testedD)

        pat_rec$patID <- dos.rec.i.result$patID_nxt
        pat_rec$cycle <- dos.rec.i.result$cycle_nxt
        pat_rec$dose <- dos.rec.i.result$rec_dose_nxt
      } else{
        ### no dose modification ###
        pat_rec$cycle <- cycle_act[act.index] + 1
        pat_rec$patID <- patID_act[act.index]
        pat_rec$dose  <- dose_act[act.index]
      }
    }
  }else if(n.cohort == Max.cohort & length(act.index)==0 & effcy.flag == T){
    #############################
    # model fitting for stage 3 #
    #############################
    cat("Stage III \n")
    retrieve_param = c("beta_dose", "beta_other", "alpha", "gamma")
    post_samples <- phase1stage2(patlist = patlist,
                                 ctrl_param = ctrl_param,
                                 n.iters = 5000, burn.in = 5000,
                                 retrieve_param = retrieve_param,
                                 n.chains = n.chains,
                                 # need at least two MCMC chains to check the convergence of MCMC chain
                                 dose_flag = dose_flag)

    ## Check MCMC convergency ##
    if(n.chains > 1){
      diag.converge <- c()
      for(k in 1:length(retrieve_param)){
        diag.converge <- c(diag.converge,
                           gelman.diag(as.mcmc.list(post_samples[[k]]),
                                       autoburnin = F)$psrf[ , 1])
      }
      if(any(diag.converge > 1.1)){cat("\n MCMC fail to converge")}
    }

    ## safe dose determination
    allow.doses <- stg1.safe.dos(post_samples = post_samples,
                                 doses = doses, cycles = cycles,
                                 thrd1 = thrd1, thrd2 = thrd2,
                                 p_tox1 = p_tox1, p_tox2 = p_tox2,
                                 ICD.flag = ICD.flag)

    if (length(allow.doses) == 0) {
      cat("\n early stop, no allowable dose level")
      doseA <- NA
      res <- list(doseA = doseA, pat_rec = pat_rec)
      attr(res,'class') <- 'RunPRMD'
      return(res)
    } else {
      ##################
      # dose recommend #
      ##################
      Max_tested_doseA = max(patlist$dose[which(patlist$cycle == 1)]) ## calculate the maximum tested dose level for cycle1

      doseA <- stg3.eff.rec(post_samples = post_samples,
                            allow.doses = allow.doses,
                            Max_tested_doseA = Max_tested_doseA,
                            proxy.thrd = proxy.thrd)
    }
  } else {
    cat("Stage II \n")
    ### stage 2 ###
    ## model fitting ##
    post_samples <- phase1stage2(patlist = patlist,
                                 ctrl_param = ctrl_param,
                                 n.iters = 5000, burn.in = 5000,
                                 retrieve_param =
                                   c("beta_dose", "beta_other", "alpha", "gamma"),
                                 n.chains = n.chains,
                                 dose_flag = dose_flag)
    ## safe dose determination ##
    allow.doses <- stg1.safe.dos(post_samples = post_samples,
                                 doses = doses, cycles = cycles,
                                 thrd1 = thrd1, thrd2 = thrd2,
                                 p_tox1 = p_tox1, p_tox2 = p_tox2,
                                 ICD.flag = ICD.flag)
    if (length(allow.doses) == 0) {
      cat("early stop, no allowable dose level\n\n")
      doseA <- NA
      res <- list(doseA = doseA, pat_rec = pat_rec)
      attr(res,'class') <- 'RunPRMD'
      return(res)
    }

    Max_tested_doseA = max(patlist$dose[which(patlist$cycle == 1)]) ## calculate the maximum tested dose level for cycle1

    ## dose recommendation for new cohort##
    doseA <- stg2.eff.rec(post_samples = post_samples,
                          allow.doses = allow.doses,
                          Max_tested_doseA = Max_tested_doseA,
                          proxy.thrd = proxy.thrd)

    ################################################
    # update the active patient for the next cycle #
    ################################################

    if(IED.flag == T){

      ### allow dose modification ###

      ## if more than one dose is assigned in the study
      uniq_ID <- unique(patlist$PatID)

      ## safe dose for active patients for next cycle ##
      dos.rec.i.result <-
        stg1.dos.rec.i(post_samples = post_samples, uniq_ID = uniq_ID,
                       patID_act = patID_act, cycle_act = cycle_act,
                       rec_dose_act = dose_act,
                       Max_tested_doseA = Max_tested_doseA,
                       doses = doses, cycles = cycles, c1 = tox.target,
                       p1 = ICD_thrd, DLT.drop.flag = DLT.drop.flag,
                       y.dlt = dlt_act, testedD = testedD)

      pat_rec$cycle <- dos.rec.i.result$cycle_nxt
      pat_rec$patID <- dos.rec.i.result$patID_nxt
      dos.rec.i.IED <- stg2.eff.rec.i(
        post_samples, cycle_nxt = dos.rec.i.result$cycle_nxt,
        rec_dose_nxt = dos.rec.i.result$rec_dose_nxt,
        proxy.thrd = proxy.thrd)
      pat_rec$dose <- dos.rec.i.IED$rec_dose_nxt
    } else {
      ### no dose modification ###
      pat_rec$cycle <- cycle_act[act.index]
      pat_rec$patID <- patID_act[act.index] + 1
      pat_rec$dose <- dose_act[act.index]
    }
  }
  cat("Recommend dose for new cohort: ", doseA, "\n")
  cat("\nFor patients: \n", pat_rec$patID, "\non cycle: \n", pat_rec$cycle,
      "\nWe suggest dose levels: \n", pat_rec$dose, "\n")
  res <- list(patlist = patlist, doseA = doseA, pat_rec = pat_rec,
              effcy.flag = effcy.flag, doses = doses,
              cycles = cycles)
  attr(res,'class') <- 'RunPRMD'
  return(res)
}



