% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimPRMD.R
\name{SimPRMD}
\alias{SimPRMD}
\title{Simulation for a Multi-Stage Phase I Dose-Finding Design}
\usage{
SimPRMD(
  seed = 1234,
  numTrials = 100,
  doses = 1:6,
  cycles = 1:6,
  eff.structure = matrix(0, nrow = 6, ncol = 6),
  eff.Sigma = diag(6),
  eff.sd_trans = 1.5,
  tox.target = 0.28,
  p_tox1 = 0.2,
  p_tox2 = 0.2,
  trialSize = 36,
  chSize = 3,
  thrd1 = 0.28,
  thrd2 = 0.28,
  proxy.thrd = 0.1,
  tox.matrix = NULL,
  wm = matrix(c(0, 0.5, 0.75, 1, 1.5, 0, 0.5, 0.75, 1, 1.5, 0, 0, 0, 0.5, 1), byrow = T,
    ncol = 5),
  toxmax = 2.5,
  toxtype = NULL,
  intercept.alpha = NULL,
  coef.beta = NULL,
  cycle.gamma = NULL,
  param.ctrl = list(),
  n.iters = 10000,
  burn.in = 5000,
  thin = 2,
  n.chains = 1,
  effcy.flag = T,
  ICD.flag = T,
  DLT.drop.flag = T,
  testedD = T,
  IED.flag = T,
  ICD_thrd = 0.3
)
}
\arguments{
\item{seed}{The seed of R's random number generator. Default is 1234}

\item{numTrials}{An integer specifying the number of simulations}

\item{doses}{A vector of doses that users are going to explore. Default is
1:6, where dose 1 through dose 6 are being tested.}

\item{cycles}{A vector of cycles that the treatment plans to go through.
Default is 1:6, where patients will experience up to 6 cycles of the
treatment}

\item{eff.structure}{A matrix provides the mean of the multivariate
Gaussian distribution in efficacy data generation. Specifically, the
\eqn{(i, j)}th element represents the mean value of \eqn{i}th dose level
and \eqn{j}th cycle of the Gaussian distribution for efficacy data
generation. Default is a 6 by 6 zero matrix}

\item{eff.Sigma}{The covariance matrix of the multivariate Guassian
distribution in efficacy data generation. See details below.}

\item{eff.sd_trans}{A positive number controls the skewness of the
distribution of the efficacy response. Default is 1.5. See details below.}

\item{tox.target}{The target toxicity of the treatment. Default is 0.28.
See details below.}

\item{p_tox1}{The probability cutoff for cycle 1 toxicity. Default is 0.2.
See details below.}

\item{p_tox2}{The probability cutoff for later cycles toxicity beyond
cycle 1. Default is 0.2. See Details below.}

\item{trialSize}{The maximum sample size for trial simulation. Default is
36. Must be the multiple of cohort size, represented by chSize}

\item{chSize}{The cohort size of patients recruited. Default is 3.}

\item{thrd1}{An upper bound of toxicity for cycle 1 of the treatment.
Default is 0.28. See Details below.}

\item{thrd2}{An upper bound of toxicity for late cycles of the treatment,
beyond cycle 1. Default is 0.28. See Details below}

\item{proxy.thrd}{A distance parameter to define efficacious doses. Any
dose whose predicted efficacy is within proxy.thrd away from the largest
one among the safe doses will be declared an efficacious dose.}

\item{tox.matrix}{Optional. A four-dimension array specifying the
probabilities of the occurrences of certain grades for certain types of
toxicities, at each dose level and cycle under consideration. Dimension 1
refers to doses; dimension 2 corresponds to cycles of the treatment;
dimension 3 regards the types of toxicities while dimension 4 relates to
grades. If null, which is default choice, the arguments toxtype,
intercept.alpha, coef.beta, cycle.gamma must be provided to simulate this
array.}

\item{wm}{Clinical weight matrix, where toxicity types define the rows
while the toxicity grades define the columns. Usually solicited from
physicians.}

\item{toxmax}{The normalization constant used in computing nTTP score. For
details, see Ezzalfani et al(2013).}

\item{toxtype}{Only specified when tox.matrix is null. This argument, a
character vector, specifies toxicity types considered in the trial.}

\item{intercept.alpha}{Only specified when tox.matrix is null. A four
element numeric vector specifying the intercepts for the cumulative
probabilities of the occurrences of grades 0-4 of toxicities in
proportional odds model. See Details below.}

\item{coef.beta}{Only specified when tox.matrix is null. A n numeric vector
specifying the slope for dose in proportional odds model for n types of
toxicities. See Details below}

\item{cycle.gamma}{Only specified when tox.matrix is null. A scalar
controlling the cycle effect in simulation in proportional odds model.
See Details below}

\item{param.ctrl}{A list specifying the prior distribution for the
  parameters. \describe{\item{p1_beta_intercept}{the prior mean of
  intercept of toxicity model assuming a normal prior}

  \item{p2_beta_intercept}{the precision (inverse of variance) of intercept
  of toxicity model assuming a normal prior}

  \item{p1_beta_cycle}{the prior mean of cycle effect of toxicity model
  assuming a normal prior}

  \item{p2_beta_cycle}{the precision (inverse of variance) of cycle effect
  of toxicity model assuming a normal prior}

  \item{p1_beta_dose}{the prior minimum of dose effect of toxicity model
  assuming a uniform prior}

  \item{p2_beta_dose}{the prior maximum of dose effect of toxicity model
  assuming a uniform prior}

  \item{p1_alpha}{the prior mean vector of the parameters from efficacy
  model assuming a multivariate normal prior}

  \item{p2_alpha}{the prior precision matrix (inverse of covariance matrix)
  of the parameters from efficacy model assuming a multivariate normal
  prior}

  \item{p1_gamma0}{the prior mean of association parameter \eqn{\gamma}
  (See Du et al(2017)) of two submodels of the joint model assuming a
  normal prior}

  \item{p2_gamma0}{the prior precision (inverse of variance) of association
  parameter \eqn{\gamma} of two submodels of the joint model assuming a
  normal prior. } Default is non-informative priors. }}

\item{n.iters}{Total number of MCMC simulations. Default is 10,000.}

\item{burn.in}{Number of burn=ins in the MCMC simulation. Default is 5,000.}

\item{thin}{Thinning parameter. Default is 2.}

\item{n.chains}{No. of MCMC chains in Bayesian model fitting. Default is 1}

\item{effcy.flag}{Whether we include efficacy response in modeling or not?}

\item{ICD.flag}{Whether we allow dose changing for cycle > 1 in stage 1
model or not? Default is TRUE. See details below}

\item{DLT.drop.flag}{Whether the patients should suspend the treatment when
observing DLT. Default is TRUE}

\item{testedD}{Default is TRUE. Whether we only allow ICD or IED among
cycle 1 tested dose level}

\item{IED.flag}{Default is TRUE. Whether we allow dose changing for cycle >
1 in stage 2 model or not?}

\item{ICD_thrd}{The cut-off point of the posterior toxicity probability in
defining ICD. Default is 0.3. See details below.}
}
\value{
\item{senerio_sum}{contains \code{mnTTP.M} the matrix of mean nTTP
  for each dose and cycle and \code{pDLT.M} matrix of probability of
  observing DLT for each dose and cycle} \item{eff_sum}{When
  \code{effcy.flag == TRUE}, contains \code{eff.M} the mean efficacy for
  each dose and cycle and \code{err.cor.ls} A list with a length of dose
  levels numbers recording the marginal correlation matrix across cycles of
  efficacy data for each dose level} \item{list_simul}{A list of length
  numTrials. Each element includes \code{patlist} which records all the
  treatment and outcome information; \code{dose_aloca} which shows the
  cycle 1 dose allocation; \code{doseA} which saves the recommended dose
  level for cycle 1 at the end of the phase I simulation, equals "early
  break" if the trial was stop before finishing the trial; \code{n.cohort}
  indicates the last cohort in the trial; \code{pp.nTTPM} gives the
  posterior probability of nTTP less than target toxicity \code{tox.target}
  for all dose level any cycles and \code{message} saves the message of
  each trial.} \item{chSize}{The input argument \code{chSize}}
  \item{sim.time}{Time cost in simulation} \item{doses}{The input argument
  \code{doese}} \item{cycles}{The input argument \code{cycles}}
  \item{effcy.flag}{The input argument \code{effcy.flag}}
  \item{proxy.thrd}{The input argument \code{proxy.thrd}}
  \item{DLT.drop.flag}{The input argument \code{DLT.drop.flag}}
}
\description{
A function to implement simulations for a multi-stage phase 1 dose-finding
design incorporating a longitudinal continuous efficacy outcome and
toxicity data from multiple treatment cycles. The available models include
1-stage model with/without individualized dose modification, 3-stage model
with/without individualized dose modification, 3-stage model with
individualized dose modification on stage II and 3-stage model with
individualized dose modification on stage I and dose modification on stage
II.
}
\details{
The user can simulation efficacy response with different
  dose-efficacy and cycle-efficacy pattern using argument
  \code{eff.structure}, \code{eff.Sigma} and \code{eff.sd_trans}. The
  sampling process of efficacy response start from generating sample \eqn{z
  = {z1, \ldots, zd} } from multivariate Gaussian distribution \deqn{z ~
  MVN(\mu, V)}, where \eqn{\mu} and \eqn{V} are specified by
  \code{eff.structure} and \code{eff.Sigma}, respectively. Define
  \eqn{\phi} be the density of \eqn{N(0, \sigma^2)} with CDF \eqn{\Phi},
  and \eqn{\sigma^2} is set by \code{eff.sd_trans}. Then the efficacy
  response is calculated by taking the CDF of \eqn{z}: \deqn{x={x1, \ldots,
  xd} = \Phi(z) = { \Phi(z1), \ldots, \Phi(zd)}} is the generated efficacy
  response. Notice here the variance parameter \eqn{\sigma^2_{trans}}
  controls the variance of the generated efficacy.

The user can simulate longitudinal efficacy response with
  different dose-efficacy and cycle-efficacy pattern using argument
  \code{eff.structure}, \code{eff.Sigma} and \code{eff.sd_trans}. The
  sampling process of efficacy response starts from generating \eqn{z =
  {z1, \ldots, zd} } from multivariate Gaussian distribution \deqn{z ~
  MVN(\mu, V)}, where \eqn{\mu} and \eqn{V} are specified by
  \code{eff.structure} and \code{eff.Sigma}, respectively. Define
  \eqn{\phi} be the density of \eqn{N(0, \sigma^2)} with CDF \eqn{\Phi},
  where \eqn{\sigma^2} is set by \code{eff.sd_trans}. Then the efficacy
  measure is generated by taking the CDF of \eqn{z}: \deqn{x={x1, \ldots,
  xd} = \Phi(z) = { \Phi(z1), \ldots, \Phi(zd)}}. Notice here the variance
  parameter \eqn{\sigma^2_{trans}} controls the variance of the generated
  efficacy.

  \code{p_tox1}, \code{p_tox2}, \code{thrd1} and \code{thrd2} are used to
  define allowable (safe) doses the probability conditions for cycle 1:
  \deqn{P(nTTP1 < thrd1) > p_tox1} and for cycle > 1: \deqn{p(nTTP2 <
  thrd2) > p_tox2} , where \eqn{nTTP1} and \eqn{nTTP2} denote the posterior
  estimate of nTTP for cycle 1 and the average of cycle > 1. When we
  implement model with individualized dose modification, we only check the
  condition for cycle 1 for defining allowable (safe) doses.

  \code{ICD_thrd} are used to find ICD. ICD is defined as the maximum dose
  which satisfy the condition \deqn{P(nTTPi <  target.tox) > ICD_thrd} ,
  where \eqn{nTTPi} is the individualized posterior predicted nTTP score.
  The individualized dose modification for next cycle will not escalate
  more than 1 dose from the current dose.
}
\examples{
data("prob")      # load prob.RData from package phaseI, Details see "?prob"
data("eff")       # load eff.RData from package phaseI. Details see "?eff"

eff.structure = eff$Dose_Cycle_Meff[2, 2, , ]
eff.Sigma = eff$Sigma
eff.sd_trans = eff$sd_trans

wm <- matrix(c(0, 0.5, 0.75, 1, 1.5,
               0, 0.5, 0.75, 1, 1.5,
               0, 0, 0, 0.5, 1),
             byrow = TRUE, ncol
              = 5)                          # weighted matrix for toxicity matrix
                                            # nrow = No.of type; ncol = No. of grade
toxmax <- 2.5
tox.matrix <- prob["MTD4", "flat", , , , ]


#------- a flat dose-toxicity, dose-efficacy, cycle-efficacy pattern------#
\donttest{
simul1 <- SimPRMD(numTrials = 1, tox.matrix = tox.matrix,
                  eff.structure = eff.structure, eff.Sigma = eff.Sigma,
                  eff.sd_trans = eff.sd_trans, wm = wm, toxmax = toxmax,
                  trialSize = 36)
}
#------- a flat dose-toxicity pattern model ------#
\donttest{
simul2 <- SimPRMD(numTrials = 1, toxtype = c("H", "L", "M"),
                  intercept.alpha = c(1.9, 2.3, 2.6, 3.1),
                  coef.beta = c(-0.3, -0.2, -0.25),
                  cycle.gamma = 0, tox.target = 0.23,
                  thrd1 = 0.23, thrd2 = 0.23, p_tox1 = 0.2, p_tox2 = 0.2,
                  ICD.flag = FALSE, IED.flag = FALSE, effcy.flag = TRUE)

summary(simul2)
plot(simul2)
}

}
