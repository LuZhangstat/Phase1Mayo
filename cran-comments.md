## Test environments
* local x86_64-w64-mingw32 (64-bit), R 3.5.1
* ubuntu 12.04 (on travis-ci), R 3.5.1


## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking R code for possible problems ... NOTE
  eff_summary: no visible binding for global variable 'eff'
  eff_summary: no visible binding for global variable 'Var1'
  eff_summary: no visible binding for global variable 'Var2'
  eff_summary: no visible binding for global variable 'value'
  patlist.display: no visible binding for global variable '.'
  patlist.display: no visible binding for global variable 'PatID'
  phase1stage1 : model.file: no visible binding for global variable
    'beta_other'
  phase1stage1 : model.file: no visible binding for global variable
    'beta_dose'
  phase1stage1 : model.file: no visible binding for global variable 'N1'
  phase1stage1 : model.file: no visible global function definition for
    'inprod'
  phase1stage1 : model.file: no visible binding for global variable 'N2'
  phase1stage2 : model.file: no visible binding for global variable
    'beta_other'
  phase1stage2 : model.file: no visible binding for global variable
    'beta_dose'
  phase1stage2 : model.file: no visible binding for global variable 'N1'
  phase1stage2 : model.file: no visible global function definition for
    'inprod'
  phase1stage2 : model.file: no visible binding for global variable 'N2'
  phase1stage2 : model.file: no visible binding for global variable 'N3'
  phase1stage2 : model.file: no visible binding for global variable 'rho'
  plot.RunPRMD: no visible binding for global variable 'dose'
  plot.RunPRMD: no visible binding for global variable 'nTTP'
  plot.RunPRMD: no visible binding for global variable 'effcy'
  plot.SimPRMD: no visible binding for global variable 'dose'
  plot.SimPRMD: no visible binding for global variable 'q.2.5'
  plot.SimPRMD: no visible binding for global variable 'q.97.5'
  plot.SimPRMD : <anonymous>: no visible binding for global variable
    'nttp.c'
  plot.SimPRMD: no visible binding for global variable 'm.nttp'
  plot.SimPRMD: no visible binding for global variable 'group'
  plot.SimPRMD: no visible binding for global variable 'perc'
  print.summary.RunPRMD: no visible binding for global variable '.'
  print.summary.RunPRMD: no visible binding for global variable 'PatID'
  Undefined global functions or variables:
    . N1 N2 N3 PatID Var1 Var2 beta_dose beta_other dose eff effcy group
    inprod m.nttp nTTP nttp.c perc q.2.5 q.97.5 rho value

  All notes are caused by "no visible binding". 
  The "no visible binding" notes for function eff_summary, plot.SimPRMD are caused by options in ggplot function.
  The "no visible binding" notes for function patlist.display and print.summary.RunPRMD are caused by function "SummarySE" which calls ddply.
  The "no visible binding" notes for function phase1stage1 and phase1stage2 are caused because we need to write model in jags language,
  and pass option for jags model in R code. 
  
