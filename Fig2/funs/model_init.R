# File:         model_init.R
# Author:       Kelly Chang
#               Zhihua Li
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Helper R function to initialize a drug binding model and run
#               voltage clamp simulations at different drug concentrations.
#

model_init<-function(modelname, states, pars1, pnames, fulltimes,outdir1,fitpars){

  
    force(modelname)
    
    pidx<-match(names(pars1), pnames, nomatch=0)
	pars1[pidx!=0]<-fitpars[pidx]
    sweeptimes<-c()
    sweepevents<-NULL
    mymodel<-list()
    timepoints=fulltimes
    initstates=states
    pars1["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
	pars1["timeout"]<-5
	
	pidxs<-match(names(initstates), pnames, nomatch=0)
	initstates[pidxs!=0]<-fitpars[pidxs]

	modeldir<-"Data/"
	source(paste0(modeldir,"delaypars.R"))
	pidxs2=match(names(pars), names(pars1), nomatch=0)
	pars[pidxs2!=0]<-pars1[pidxs2]

	
	try({out <- dede(initstates, timepoints, "derivs", pars, dllname=modelname,initfunc="initmod",rtol=1e-3,atol=1e-6,nout=3,)});		
		
	if(!exists("out")||inherits(out,"try-error")||length(out[,1])!=length(timepoints) || !all(out[,1]==timepoints) || any(is.nan(out)))
		return(c())
	
	drugsweeps=out
	drugsweeps
}
  