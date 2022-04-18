# File:         model_init.R
# Author:       Kelly Chang
#               Zhihua Li
# Date:         Oct 2017
# Version:      1.0
# 
# Description:  Helper R function to initialize a drug binding model and run
#               voltage clamp simulations at different drug concentrations.
#

model_init<-function(modelname, states, pars, pnames, fulltimes, nbeats,outdir1,concid,events1){
#	print("ffff")
#	
	#print(paste0("stratcoin",concid))
	
#  print(pars)
	#events1=c()
	#print(source("Events.R"))
	#events1$time=3
	#print(events1)
	#print(events1$time[1])
	#events1$time[1]=13
	#print(source("Events.R"))
	

    force(modelname)
    force(nbeats)
    pidx<-match(names(pars), pnames, nomatch=0)
    sweeptimes<-c()
    sweepevents<-NULL
    #initstates=states
    mymodel<-list()
    timepoints=fulltimes
	
    initstates=states
    #print(states)
    mymodel$run_simulation<-function(initstates, pars, timepoints, events=NULL){
     # print(initstates)
        #pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
      if (concid<5) {
		  #print(paste0("concid",concid))
        try({out <- ode(initstates, timepoints,func=hergmod1, parms=pars,
            method="lsoda",atol = 1e-6, rtol = 1e-6)});
		}
		#write.csv(out,sprintf("outpre11_%s.csv",concid))
		
	#	print(paste0("dddconcidxxxxxx",concid))
		
	if (concid==5) {
		#print(paste0("dddconcidddddd",concid))
		timepointspre=seq(0,3600)
		try({out_pre <- ode(initstates, timepointspre,func=hergmod1, parms=pars,
							method="lsoda",atol = 1e-6, rtol = 1e-6)});
		timepoints=fulltimes
		#write.csv(out_pre,sprintf("outpre_%s.csv",concid))
		
		#timepoints=c(seq(-3600,-1),timepoints)
		print("preeeeeeeeeee")
		print(out_pre[3600,])
		initstates["RL"]=out_pre[3600,"RL"]
		initstates["R"]=out_pre[3600,"R"]
		initstates["RN"]=out_pre[3600,"RN"]
		initstates["L"]=out_pre[3600,"L"]
		initstates["N"]=10000000
		#print(pars)
		try({out <- ode(initstates, timepoints,func=hergmod1, parms=pars,
							method="lsoda",atol = 1e-6, rtol = 1e-6)});	
#		try({out <- ode(initstates, timepoints,func=hergmod1, parms=pars,
#			method="lsoda",atol = 1e-6, rtol = 1e-6,events=list(data=events1))});
			print("^^^^^^^^^^")
			print(out[3600,])
			
	}
	#write.csv(out,sprintf("outmain_%s.csv",concid))
	if (concid==5 & exists("out")) {
		print(out[3600,])
	}
        if(!exists("out")||inherits(out,"try-error")||length(out[,1])!=length(fulltimes) || !all(out[,1]==fulltimes) || any(is.nan(out)))
            return(c())
      #print("oooooooooooooooooooo")
	  out
    }
	
    sweeptimes<-fulltimes
    #sweepevents<-list(data=eventdata)
    #ctlsweeps<-mymodel$run_sweeps(ctlstates, ctlpars)
    # if(length(ctlsweeps)==0)
    #     stop("Solving error during sweep initialization!")

    #mymodel$controlsweeps<-function() ctlsweeps

    mymodel$run_drug<-function(ind,conc,states){
        states<-initstates
        pars[4:6]=ind
#print(ind)
#print(pars)
	
        out<-mymodel$run_simulation(states, pars, fulltimes)
        if(length(out)==0) {print ("root of shet here")}
           # return(out)
        #states[sidx!=0]<-out[nrow(out),sidx]
        #states<-out[nrow(out),sidx]
        
         # drugsweeps<<-mymodel$run_sweeps(states, pars,sweeptimes)
         drugsweeps=out
    }
#	print("1oooooooooooooooooooo")
	
    mymodel
}
