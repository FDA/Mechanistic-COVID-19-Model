

model_init<-function(modelname, states, pars, pnames, fulltimes, eventdata=NULL){
    force(modelname)
    
    pidx<-match(names(pars), pnames, nomatch=0)
    

    mymodel<-list()
    mymodel$states<- states
    mymodel$pars<- pars
    mymodel$fulltimes<- fulltimes
    mymodel$pidx<-pidx
    mymodel$run_simulation<-function(initstates, pars, fulltimes, events=NULL){
    
        pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
        try({out <- dede(initstates, fulltimes, "derivs", pars, dllname=modelname,
            initfunc="initmod", nout=3, rtol=1e-3, atol=1e-6, method="lsoda",
            events=events)});
        if(!exists("out")||inherits(out,"try-error")||length(out[,1])!=length(fulltimes) || !all(out[,1]==fulltimes) || any(is.nan(out)))
            return(c())
        colnames(out)[length(colnames(out))-1] <- "TcLung" #actually 0.15*this value is TcLung
        colnames(out)[length(colnames(out))] <- "APCLymphnode"
        out
    }

    mymodel$run_extrasimulation<-function(ind){  #ind has an extra variable of initial virus
    
       
          vidx<-names(states)=="FreeInfluenza"
          states[vidx] <- ind[length(ind)]
          pidx<- mymodel$pidx
          pars[pidx!=0]<-ind[pidx]
        pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
        try({out <- dede(states, fulltimes, "derivs", pars, dllname=modelname,
            initfunc="initmod", nout=3, rtol=1e-3, method="lsoda",
            events=NULL)});
        if(!exists("out")||inherits(out,"try-error")||length(out[,1])!=length(fulltimes) || !all(out[,1]==fulltimes) || any(is.nan(out)))
            return(c())
        colnames(out)[length(colnames(out))-1] <- "TcLung" #actually 0.15*this value is TcLung
        colnames(out)[length(colnames(out))] <- "APCLymphnode"
        out
    }
    
   

    
   

   

    mymodel
}
