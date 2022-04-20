if (caseMS=="severReady") {
	xunc=xuncS
	round_correct=0
	deadpop=number_p*0
	
	
	c1day=0
	c2day=13
	set.seed(seedn)
	r=rnorm(parnum_init*2, mean = 5.5, sd = 3)
	r1=r[r>1 & r<14]
	set.seed(seedn)
	incdayV=sample(r1,parnum_init)
	
	if (ip==0) {
		
		BCinf_parsA=BCinf_pars0$V2
		incDay=5.5
	}else{
		BCinf_parsA=0
		BCinf_parsA=as.numeric(allpar[ip,])
		incDay=incdayV[ip]

	}
	
	names(BCinf_parsA)=BCinf_pars0$V1
	BCinf_parsnotsort=BCinf_parsA[1:35]
	
	
	
	CSpar=BCinf_parsnotsort
	
	virusdiss=runif(parnum_init,.5*CSpar["FreeVirus"],2*CSpar["FreeVirus"])
	
	
	
	source(paste0(modeldir,"delaypars.R"))
	pidxs2=match(names(pars), names(BCinf_parsnotsort), nomatch=0)
	pars[pidxs2!=0]<-BCinf_parsnotsort[pidxs2]
	tki=576
}



if (caseMS=="mildReady") {
	xunc=xuncM
	round_correct=0
	deadpop=number_p*0
	c1day=1
	c2day=12
	set.seed(seedn)
	r=rnorm(parnum_init*2, mean = 7, sd = 3)
	r1=r[r>1 & r<14]
	set.seed(seedn)
	incdayV=sample(r1,parnum_init)
	
	if (ip==0) {
		BCinf_parsA=BCinf_pars0$V2
		incDay=7
	}else{
		BCinf_parsA=0
		BCinf_parsA=as.numeric(allpar[ip,])
		incDay=incdayV[ip]
		
	}
	
	names(BCinf_parsA)=BCinf_pars0$V1
	BCinf_parsnotsort=BCinf_parsA[1:35]
	CSpar=BCinf_parsA[36:44]
	pidxs2=match(names(BCinf_parsnotsort), names(CSpar), nomatch=0)
	BCinf_parsnotsort[pidxs2!=0]<-CSpar[pidxs2]
	
	virusdiss=runif(parnum_init,.5*CSpar["FreeVirus"],2*CSpar["FreeVirus"])
	
	
	source(paste0(modeldir,"delaypars.R"))
	pidxs2=match(names(pars), names(BCinf_parsnotsort), nomatch=0)
	pars[pidxs2!=0]<-BCinf_parsnotsort[pidxs2]
	tki=504
	
}

