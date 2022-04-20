options(warn = -1)
source("random.R")
library(grid)
numsiz=16
labsiz=18
BCinf_pars0=read.csv(paste0("Inputs_parameters/","viruscytokines_pars.txt"),sep = "",header = F,stringsAsFactors =FALSE)
out_plot_paperall=Bioout_plot_paperall=c()
All_adm_t=c()
DamP=15
allpar=par_rand1
if (i1==0) {totlaip=seq(0,parnum_init)} else {totlaip=parnum}
ipun=-1
cross5per1=c()
Lip=1
usedip=c()
for (ip in totlaip) {
	ipun=ipun+1
	source("parameters.R")	
	
	BCinf_pars=pars
	BCinf_pars["timeout"]=30
	BCinf_pars["starttime"]=0
	BCinf_pars["unnamed"]=1
	initstates=states
	if (ip==0) {initstates["FreeVirus"]=CSpar["FreeVirus"]}else{initstates["FreeVirus"]=virusdiss[ip]}
	
	
	VIRUS=read.csv(paste0("Data_Covid/Viral_Titer/",caseMS,".csv"))
	if (caseMS=="mildReady") {
		VIRUS=VIRUS[-1:-2,]
		VIRUS=VIRUS[VIRUS$Rx<23,]
		addday=7
		adddaysim=incDay
		
	}
	
	if (caseMS=="severReady") {
		VIRUS=VIRUS[VIRUS$Rx<23,]
		addday=5.5
		adddaysim=incDay
		
	}
	
	VIRUS$time=(VIRUS$Rx*24+addday*24)
	IL6=read.csv(paste0("Data_Covid/IL6/",caseMS,".csv")); IL6$y1mean=IL6$y1mean*75/90 
	IL6=read.csv(paste0("Data_Covid/IL6/",caseMS,".csv"))
	IL6$y1mean=IL6$y1mean*75/90 
	IL6$y1min=IL6$y1min*75/90 
	IL6$y1max=IL6$y1max*75/90 
	IL6$y1mean=IL6$y1mean/75
	IL6$y1min=IL6$y1min/75
	IL6$y1max=IL6$y1max/75
	IL6$time=(IL6$day+addday)*24
	TNAIVE=read.csv(paste0("Data_Covid/lymphocyte/",caseMS,".csv"))
	TNAIVE$time=(TNAIVE$day+addday)*24
	TNAIVE_percent=(2*10^9-10^9*TNAIVE["ymean"])/(2*10^9)*100
	TNAIVE_percent1plot=cbind(TNAIVE$time,TNAIVE_percent)
	ErrbarTn=(2*10^9-10^9*TNAIVE[c("ymin","ymax")])/(2*10^9)*100
	ErrbarTn$time=TNAIVE$time
	names(TNAIVE_percent1plot)=c("time","Tnaivepercent")
	
	
	isWindows="FALSE"	
	plotlist = list()
	mymodel<-list()
	source("PK.R")
	out_plot_paper=c()
	usedip[Lip]=ip
	ipcoef=runif(1000*(ipun+1),.98,1.02)
	ipcoef1=sample(ipcoef,length(pars))
	if (ip%in%usedip & i1>0) {pars=ipcoef1*pars;}
	Lip=Lip+1
	pars["timeout"]=30
	pars["unnamed"]=1
	pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
	
	if (ModelPlatform=="R") {
	source("Model/delaymodelfun.R")	
	try({out_plot_paper=dede(y=initstates,times=timepoints,func=modelfun,method="lsoda",parms=pars,atol = 1e-6,rtol=1e-3,events=list(data=eventdata))})
} else{	
	isWindows<-Sys.info()[["sysname"]]=="Windows"
	
	if (isWindows) {print("for running in Windows ModelPlatform should be=  'R'");WindowsError}
		extension<-ifelse(isWindows, ".dll", ".so")
		dyn.load(paste0(modeldir,modelname,extension))
	#}
try({out_plot_paper <- dede(initstates, timepoints, "derivs", pars, dllname=modelname,initfunc="initmod",rtol=1e-3,atol=1e-6,nout=3,events=list(data=eventdata))});		
}

	out_plot_paper=data.frame(out_plot_paper)
	out_plot_paper$ip=ip	
	out_plot_paper$ipun=ipun
	out_plot_paper$initVirus=(out_plot_paper$FreeVirus)[out_plot_paper$time==20]
	out_plot_paper$maxIL6=max(out_plot_paper$IL6)
	out_plot_paper$damage=max(4E8-out_plot_paper$ResistantEpithelium-out_plot_paper$EpitheliumLung)/4E8*100
	out_plot_paper$damage0=(4E8-out_plot_paper$ResistantEpithelium-out_plot_paper$EpitheliumLung)/4E8*100
	out_plot_paper$maxFreeVirus=max(out_plot_paper$FreeVirus)
	Il6KP=out_plot_paper$EpiKilledbyIL6/(out_plot_paper$EpiKilledbyIFN+out_plot_paper$EpiKilledbyTc+out_plot_paper$EpiKilledbyVirus+out_plot_paper$EpiKilledbyIL6)
	VKp=out_plot_paper$EpiKilledbyVirus/(out_plot_paper$EpiKilledbyIFN+out_plot_paper$EpiKilledbyTc+out_plot_paper$EpiKilledbyVirus+out_plot_paper$EpiKilledbyIL6)
	out_plot_paper$il6killedpercent=max(Il6KP[Il6KP!="Inf" & Il6KP!="NaN"])
	out_plot_paper$Vkilledpercent=max(VKp[VKp!="Inf" & VKp!="NaN"])
	daybase=min(eventtimes)
	reme1=daybase%%2
	if (reme1==1){daybase=daybase+1}else{daybase=daybase}
	out_plot_paper$eventtimes=daybase
	if (i1==0) {out_plot1_opt0Biolog=out_plot_paper}
	if (i1==0) {out_plot1_opt0Biolog_opt_alltime=data.frame(out_plot_paper)}
	if (i1==0) {out_plot1_opt0_opt_alltime=data.frame(out_plot_paper)}
	if (ip ==0 & i1==1) {out_plot1_opt0Biolog=out_plot_paper}
	if (ip ==0 & i1==1) {out_plot1_opt0Biolog_opt_alltime=data.frame(out_plot_paper)}
	if (ip ==0 & i1==2) {out_plot1_opt0_opt_alltime=data.frame(out_plot_paper)}
	
	Bioout_plot_paper=out_plot_paper
	out_plot_paper[,"time"]=out_plot_paper[,"time"]-daybase 
	out_plot_paper=out_plot_paper[out_plot_paper[,"time"]>=0,]
	if (ip ==0 ) {out_plot1_opt0=out_plot_paper}
	out_plot_paper$minVir=min(out_plot_paper$FreeVirus)
	if (i1!=0) {
		allt1=c(ip,min(eventtimes))
		All_adm_t=rbind(All_adm_t,allt1)
		Bioout_plot_paperall=rbind(Bioout_plot_paperall,Bioout_plot_paper)
		Bioout_plot_paperall=data.frame(Bioout_plot_paperall)
		out_plot_paperall=rbind(out_plot_paperall,out_plot_paper)
	}
	
	if (ip!=0) {
		out_plot1_paper_ip=out_plot_paper
		out_plot1_paper_ip$Damageper=(4E8-out_plot1_paper_ip$ResistantEpithelium-out_plot1_paper_ip$EpitheliumLung)/4E8*100
		recovery_threshold=5
		out_plot1_paper_ip_m1=out_plot1_paper_ip$Damageper-recovery_threshold
		updn <- c(diff(sign(out_plot1_paper_ip_m1)))
		cross5per <- which(updn == -2)+1	
		cross5pert0=out_plot1_paper_ip$time[cross5per[1]]
		
		minDamageBeforInj=min(out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time==0]) 
		if (minDamageBeforInj<5) {cross5pert0=-10000}
		if (is.na(cross5pert0)) {cross5pert0=10000}
		
		damageD0=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*0][1]
		damageD1=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*1][1]
		damageD10=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*10][1]
		damageD11=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*11][1]
		damageD12=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*12][1]
		damageD15=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*15][1]
		damageD28=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*28][1]
		
		
		cross5per1=rbind(cross5per1,c(ip,ipun,cross5pert0,min(eventtimes),damageD0,damageD1,damageD10,damageD11,damageD12,damageD15,damageD28))
		
	}
	
} 
colnames(cross5per1)=c("ipr","ip","Recoveredt","minevent","damageD0","damageD1","damageD10","damageD11","damageD12","damageD15","damageD28")
cross5per1=data.frame(cross5per1)


if (i1!=0) {
	colnames(All_adm_t)=c("patient","hr")
	out_plot_paperall=out_plot_paperall[out_plot_paperall$ip!=0,]
	out_plot1_paper=data.frame(out_plot_paperall)
	out_plot2_paper=data.frame(out_plot_paperall)
	
	
	out_plot1=data.frame(out_plot1_opt0Biolog)
	out_plot2=data.frame(out_plot1_opt0)
	IL6q=IgMq=IgGq=Tnaiveq=IFNA1q=FreeVirusq=REcell=Damageper1=c()
	out_plot1_paper$REcell=out_plot1_paper$ResistantEpithelium+out_plot1_paper$EpitheliumLung
	for (ti1 in seq(0,(max(out_plot1_paper[,"time"])-24),2)) {
		tout_plot1=out_plot1_paper[out_plot1_paper[,1]==ti1,]
		tq=(lapply(tout_plot1,quantile,probs=c(0.05,.95),na.rm=TRUE))
		
		REcell=rbind(REcell,c(ti1,tq$REcell))
		IL6q=rbind(IL6q,c(ti1,tq$IL6))
		FreeVirusq=rbind(FreeVirusq,c(ti1,tq$FreeVirus))
		IFNA1q=rbind(IFNA1q,c(ti1,tq$IFNA1))
		Tnaiveq=rbind(Tnaiveq,c(ti1,tq$Tnaive))
		IgGq=rbind(IgGq,c(ti1,tq$IgG))
		IgMq=rbind(IgMq,c(ti1,tq$IgM))
		
	}
	
	BIL6q=BIgMq=BIgGq=BTnaiveq=BIFNA1q=BFreeVirusq=BREcell=BEpitheliumLungq=c()
	Bioout_plot_paperall$FreeVirus=log10(Bioout_plot_paperall$FreeVirus)
	
	for (ti1 in seq(0,max(out_plot1[,"time"]),2)) { 
		
		tout_plot1=Bioout_plot_paperall[Bioout_plot_paperall[,1]==ti1,]
		tq=(lapply(tout_plot1,quantile,probs=c(0.1,.9),na.rm=TRUE))
		
		BIL6q=rbind(BIL6q,c(ti1,tq$IL6))
		BFreeVirusq=rbind(BFreeVirusq,c(ti1,tq$FreeVirus))
		BIFNA1q=rbind(BIFNA1q,c(ti1,tq$IFNA1))
		BTnaiveq=rbind(BTnaiveq,c(ti1,tq$Tnaive))
		BIgGq=rbind(BIgGq,c(ti1,tq$IgG))
		BIgMq=rbind(BIgMq,c(ti1,tq$IgM))
		BEpitheliumLungq=rbind(BEpitheliumLungq,c(ti1,tq$EpitheliumLung))
		
		
	}
	
}


Damageper1=c()

cross5per_step2=cross5per1
if (i1!=2) { 
	cross5per1=cross5per1[cross5per1$Recoveredt>0,]
}
maxtime=29*24


#if (i1==1) {write.csv(cross5per1,sprintf("%s_i1is1_cross5per1_%s.csv",outdir2,caseMS))}
#if (i1==2) {write.csv(cross5per1,sprintf("%s_i1is2_cross5per1_%s.csv",outdir2,caseMS))}
if (i1==0) {
	if (caseMS=="mildReady") { 
#		write.csv(cross5per_step2,sprintf("%s_cross5per1_%s.csv",outdir2,caseMS))
		Mild_Control=read.csv("Data_Covid/Recovery/Mild_Control.csv")
		Mild_RMD=read.csv("Data_Covid/Recovery/Mild_RMD.csv")
		
		Mild_Control_step=read.csv("Data_Covid/Recovery/step_Mild_Control.csv")
		Mild_RMD_step=read.csv("Data_Covid/Recovery/step_Mild_RMD.csv")
		
		Mild_Control_meanstep=read.csv("Data_Covid/Recovery/meanstep_Mild_Control.csv")
		Mild_RMD_meanstep=read.csv("Data_Covid/Recovery/meanstep_Mild_RMD.csv")
		
		d_diff=diff(Mild_Control$opt)
		
	}
	if ( caseMS=="severReady") { 
#		write.csv(cross5per_step2,sprintf("%scross5per1_%s.csv",outdir2,caseMS))
		Sever_Control=read.csv("Data_Covid/Recovery/Sever_Control.csv")
		Sever_RMD=read.csv("Data_Covid/Recovery/Sever_RMD.csv")
		
		Sever_Control_step=read.csv("Data_Covid/Recovery/step_Sever_Control.csv")
		Sever_RMD_step=read.csv("Data_Covid/Recovery/step_Sever_RMD.csv")
		
		Sever_Control_meanstep=read.csv("Data_Covid/Recovery/meanstep_Sever_Control.csv")
		Sever_RMD_meanstep=read.csv("Data_Covid/Recovery/meanstep_Sever_RMD.csv")
		d_diff=diff(Sever_Control$opt)
		
	}
	cross5per1$Roundedday=round(cross5per1$Recoveredt/24)
	cross5per1$NRoundedday=(cross5per1$Recoveredt/24)
	
	
	
	wantip1=c()
	
	version="Nsimple"
	if (version=="simple") {
		for (id1 in 1:11) {
			whip=which(abs(cross5per1$NRoundedday-id1)==min(abs(cross5per1$NRoundedday-id1)))
			howmany_want=round(d_diff[id1]*number_p)
			
			wantip=rep(cross5per1$ip[whip[1]],howmany_want)
			wantip1=c(wantip1,wantip)
			
		}
	}else{
		Er_round=0
		for (id1 in 1:27) {
			whip=which(cross5per1$Roundedday==id1)
			howmany_want=round(d_diff[id1]*number_p+Er_round)+round_correct 
			
			Er_round=(d_diff[id1]*number_p)-round(d_diff[id1]*number_p)
			
			if (length(whip)==0 | is.na(howmany_want)) {next}
			needP=howmany_want-length(whip)
			if (needP==0) {rep_ip=cross5per1$ip[whip]}
			if (needP>0) {rep_ip=c(cross5per1$ip[whip],rep(cross5per1$ip[whip[1]],needP))}
			if (needP<0) {rep_ip=cross5per1$ip[whip[1:howmany_want]]}
			
			wantip=rep_ip
			wantip1=c(wantip1,wantip)
		}
	}
	ucov_PN=number_p-length(wantip1)
	cross5per1=cross5per1[!is.na(cross5per1$NRoundedday),]
	
	simpleunc="NO"
	if (simpleunc=="yes") {
		ucov_P0=cross5per1$ip[!cross5per1$ip%in%wantip1]
		ucov_P=ucov_P0[ucov_P0!="NA"]
		set.seed(seedn)
		ucov_Pt1=sample(ucov_P,(ucov_PN-1),replace=TRUE)
		parnum=c(0,wantip1,ucov_Pt1)
		
	}else{
		
		set.seed(seedn)
		Uc2=round(xunc*number_p)
		Uc1=round(number_p-length(wantip1)-Uc2)*0
		deadman=0
		if (caseMS=="severReady") {deadman=0}
		
		Uc2=number_p-length(wantip1)-deadman
		if(Uc2<0) {Uc2=0}
		
		ucov_P=cross5per1$ip[cross5per1$NRoundedday>c1day & cross5per1$NRoundedday<c2day]
		set.seed(seedn)
		if (caseMS=="severReady"){
			
			ucov_P2=cross5per1$ip[cross5per1$NRoundedday>29  ]
			
			
		}
		if (caseMS=="mildReady"){
			
			ucov_P2=cross5per1$ip[cross5per1$NRoundedday>29 ]
		}
		
		if (length(ucov_P2)==0) {ucov_P2=ucov_P;}
		ucov_Pt2=sample(ucov_P2,Uc2,replace=TRUE) 
		if (Uc1>0) {
			ucov_Pt1=sample(ucov_P,Uc1,replace=TRUE)
			parnum=c(0,wantip1,ucov_Pt1,ucov_Pt2)
		}else{
			parnum=c(0,wantip1,ucov_Pt2)
		}
	}
	next 
}

for (ita in seq(0,47,1)) {
	Damageper0=sum(round(cross5per1$Recoveredt/24)<=ita)/(number_p+deadpop)*100
	Damageper1=rbind(Damageper1,c(ita*24,Damageper0))
}

Damageper2=data.frame(Damageper1)
colnames(Damageper2)=c("time","Recovered")
if (i1==1) {Damageper2placebo=Damageper2}
if (i1==2) {Damageper2RMD=Damageper2}



colnames(REcell)=c("time","qmin","qmax")
colnames(IL6q)=c("time","qmin","qmax")
colnames(FreeVirusq)=c("time","qmin","qmax")
colnames(IFNA1q)=c("time","qmin","qmax")
colnames(Tnaiveq)=c("time","qmin","qmax")
colnames(IgGq)=c("time","qmin","qmax")
colnames(IgMq)=c("time","qmin","qmax")

colnames(BEpitheliumLungq)=c("time","qmin","qmax")
colnames(BIL6q)=c("time","qmin","qmax")
colnames(BFreeVirusq)=c("time","qmin","qmax")
colnames(BIFNA1q)=c("time","qmin","qmax")
colnames(BTnaiveq)=c("time","qmin","qmax")
colnames(BIgGq)=c("time","qmin","qmax")
colnames(BIgMq)=c("time","qmin","qmax")



REcell=data.frame(REcell)
IL6q=data.frame(IL6q)
FreeVirusq=data.frame(FreeVirusq)
IFNA1q=data.frame(IFNA1q)
Tnaiveq=data.frame(Tnaiveq)
IgGq=data.frame(IgGq)
IgMq=data.frame(IgMq)

if (i1==1 ) {
	BEpitheliumLungq1=data.frame(BEpitheliumLungq)
	BIL6q1=data.frame(BIL6q)
	BFreeVirusq1=data.frame(BFreeVirusq)
	BIFNA1q1=data.frame(BIFNA1q)
	BTnaiveq1=data.frame(BTnaiveq)
	BIgGq1=data.frame(BIgGq)
	BIgMq1=data.frame(BIgMq)
}

if (i1==1) {
	REcell_NoRMD=REcell
	FreeVirusq_NoRMD=FreeVirusq
	
	
}



print("-------------------------***-------------------")

if (caserun=="Optimum") {
	source("opt.R")
} 
if (caserun=="Population" & i1==1) {
	source("population.R")
	
}



