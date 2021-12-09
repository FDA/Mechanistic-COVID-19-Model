print("next1")
source("random.R")
library(grid)
numsiz=16
labsiz=18
pars_Changed_path=outdiropt
BCinf_pars0=read.csv(paste0(pars_Changed_path,"viruscytokines_pars.txt"),sep = "",header = F,stringsAsFactors =FALSE)
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
	print(ip)
	if (ip==500 | ip==1000 | ip==1500 | ip==2000) {print(ip)}
	ipun=ipun+1
	source("parameters.R")	
#-------------------------
	BCinf_pars=pars
	BCinf_pars["timeout"]=30
	BCinf_pars["starttime"]=0
	BCinf_pars["unnamed"]=1
	
	initstates=states
	if (ip==0) {initstates["FreeInfluenza"]=CSpar["FreeInfluenza"]}else{initstates["FreeInfluenza"]=virusdiss[ip]}
#------------------------------------------------------------------------------
# READ EXPERIMENTAL DATA OF COVID-------------------------------
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
#
	VIRUS$time=(VIRUS$Rx*24+addday*24)
	IL6=read.csv(paste0("Data_Covid/IL6/",caseMS,".csv")); IL6$y1mean=IL6$y1mean*75/90 #nASAL IS USED FOR FITTING PLASMA*75=NASAL
	IL6=read.csv(paste0("Data_Covid/IL6/",caseMS,".csv"))
	IL6$y1mean=IL6$y1mean*75/90 #
	IL6$y1min=IL6$y1min*75/90 #
	IL6$y1max=IL6$y1max*75/90 #
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
#-----------eXPERIMENTAL ----------------------------
#--- load ODE model, states, and parameters
	isWindows="FALSE"	
	plotlist = list()
	mymodel<-list()
	source("PK.R")
	out_plot_paper=c()
	usedip[Lip]=ip
	ipcoef=runif(1000*(ipun+1),.98,1.02)
	ipcoef1=sample(ipcoef,length(pars))
	if (ip%in%usedip & i1>0) {pars=ipcoef1*pars;print(ipcoef1)}
	Lip=Lip+1
	pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
	pars["timeout"]=30
	pars["unnamed"]=1
#--------------------------------
	try({out_plot_paper <- dede(initstates, timepoints, "derivs", pars, dllname=modelname,initfunc="initmod",rtol=1e-3,atol=1e-6,nout=3,events=list(data=eventdata))});		
	out_plot_paper=data.frame(out_plot_paper)
	out_plot_paper$ip=ip	
	out_plot_paper$ipun=ipun
	out_plot_paper$initVirus=(out_plot_paper$FreeInfluenza)[out_plot_paper$time==20]
	out_plot_paper$maxIL6=max(out_plot_paper$IL6)
	out_plot_paper$damage=max(4E8-out_plot_paper$ResistantEpithelium-out_plot_paper$EpitheliumLung)/4E8*100
	out_plot_paper$damage0=(4E8-out_plot_paper$ResistantEpithelium-out_plot_paper$EpitheliumLung)/4E8*100
	out_plot_paper$maxFreeInfluenza=max(out_plot_paper$FreeInfluenza)
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
	out_plot_paper[,"time"]=out_plot_paper[,"time"]-daybase #one day befor first injection
	out_plot_paper=out_plot_paper[out_plot_paper[,"time"]>=0,]
	if (ip ==0 ) {out_plot1_opt0=out_plot_paper}
	out_plot_paper$minVir=min(out_plot_paper$FreeInfluenza)
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
		minDamageBeforInj=min(out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time==0]) #day 0 is the inection day
		if (minDamageBeforInj<5) {cross5pert0=-10000}

		if (is.na(cross5pert0)) {cross5pert0=10000}
		#-----------lets check what is the damage at day 28 -----------------------
	    damageD0=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*0][1]
		damageD1=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*1][1]
		damageD10=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*10][1]
		damageD11=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*11][1]
		damageD12=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*12][1]
		damageD15=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*15][1]
		damageD28=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*28][1]
	    #--------------------------------------------------------------
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
	IL6q=IgMq=IgGq=Tnaiveq=IFNA1q=FreeInfluenzaq=REcell=Damageper1=c()
	out_plot1_paper$REcell=out_plot1_paper$ResistantEpithelium+out_plot1_paper$EpitheliumLung
	for (ti1 in seq(0,(max(out_plot1_paper[,"time"])-24),2)) {
		tout_plot1=out_plot1_paper[out_plot1_paper[,1]==ti1,]
		tq=(lapply(tout_plot1,quantile,probs=c(0.05,.95),na.rm=TRUE))
		REcell=rbind(REcell,c(ti1,tq$REcell))
		IL6q=rbind(IL6q,c(ti1,tq$IL6))
		FreeInfluenzaq=rbind(FreeInfluenzaq,c(ti1,tq$FreeInfluenza))
		IFNA1q=rbind(IFNA1q,c(ti1,tq$IFNA1))
		Tnaiveq=rbind(Tnaiveq,c(ti1,tq$Tnaive))
		IgGq=rbind(IgGq,c(ti1,tq$IgG))
		IgMq=rbind(IgMq,c(ti1,tq$IgM))
	}
	BIL6q=BIgMq=BIgGq=BTnaiveq=BIFNA1q=BFreeInfluenzaq=BREcell=BEpitheliumLungq=c()
	Bioout_plot_paperall$FreeInfluenza=log10(Bioout_plot_paperall$FreeInfluenza)
	for (ti1 in seq(0,max(out_plot1[,"time"]),2)) { 
		tout_plot1=Bioout_plot_paperall[Bioout_plot_paperall[,1]==ti1,]
		tq=(lapply(tout_plot1,quantile,probs=c(0.1,.9),na.rm=TRUE))
		BIL6q=rbind(BIL6q,c(ti1,tq$IL6))
		BFreeInfluenzaq=rbind(BFreeInfluenzaq,c(ti1,tq$FreeInfluenza))
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
if (i1==1) {write.csv(cross5per1,sprintf("%s_i1is1_cross5per1_%s.csv",outdir2,caseMS))}
if (i1==2) {write.csv(cross5per1,sprintf("%s_i1is2_cross5per1_%s.csv",outdir2,caseMS))}
if (i1==0) {
	if (caseMS=="mildReady") { 
		write.csv(cross5per_step2,sprintf("%s_cross5per1_%s.csv",outdir2,caseMS))
		Mild_Control=read.csv("Data_Covid/Recovery/Mild_Control.csv")
		Mild_RMD=read.csv("Data_Covid/Recovery/Mild_RMD.csv")
		
		Mild_Control_step=read.csv("Data_Covid/Recovery/step_Mild_Control.csv")
		Mild_RMD_step=read.csv("Data_Covid/Recovery/step_Mild_RMD.csv")
		
		Mild_Control_meanstep=read.csv("Data_Covid/Recovery/meanstep_Mild_Control.csv")
		Mild_RMD_meanstep=read.csv("Data_Covid/Recovery/meanstep_Mild_RMD.csv")
		
		d_diff=diff(Mild_Control$opt)

	}
	if ( caseMS=="severReady") { 
		write.csv(cross5per_step2,sprintf("%scross5per1_%s.csv",outdir2,caseMS))
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
			print(howmany_want)
			wantip=rep(cross5per1$ip[whip[1]],howmany_want)
			wantip1=c(wantip1,wantip)
		}
	}else{
		Er_round=0
		for (id1 in 1:27) {
			whip=which(cross5per1$Roundedday==id1)
			howmany_want=round(d_diff[id1]*number_p+Er_round)+round_correct #number of patinet recovers in each day based on number_p
			
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
		print(ucov_P)
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
		ucov_P=cross5per1$ip[cross5per1$NRoundedday>c1day & cross5per1$NRoundedday<c2day]
		set.seed(seedn)
	if (caseMS=="severReady"){
			ucov_P2=cross5per1$ip[cross5per1$NRoundedday>29  ]
	}
	if (caseMS=="mildReady"){
		ucov_P2=cross5per1$ip[cross5per1$NRoundedday>29 ]
	}
		if (length(ucov_P2)==0) {ucov_P2=ucov_P; print("length ucov_P2 is zero");}
		ucov_Pt2=sample(ucov_P2,Uc2,replace=TRUE) #they never recovered
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
colnames(FreeInfluenzaq)=c("time","qmin","qmax")
colnames(IFNA1q)=c("time","qmin","qmax")
colnames(Tnaiveq)=c("time","qmin","qmax")
colnames(IgGq)=c("time","qmin","qmax")
colnames(IgMq)=c("time","qmin","qmax")
#
colnames(BEpitheliumLungq)=c("time","qmin","qmax")
colnames(BIL6q)=c("time","qmin","qmax")
colnames(BFreeInfluenzaq)=c("time","qmin","qmax")
colnames(BIFNA1q)=c("time","qmin","qmax")
colnames(BTnaiveq)=c("time","qmin","qmax")
colnames(BIgGq)=c("time","qmin","qmax")
colnames(BIgMq)=c("time","qmin","qmax")



REcell=data.frame(REcell)
IL6q=data.frame(IL6q)
FreeInfluenzaq=data.frame(FreeInfluenzaq)
IFNA1q=data.frame(IFNA1q)
Tnaiveq=data.frame(Tnaiveq)
IgGq=data.frame(IgGq)
IgMq=data.frame(IgMq)
#--------bio
if (i1==1 ) {
	BEpitheliumLungq1=data.frame(BEpitheliumLungq)
	BIL6q1=data.frame(BIL6q)
	BFreeInfluenzaq1=data.frame(BFreeInfluenzaq)
	BIFNA1q1=data.frame(BIFNA1q)
	BTnaiveq1=data.frame(BTnaiveq)
	BIgGq1=data.frame(BIgGq)
	BIgMq1=data.frame(BIgMq)
}

if (i1==1) {
	REcell_NoRMD=REcell
	FreeInfluenzaq_NoRMD=FreeInfluenzaq
	
	
}


#-------------------------***-------------------
point_size=3
line_size=1.5
out_plot1$yPred_TNAIVE_perc=(10^5-out_plot1$Tnaive)/10^5*100
BTnaiveq2=BTnaiveq1
BTnaiveq2$qmin=(10^5-BTnaiveq1$qmax)/10^5*100
BTnaiveq2$qmax=(10^5-BTnaiveq1$qmin)/10^5*100
p<-ggplot()
p<-p+geom_point(data=TNAIVE_percent1plot, aes( x=time,y=Tnaivepercent), size=point_size, alpha=1, color="black")
p<-p+geom_line(data=out_plot1, aes(x=time, y=yPred_TNAIVE_perc), size=line_size, alpha=1,linetype = "solid")
p<-p+geom_errorbar(data=ErrbarTn,aes(x=time,ymin=ymin, ymax=ymax),color="black",linetype = "solid",width=20,
		position=position_dodge(20))#, width=.2,position=position_dodge(.9))
p<-p+geom_ribbon(data=BTnaiveq2, aes(x=time, ymin=qmax,ymax=qmin),alpha=0.3)
p<-p+ylab(paste0("Decrease % of \n T lymphocytes"))+ylim(0,100)+xlim(0,576)+ theme_bw()
p<-p+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
BIL6q=data.frame(BIL6q)
BIL6q$qmin=BIL6q$qmin/75
BIL6q$qmax=BIL6q$qmax/75
BIL6q$time=BIL6q$time
out_plot1$IL6=out_plot1$IL6/75
IL6$y1min[IL6$y1min<0]=0
p4<-ggplot()
p4<-p4+geom_point(data=IL6, aes(x=time, y=y1mean), size=point_size, alpha=1, color="black")
p4<-p4+geom_line(data=out_plot1, aes(x=time, y=IL6), size=line_size, alpha=1,linetype = "solid")
p4<-p4+geom_ribbon(data=BIL6q, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)
p4<-p4+geom_errorbar(data=IL6,aes(x=time,ymin=y1min, ymax=y1max),width=20,
position=position_dodge(20),col="black",na.rm=TRUE)#)
p4<-p4+ylab(paste0("IL6 \n (pg/ml)"))+ylim(0,175)+xlim(0,576)+ theme_bw()
p4<-p4+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
VIRUS$FreeInfluenza1=log10(VIRUS$yP)
out_plot1$FreeInfluenza=log10(out_plot1$FreeInfluenza)
out_plot2$FreeInfluenza=log10(out_plot2$FreeInfluenza)
p1<-ggplot()
p1<-p1+geom_line(data=out_plot1, aes(x=time, y=FreeInfluenza), size=line_size, alpha=1,linetype = "solid")
p1<-p1+geom_point(data=VIRUS, aes(x=time, y=FreeInfluenza1), size=point_size, alpha=1, color="black")
p1<-p1+geom_ribbon(data=BFreeInfluenzaq1, aes(x=time, ymin=(qmin),ymax=(qmax)),alpha=0.3)
p1<-p1+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title=element_text(size=15,face="bold"))+ theme_bw()
p1<-p1+ylab("Virus Titer \n (log[copy num/ml])")+ylim(0,9)+xlim(0,576)+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
IFNAalpha0<-read.csv(paste0("Data_Covid/IFNAalpha/","IFNA_Fig1A_withmaxmin",".csv"),sep = ",") #read data
IFNAalpha=data.frame(cbind(unique(IFNAalpha0$Day),unique(IFNAalpha0$INFA_mean_pgperml),unique(IFNAalpha0$min),unique(IFNAalpha0$max))) #Data in each day 
colnames(IFNAalpha)=c("Day","INFA_mean_pgperml","Imin","Imax") #paper was based on fg/ml , i changed it
IFNAalpha=IFNAalpha[IFNAalpha$Day!=8,]
IFNAalpha$INFA_mean_pgperml=12*IFNAalpha$INFA_mean_pgperml    
IFNAalpha$Imin=12*IFNAalpha$Imin
IFNAalpha$Imax=12*IFNAalpha$Imax
#--------Sever--
IFA_sever=read.csv("Data_Covid/IFNAalpha/IFA_sever.csv")
IFA_sever=data.frame(IFA_sever)
IFA_sever$IFN_S=IFA_sever$IFN_S*12*.1
#------end sever-------
IFNAalpha=IFNAalpha[IFNAalpha$Day<23,] 		#how many days after onset is keep (i used 15day) 
IFNAalpha$Time=(IFNAalpha$Day+addday)*24
IFA_sever=IFA_sever[IFA_sever$Day<23,] 		#how many days after onset is keep (i used 15day) 
IFA_sever$Time=(IFA_sever$Day+addday)*24
if (IFNAalpha$Imax[2]>200) {IFNAalpha$Imax[2]=200}
p5<-ggplot()
if (caseMS=="mildReady") {
	p5<-p5+geom_point(data=IFNAalpha, aes(x=Time, y=INFA_mean_pgperml), size=point_size, alpha=1, color="black")
	p5<-p5+geom_errorbar(data=IFNAalpha,aes(x=Time,ymin=Imin, ymax=Imax),color="black",linetype = "solid",na.rm=TRUE, width=20,
			position=position_dodge(20))#, width=.2,position=position_dodge(.9))
}
if (caseMS=="severReady") {
	p5<-p5+geom_point(data=IFA_sever, aes(x=Time, y=IFN_S), size=point_size, alpha=1, color="black")
}
lim=c(-5,200)
p5<-p5+geom_line(data=out_plot1, aes(x=time, y=IFNA1), size=line_size, alpha=1,linetype = "solid")
p5<-p5+geom_ribbon(data=BIFNA1q1, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)
p5<-p5+ylab(paste0("IFN Alpha \n (pg/ml)"))+ylim(0,200)+xlim(0,576)+scale_y_continuous(expand = c(0, 0), limits = lim) +
		theme_bw()
p5<-p5+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p6<-ggplot()
p6<-p6+geom_line(data=out_plot1, aes(x=time, y=TcLymphnode), size=0.85, alpha=0.8,linetype = "dashed")
p6<-p6+geom_ribbon(data=BIFNA1q1, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)
p6<-p6+ylab(paste0("TcLymphnode"))+ theme_bw()
p6<-p6+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p7<-ggplot()
p7<-p7+geom_line(data=out_plot1, aes(x=time, y=EpitheliumLung ), size=0.85, alpha=0.8,linetype = "solid")
p7<-p7+geom_ribbon(data=BEpitheliumLungq1, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)
p7<-p7+ylab(paste0("EpitheliumLung"))+ theme_bw()
p7<-p7+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p8<-ggplot(out_plot1, aes(x=time,y=ResistantEpithelium ))
K11=pars["K11"]
K26=pars["K26"]
K288=pars["K288"]
K18=pars["K18"]
#p4<-p4+geom_point(data=IL6, aes(x=time, y=y1mean), size=line_size, alpha=.5, color="green")
p8<-p8+geom_line(data=out_plot1, aes(x=time, y=K11*FreeInfluenza*EpitheliumLung ), size=0.85, alpha=0.8,linetype = "solid",color="red")
p8<-p8+geom_line(data=out_plot1, aes(x=time, y=K26*IFNA1*EpitheliumLung ), size=0.85, alpha=0.8,linetype = "solid",color="black")
p8<-p8+geom_line(data=out_plot1, aes(x=time, y=K288*IL6*EpitheliumLung ), size=0.85, alpha=0.8,linetype = "solid",color="green")
p8<-p8+geom_line(data=out_plot1, aes(x=time, y=K18*FreeInfluenza ), size=0.85, alpha=0.8,linetype = "solid",color="pink")
p8<-p8+ylab(paste0("ResistantEpithelium"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
if (caseMS=="severReady") {
NSIgM_data<-read.csv(paste0("Data_Covid/IgMIgG/","SIgM_Ready",".csv"),sep = ",") #read data non-S IgM
NSIgM_data=NSIgM_data[NSIgM_data$day<23,] 		#how many days after onset is keep (i used 15day) 
NSIgM_data$Time=(NSIgM_data$day+addday)*24    #addday for mild case is 7 days. days after infection
NSIgG_data<-read.csv(paste0("Data_Covid/IgMIgG/","SIgG_Ready",".csv"),sep = ",") #read data non-S IgG
NSIgG_data=NSIgG_data[NSIgG_data$day<23,] 		#how many days after onset is keep (i used 15day) 
NSIgG_data$Time=(NSIgG_data$day+addday)*24    #addday for mild case is 7 days. days after infection
}
if (caseMS=="mildReady") {
	NSIgM_data<-read.csv(paste0("Data_Covid/IgMIgG/","NSIgM_Ready",".csv"),sep = ",") #read data non-S IgM
	NSIgM_data=NSIgM_data[NSIgM_data$day<23,] 		#how many days after onset is keep (i used 15day) 
	NSIgM_data$Time=(NSIgM_data$day+addday)*24    #addday for mild case is 7 days. days after infection
	NSIgG_data<-read.csv(paste0("Data_Covid/IgMIgG/","NSIgG_Ready",".csv"),sep = ",") #read data non-S IgG
	NSIgG_data=NSIgG_data[NSIgG_data$day<23,] 		#how many days after onset is keep (i used 15day) 
	NSIgG_data$Time=(NSIgG_data$day+addday)*24    #addday for mild case is 7 days. days after infection
	
}

p9<-ggplot()
p9<-p9+geom_point(data=NSIgM_data, aes(x=Time, y=mean), size=point_size, alpha=1, color="black")
p9<-p9+geom_errorbar(data=NSIgM_data,aes(x=Time,ymin=min, ymax=max),color="black",linetype = "solid",na.rm=TRUE,width=20,
		position=position_dodge(20))#, width=.2,position=position_dodge(.9))
p9<-p9+geom_line(data=out_plot1, aes(x=time, y=IgM), size=line_size, alpha=1,linetype = "solid")
p9<-p9+geom_ribbon(data=BIgMq1, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)+ylim(0,690)+xlim(0,576)

p9<-p9+ylab(paste0("Titer of IgM \n "))+ theme_bw()
p9<-p9+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title=element_text(size=20,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())


p10<-ggplot()
p10<-p10+geom_point(data=NSIgG_data, aes(x=Time, y=mean), size=point_size, alpha=1, color="black")
p10<-p10+geom_line(data=out_plot1, aes(x=time, y=IgG), size=line_size, alpha=1,linetype = "solid")
p10<-p10+geom_errorbar(data=NSIgG_data,aes(x=Time,ymin=min, ymax=max),color="black",linetype = "solid",na.rm=TRUE,width=20,
		position=position_dodge(20))#, width=.2,position=position_dodge(.9))
p10<-p10+geom_ribbon(data=BIgGq1, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)
p10<-p10+ylab(paste0("Titer of IgG \n "))+ylim(0,69)+xlim(0,576)+ theme_bw()
p10<-p10+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title=element_text(size=20,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())

plotlist[[1]] = p
plotlist[[2]] = p1
plotlist[[3]] = p4
plotlist[[4]] = p5
plotlist[[5]] = p6
library(grid)
pall=grid.arrange(arrangeGrob(p1,p5,p4,p,p9,p10, ncol=1))
biologicalfig="biological_figs6"
system(paste0("mkdir -p ",outdir2,biologicalfig))		
ggsave(paste0(outdir2,biologicalfig,"/",byplot,caseMS,"_7.png"),pall,height=16, width=8)
wsi=4
ggsave(paste0(outdir2,biologicalfig,"/",caseMS,"_Tnaive.png"), p, width=wsi, height=4)
ggsave(paste0(outdir2,biologicalfig,"/",caseMS,"_IL6.png"), p4, width=wsi, height=4)
ggsave(paste0(outdir2,biologicalfig,"/",caseMS,"_Virus.png"), p1, width=wsi, height=4)
ggsave(paste0(outdir2,biologicalfig,"/",caseMS,"_INFA.png"), p5, width=wsi, height=4)
ggsave(paste0(outdir2,biologicalfig,"/",caseMS,"_IgM.png"), p9, width=wsi, height=4)
ggsave(paste0(outdir2,biologicalfig,"/",caseMS,"_IgG.png"), p10, width=wsi, height=4)
PVIF<-ggplot()
PVIF<-PVIF+geom_line(data=out_plot1, aes(x=time, y=FreeInfluenza), size=line_size, alpha=1,linetype = "solid")
PVIF<-PVIF+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title=element_text(size=15,face="bold"))+ theme_bw()
PVIF<-PVIF+ylab("Virus Titer (log[copy num/ml]")+ylim(0,9)+xlim(0,576)+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=13,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
PVIF1<-ggplot()
lim=c(-5,200)
PVIF1<-PVIF1+geom_line(data=out_plot1, aes(x=time, y=IFNA1), size=line_size, alpha=1,linetype = "solid",color="darkgray")
PVIF1<-PVIF1+ylab(paste0("INF Alpha (pg/ml)"))+ylim(0,200)+xlim(0,576)+scale_y_continuous(expand = c(0, 0), limits = lim) +
		theme_bw()
PVIF1<-PVIF1+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=13,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
PVIF1<-PVIF1+scale_y_continuous(position = "right")
ggsave(paste0(outdir2,biologicalfig,"/",byplot,caseMS,"_V.png"),PVIF,height=3.58, width=3.38)
ggsave(paste0(outdir2,biologicalfig,"/",byplot,caseMS,"_IF.png"),PVIF1,height=3.58, width=3.38)
#----influenza Virus vs IFN pick -------------------------------------------------------
expINF=read.csv("Data_Covid/InfluenzaPick/expINF.csv",header=F) #no need
IFNP=read.csv("Data_Covid/InfluenzaPick/IFNP.csv",header=F)
expV=read.csv("Data_Covid/InfluenzaPick/expV.csv",header=F) #no need 
VirusP=read.csv("Data_Covid/InfluenzaPick/VirusP.csv",header=F)
iPVIF<-ggplot()
iPVIF<-iPVIF+geom_line(data=VirusP, aes(x=V1, y=V2), size=line_size, alpha=1,linetype = "solid")
iPVIF<-iPVIF+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title=element_text(size=15,face="bold"))+ theme_bw()
iPVIF<-iPVIF+ylab("Virus Titer (log[copy num/ml]")+ylim(0,6)+xlim(0,192)+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=13,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
iPVIF1<-ggplot()
lim=c(-5,200)
iPVIF1<-iPVIF1+geom_line(data=IFNP, aes(x=V1, y=V2), size=line_size, alpha=1,linetype = "solid",color="darkgray")
iPVIF1<-iPVIF1+ylab(paste0("INF Alpha (pg/ml)"))+ylim(0,300)+xlim(0,192)+scale_y_continuous(expand = c(0, 0), limits = lim) +
		theme_bw()
iPVIF1<-iPVIF1+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=13,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank(),panel.grid.minor = element_blank(),panel.grid.major = element_blank())
iPVIF1<-iPVIF1+scale_y_continuous(position = "right")
ggsave(paste0(outdir2,biologicalfig,"/",byplot,caseMS,"_Vi.png"),iPVIF,height=3.58, width=3.38)
ggsave(paste0(outdir2,biologicalfig,"/",byplot,caseMS,"_IFi.png"),iPVIF1,height=3.58, width=3.38)
totalkilled<- 4E8 - out_plot1[,"EpitheliumLung"] -out_plot1[,"ResistantEpithelium"]
IL6effecterror<-min((out_plot1[,"EpiKilledbyIL6"]/totalkilled)[dim(out_plot1)[1]]-0.4,0)
print("IL6effecterror")
print(IL6effecterror)
print("EpiKilledbyIL6/totalkilled")
print((out_plot1[,"EpiKilledbyIL6"]/totalkilled)[dim(out_plot1)[1]])
infectedrate<-(4e8 - out_plot1[,colnames(out_plot1)=="EpitheliumLung"] - out_plot1[, colnames(out_plot1)=="ResistantEpithelium"])/4e8
infectionerror<- max(max(infectedrate)[1] - 0.7,0)^2 + min(max(infectedrate)[1] - 0.4,0)^2
print("infectedrate")
print(max(infectedrate))
print(min(infectedrate))
ppk=ggplot() +
		geom_line(data=out_plot2,aes(x=time, y=CR/10^6)) +
		ggtitle(paste0(caseMS))+theme(legend.position = "none") +
		labs(x ="time (hr)", y = "Central comp (mg)")#+xlim(0, (st_day+num_day)*24)

if (i1==1) {out_plot1i1=out_plot2}
if (i1==2) {out_plot1i2=out_plot2}
if (i1==1) {out_plot1i1L=out_plot_paperall}
if (i1==2) {out_plot1i2L=out_plot_paperall}

if (i1==2) {

	p1<-ggplot()

	p1<-p1+geom_line(data=out_plot1i1, aes(x=time, y=FreeInfluenza), size=0.85, alpha=0.8,linetype = "solid")
	p1<-p1+geom_line(data=out_plot1i2, aes(x=time, y=FreeInfluenza), size=0.85, alpha=0.8,linetype = "dashed",col="green")
	p1<-p1+geom_ribbon(data=FreeInfluenzaq, aes(x=time, ymin=log10(qmin),ymax=log10(qmax)),alpha=0.3)
	p1<-p1+geom_ribbon(data=FreeInfluenzaq_NoRMD, aes(x=time, ymin=log10(qmin),ymax=log10(qmax)),alpha=0.3,color="brown",fill="brown")
	
	p1<-p1+xlab(paste0("time (hr)"))
	p1<-p1+ylab("Virus log")+ggtitle(paste0(caseMS,": \n brown= without drug, Gray= with drug"))+xlim(0,15*24)#+xlim(0, (st_day+num_day)*24)
	ppk1=ggplot() +
			geom_ribbon(data=REcell, aes(x=time/24, ymin=(4E8-qmax)/4E8*100,ymax=(4E8-qmin)/4E8*100) ,alpha=0.8,fill="gray")+
			geom_hline(yintercept=5, linetype="dashed", color = "black")+
			theme(legend.position = "none") +theme_bw()+
			labs(x ="Time(day)", y = "%Damaged Cell")+xlim(0,15)+ylim(0,100)#+xlim(0, (st_day+num_day)*24)
		plotlist1=list()
	plotlist1[[1]] = ppk1
	pall1 <- grid.arrange(grobs=plotlist1,ncol=1)
	ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,".png"), pall1, width=3.38, height=3.38)

ppkjc1=ggplot() +
		geom_ribbon(data=REcell_NoRMD, aes(x=time/24, ymin=(4E8-qmax)/4E8*100,ymax=(4E8-qmin)/4E8*100) ,alpha=0.8,fill="gray")+
		geom_hline(yintercept=5, linetype="dashed", color = "black")+
		theme(legend.position = "none") + theme_bw()+
		labs(x ="Time(day)", y = "%Damaged Cell")+xlim(0,15)+ylim(0,100)#+xlim(0, (st_day+num_day)*24)

plotlist1=list()
plotlist1[[1]] = ppkjc1
pall1 <- grid.arrange(grobs=plotlist1,ncol=1)
ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,"_justcontrol.png"), pall1, width=3.38, height=3.38)
ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,"_justcontrol1.png"), ppkjc1, width=3.38, height=3.38)
	out_plot1_opt0Biolog_opt_alltime$FreeInfluenza=log10(out_plot1_opt0Biolog_opt_alltime$FreeInfluenza)
	out_plot1_opt0_opt_alltime$FreeInfluenza=log10(out_plot1_opt0_opt_alltime$FreeInfluenza)
	optp1<-ggplot()
	optp1<-optp1+geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=FreeInfluenza), size=0.85, alpha=0.8,linetype = "solid")
	optp1<-optp1+geom_line(data=out_plot1_opt0_opt_alltime, aes(x=time, y=FreeInfluenza), size=0.85, alpha=0.8,linetype = "dashed",col="green")
	optp1<-optp1+xlab(paste0("time (hr)"))
	optp1<-optp1+ylab("Virus log")+ggtitle(paste0(caseMS,": \n Black= without drug, Green= with drug"))+xlim(0,15*24)#+xlim(0, (st_day+num_day)*24)

	
	
	optppk1=ggplot() +
			geom_line(data=out_plot1_opt0Biolog_opt_alltime,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100)) +
			geom_line(data=out_plot1_opt0_opt_alltime,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100),color="green",linetype = "dashed") +
			theme(legend.position = "none") +
			labs(x ="time (hr)", y = "%Damaged Cell")+xlim(0,15*24)
	plotlist2=list()
	plotlist2[[1]] = optp1
	plotlist2[[2]] = optppk1
	pall2 <- grid.arrange(grobs=plotlist2,ncol=1)
	ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,"_opt.png"), pall2, width=6, height=8)
#-----------------
	totalkilled<- out_plot1[,"EpiKilledbyIFN"]+out_plot1[,"EpiKilledbyIL6"]+out_plot1[,"EpiKilledbyTc"]+out_plot1[,"EpiKilledbyVirus"]
	totalkilled=1
	kia1<-ggplot()+
			geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=EpiKilledbyIFN/totalkilled ), size=0.85, alpha=0.8,linetype = "solid")+
			geom_line(data=out_plot1_opt0_opt_alltime, aes(x=time, y=EpiKilledbyIFN/totalkilled ), size=0.85, alpha=0.8,linetype = "solid",color="green")+
			ylab(paste0("EpiKilledbyIFN"))
	
	kia2<-ggplot()+
			geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=EpiKilledbyIL6/totalkilled ), size=0.85, alpha=0.8,linetype = "solid")+
			geom_line(data=out_plot1_opt0_opt_alltime, aes(x=time, y=EpiKilledbyIL6/totalkilled ), size=0.85, alpha=0.8,linetype = "solid",color="green")+
			ylab(paste0("EpiKilledbyIL6"))
	kia3<-ggplot()+
			geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=EpiKilledbyTc/totalkilled ), size=0.85, alpha=0.8,linetype = "solid")+
			geom_line(data=out_plot1_opt0_opt_alltime, aes(x=time, y=EpiKilledbyTc/totalkilled ), size=0.85, alpha=0.8,linetype = "solid",color="green")+
			ylab(paste0("EpiKilledbyTc"))
	
	kia4<-ggplot()+

			geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=EpiKilledbyVirus ), size=0.85, alpha=0.8,linetype = "solid")+
			geom_line(data=out_plot1_opt0_opt_alltime, aes(x=time, y=EpiKilledbyVirus ), size=0.85, alpha=0.8,linetype = "solid",color="green")+
			
			ylab(paste0("EpiKilledbyVirus_",EC,"_sever"))
	
	plotlist2k = list()
	
	plotlist2k[[1]] = kia1
	plotlist2k[[2]] = kia2
	plotlist2k[[3]] = kia3
	plotlist2k[[4]] = kia4
#plotlist[[5]] = p6
	
	pallkk <- grid.arrange(grobs=plotlist2k,ncol=2)
	ggsave(paste0(outdir2,byplot,"Sever_Kill4.png"), pallkk, width=8, height=8)
	totalkilled<- out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyIFN"]+out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyIL6"]+out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyTc"]+out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyVirus"]
	out_plot1_opt0Biolog_opt_alltime$totalkilled_t=totalkilled[out_plot1_opt0Biolog_opt_alltime[,"time"]==tki]-totalkilled[out_plot1_opt0Biolog_opt_alltime[,"time"]==(tki-2)]
	out_plot1_opt0Biolog_opt_alltime$EpiKilledbyIFN_t=out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyIFN"][out_plot1_opt0Biolog_opt_alltime[,"time"]==tki]-out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyIFN"][out_plot1_opt0Biolog_opt_alltime[,"time"]==(tki-2)]
	out_plot1_opt0Biolog_opt_alltime$EpiKilledbyIL6_t=out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyIL6"][out_plot1_opt0Biolog_opt_alltime[,"time"]==tki]-out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyIL6"][out_plot1_opt0Biolog_opt_alltime[,"time"]==(tki-2)]
	out_plot1_opt0Biolog_opt_alltime$EpiKilledbyTc_t=out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyTc"][out_plot1_opt0Biolog_opt_alltime[,"time"]==tki]-out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyTc"][out_plot1_opt0Biolog_opt_alltime[,"time"]==(tki-2)]
	out_plot1_opt0Biolog_opt_alltime$EpiKilledbyVirus_t=out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyVirus"][out_plot1_opt0Biolog_opt_alltime[,"time"]==tki]-out_plot1_opt0Biolog_opt_alltime[,"EpiKilledbyVirus"][out_plot1_opt0Biolog_opt_alltime[,"time"]==(tki-2)]
	ki1t<-ggplot()+
			geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=EpiKilledbyIFN_t/totalkilled_t ), size=0.85, alpha=0.8,linetype = "solid")+
			ylab(paste0("EpiKilledbyIFN"))
	
	ki2t<-ggplot()+
			geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=EpiKilledbyIL6_t/totalkilled_t ), size=0.85, alpha=0.8,linetype = "solid")+
			ylab(paste0("EpiKilledbyIL6"))
	
	ki3t<-ggplot()+
			geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=EpiKilledbyTc_t/totalkilled_t ), size=0.85, alpha=0.8,linetype = "solid")+
			ylab(paste0("EpiKilledbyTc"))
	
	ki4t<-ggplot()+
			geom_line(data=out_plot1_opt0Biolog_opt_alltime, aes(x=time, y=EpiKilledbyVirus_t/totalkilled_t ), size=0.85, alpha=0.8,linetype = "solid")+
			ylab(paste0("EpiKilledbyVirus_",EC,"_sever"))
	
	plotlist2t = list()
	
	plotlist2t[[1]] = ki1t
	plotlist2t[[2]] = ki2t
	plotlist2t[[3]] = ki3t
	plotlist2t[[4]] = ki4t

	pallkt <- grid.arrange(grobs=plotlist2t,ncol=2)
	ggsave(paste0(outdir2,byplot,"Sever_Kill_t.png"), pallkt, width=8, height=8)

	All_adm_t=data.frame(All_adm_t)
	thist=ggplot(data.frame(x=All_adm_t$hr),aes(x,y=..density..)) + 
			geom_histogram(fill=NA,color="black")+
			geom_density(color="red",linetype="dotted")+
			labs(x ="adm_time(hr)", y = "density")
	
	all_V_int=out_plot1_paper$FreeInfluenza[out_plot1_paper$time==0]
	thisV=ggplot(data.frame(x=log10(all_V_int)),aes(x,y=..density..)) + 
			geom_histogram(fill=NA,color="black")+
			geom_density(color="red",linetype="dotted")+
			labs(x ="init_V", y = "density")
	
	plotlist3=list()
	plotlist3[[1]] = thist
	plotlist3[[2]] = thisV
	pall3 <- grid.arrange(grobs=plotlist3,ncol=1)
	ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,"_distribution.png"), pall3, width=6, height=8)
	optp1L<-ggplot()
	
	optp1L<-optp1L+geom_line(data=out_plot1i1L, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "solid")
	optp1L<-optp1L+geom_line(data=out_plot1i2L, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "dashed",col="red")
	optp1L<-optp1L+xlab(paste0("time (hr)"))
	optp1L<-optp1L+ylab("Virus log")+ggtitle(paste0(caseMS,": \n Black= without drug, Green= with drug"))#+xlim(0, (st_day+num_day)*24)
	
	optL=ggplot() +
			geom_line(data=out_plot1i1L,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip)) +
			geom_line(data=out_plot1i2L,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip),color="red",linetype="dashed") +
			theme(legend.position = "none") +
			labs(x ="time (hr)", y = "%Damaged Cell")#+xlim(0, (st_day+num_day)*24)
	
	plotlist4=list()
	plotlist4[[2]] = optL
	plotlist4[[1]] = optp1L
	pall4 <- grid.arrange(grobs=plotlist4,ncol=1)
	ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,"_L.png"), pall4, width=6, height=8)
	sepday=12
	out_plot1i1LSH=out_plot1i1L[out_plot1i1L$eventtimes<=sepday*24,]
	out_plot1i2LSH=out_plot1i2L[out_plot1i2L$eventtimes<=sepday*24,]
	
	out_plot1i1LL=out_plot1i1L[out_plot1i1L$eventtimes>sepday*24,]
	out_plot1i2LL=out_plot1i2L[out_plot1i2L$eventtimes>sepday*24,]
	
	optp1LSH<-ggplot()+
			geom_line(data=out_plot1i1LSH, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "solid")+
			geom_line(data=out_plot1i2LSH, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "dashed",col="red")+
			geom_line(data=out_plot1i1LL, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "solid",col="blue")+
			geom_line(data=out_plot1i2LL, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "dashed",col="green")+
			xlab(paste0("time (hr)"))+
			ylab("Virus log")+ggtitle(paste0(caseMS,": \n Black= without drug, Green= with drug"))#+xlim(0, (st_day+num_day)*24)
	
	optLLSH=ggplot() +
			geom_line(data=out_plot1i1LSH,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip)) +
			geom_line(data=out_plot1i2LSH,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip),color="red",linetype="dashed") +
			geom_line(data=out_plot1i1LL,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip),col="blue") +
			geom_line(data=out_plot1i2LL,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip))+
			theme(legend.position = "none") +
			labs(x ="time (hr)", y = "%Damaged Cell") #+xlim(0, (st_day+num_day)*24)
	
	plotlist5=list()
	
	plotlist5[[1]] = optp1LSH
	plotlist5[[2]] = optLLSH
	pall5 <- grid.arrange(grobs=plotlist5,ncol=1)
	ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,"_LOSH.png"), pall5, width=6, height=8)
	
	optp1LSH1<-ggplot()+
			geom_line(data=out_plot1i1LSH, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "solid")+
			geom_line(data=out_plot1i2LSH, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "dashed",col="red")+
			xlab(paste0("time (hr)"))+
			ylab("Virus log_short")#+ggtitle(paste0(caseMS,": \n Black= without drug, Green= with drug"))#+xlim(0, (st_day+num_day)*24)
	
	optLLSH2=ggplot() +
			geom_line(data=out_plot1i1LSH,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip)) +
			geom_line(data=out_plot1i2LSH,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip),color="red",linetype="dashed") +
			theme(legend.position = "none") +
			labs(x ="time (hr)", y = "%Damaged Cell_short")#+xlim(0, (st_day+num_day)*24)
	
	optp1LSH3<-ggplot()+
			geom_line(data=out_plot1i1LL, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "solid",col="blue")+
			geom_line(data=out_plot1i2LL, aes(x=time, y=log10(FreeInfluenza),group=ip), size=0.85, alpha=0.8,linetype = "dashed",col="green")+
			xlab(paste0("time (hr)"))+
			ylab("Virus log_long")#+ggtitle(paste0(caseMS,": \n Black= without drug, Green= with drug"))#+xlim(0, (st_day+num_day)*24)
	
	optLLSH4=ggplot() +
			geom_line(data=out_plot1i1LL,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip),col="blue") +
			geom_line(data=out_plot1i2LL,aes(x=time, y=(4E8-ResistantEpithelium-EpitheliumLung)/4E8*100,group=ip),color="green",linetype="dashed") +
			
			theme(legend.position = "none") +
			labs(x ="time (hr)", y = "%Damaged Cell_long")#+xlim(0, (st_day+num_day)*24)
	plotlist6=list()
	
	plotlist6[[1]] = optp1LSH1
	plotlist6[[2]] = optLLSH2
	plotlist6[[3]] = optp1LSH3
	plotlist6[[4]] = optLLSH4
	pall6 <- grid.arrange(grobs=plotlist6,ncol=2)
	ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,"_LOSH4.png"), pall6, width=8, height=8)

	k1<-ggplot()+
			geom_line(data=out_plot1i1L, aes(x=time, y=EpiKilledbyTc,group=ip), size=0.85, alpha=0.8,linetype = "solid")+
			geom_line(data=out_plot1i2L, aes(x=time, y=EpiKilledbyTc,group=ip), size=0.85, alpha=0.8,linetype = "solid",col="red")+
			xlab(paste0("time (hr)"))+
			ylab("EpiKilledbyTc")#+ggtitle(paste0(caseMS,": \n Black= without drug, Green= with drug"))#+xlim(0, (st_day+num_day)*24)
	
	k2=ggplot() +
			geom_line(data=out_plot1i1L,aes(x=time, y=EpiKilledbyIL6,group=ip)) +
			geom_line(data=out_plot1i2L,aes(x=time, y=EpiKilledbyIL6,group=ip),col="red") +

			theme(legend.position = "none") +
			labs(x ="time (hr)", y = "EpiKilledbyIL6")#+xlim(0, (st_day+num_day)*24)
	
	k3<-ggplot()+
			geom_line(data=out_plot1i1L, aes(x=time, y=EpiKilledbyIFN,group=ip), size=0.85, alpha=0.8,linetype = "solid")+
			geom_line(data=out_plot1i2L, aes(x=time, y=EpiKilledbyIFN,group=ip), size=0.85, alpha=0.8,linetype = "solid",col="red")+
			xlab(paste0("time (hr)"))+ylab("EpiKilledbyIFN")
	plotlist6=list()
	
	plotlist6[[1]] = k1
	plotlist6[[2]] = k2
	plotlist6[[3]] = k3
	pall6 <- grid.arrange(grobs=plotlist6,ncol=1)
	ggsave(paste0(outdir2,byplot,EC,"_","PK",caseMS,"_kill.png"), pall6, width=8, height=8)
	if (i1==1) {Damageper2placebo=Damageper2}
	if (i1==2) {Damageper2RMD=Damageper2}
	if (caseMS=="severReady") {RC_C=Sever_Control;RC_D=Sever_RMD;RC_C_step=Sever_Control_step;RC_D_step=Sever_RMD_step;RC_C_meanstep=Sever_Control_meanstep;RC_D_meanstep=Sever_RMD_meanstep;}
	if (caseMS=="mildReady") {RC_C=Mild_Control;RC_D=Mild_RMD;RC_C_step=Mild_Control_step;RC_D_step=Mild_RMD_step;RC_C_meanstep=Mild_Control_meanstep;RC_D_meanstep=Mild_RMD_meanstep;}
	RC_C[,2:4]=RC_C[,2:4]*100 #used for error bar
	RC_D[,2:4]=RC_D[,2:4]*100 
	RC_C_meanstep[,2:4]=RC_C_meanstep[,2:4]*100 #used for dash line and fiting 
	RC_D_meanstep[,2:4]=RC_D_meanstep[,2:4]*100
	RC_C_step[,2:4]=RC_C_step[,2:4]*100 #used for step plot
	RC_D_step[,2:4]=RC_D_step[,2:4]*100
	Damageper2placebo=Damageper2placebo[(Damageper2placebo$time/24)<=max(RC_C$day),]
	Damageper2RMD=Damageper2RMD[(Damageper2RMD$time/24)<=max(RC_D$day),]
	Rec<-ggplot()+
			geom_line(data=Damageper2placebo, aes(x=time/24, y=Recovered ), size=1.5, alpha=1,linetype = "solid",col="red")+
			geom_line(data=Damageper2RMD, aes(x=time/24, y=Recovered ), size=1.5, alpha=1,linetype = "solid",col="deepskyblue")+
		geom_errorbar(data=RC_D, aes(x=day, ymin=(min),ymax=(max)),color="blue",fill="blue",width=.3,
				position=position_dodge(.9),color="blue",fill="blue")+
		geom_errorbar(data=RC_C, aes(x=day, ymin=(min),ymax=(max)),color="indianred3",fill="indianred3",width=.3,
				position=position_dodge(.9),color="indianred3",fill="indianred3")+
			geom_line(data=RC_D, aes(x=day, y=opt), size=0.85, alpha=0.8,linetype = "dashed",col="blue")+
			geom_line(data=RC_C, aes(x=day, y=opt), size=0.85, alpha=0.8,linetype = "dashed",col="indianred3")+
			xlab(paste0("Time(day)"))+ylab("Recovered%")+ylim(0,105)+ggtitle(paste0(caseMS,"_",iT,"_seed",seedn))+
	theme(axis.text=element_text(size=20),axis.title=element_text(size=labsiz,face="bold"))
	ggsave(paste0(outdir2,byplot,EC,"_","Recovery",caseMS,seedn,"ooo","_kill.png"), Rec, width=8, height=8)
	Rec_JustCtrl<-ggplot()+
		geom_line(data=Damageper2placebo, aes(x=time/24, y=Recovered ), size=.6, alpha=1,linetype = "solid",col="black")+
		geom_errorbar(data=RC_C, aes(x=day, ymin=(min),ymax=(max)),width=.8,
		position=position_dodge(2),color="black",fill="black")+
		geom_point(data=RC_C, aes(x=day, y=opt), size=1.5, alpha=1,color="black")+
		xlab(paste0("Time(day)"))+ylab("Recovered%")+ylim(0,100)+xlim(0,27)+ scale_x_continuous(breaks = seq(0, 27, by = 9))+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=labsiz,face="bold"))+ theme_bw()
	ggsave(paste0(outdir2,byplot,EC,"_","Recovery_Rec_JustCtrl",caseMS,"_kill.png"), Rec_JustCtrl, width=3.38, height=3.38)
	Recnoctr<-ggplot()+
		geom_line(data=Damageper2RMD, aes(x=time/24, y=Recovered ), size=.6, alpha=1,linetype = "solid",col="black")+
		geom_errorbar(data=RC_D, aes(x=day, ymin=(min),ymax=(max)),color="black",fill="black",width=.8,
		position=position_dodge(2),color="black",fill="black")+
		geom_point(data=RC_D, aes(x=day, y=opt), size=1.5, alpha=1,color="black")+
		
		xlab(paste0("Time(day)"))+ylab("Recovered%")+ylim(0,100)+xlim(0,27)+ scale_x_continuous(breaks = seq(0, 27, by = 9))+
			theme(axis.text=element_text(size=20),axis.title=element_text(size=labsiz,face="bold"))+ theme_bw()
	ggsave(paste0(outdir2,byplot,EC,"_","Recovery",caseMS,"_killRecnoctr.png"), Recnoctr, width=3.38, height=3.38)
	
	}				