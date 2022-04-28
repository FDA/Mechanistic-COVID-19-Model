# Codes are written by: Mohammadreza Samieegohar
# Figure 5. Calibrating (a,b) and validating (c,d) the model for the primary endpoint (time to recovery) used in the remdesivir trial

if (args$drug=="severe") {caseMS="severReady"}
if (args$drug=="mild") {caseMS="mildReady"}
# Initial states and parameters and model names ---------------------------
mymodel<-list()
modelname<-"delaymymod"
modeldir<-"Model/"
source(paste0(modeldir,"delaypars.R"))
source(paste0(modeldir,"delaystates.R"))
timepoints=seq(0,max(c(((10)*24+30*24),650)),1)
Bioout_plot_paperall=cross5per1=datares1=datares2=datares3=datares4=c()
# Virtual Population Calculator Loop-----------------------------------------------------
for (ip in 1:(popNum)) {
	if(ip == 250 | ip == 500 |ip == 1000 |ip == 1500 ) {print(ip)}
	pars["timeout"]=30
	pars["unnamed"]=1
	pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]		
# Input parameters ----------------------------------------------------------------------
	if (caseMS=="mildReady") {
		Par_population0=read.csv("Inputs_parameters/Population_parameters_Mild.csv")
		pars=Par_population0[ip,13:65]
		initstates0=Par_population0[ip,67]	
	}	
	if (caseMS=="severReady") {
		Par_population0=read.csv("Inputs_parameters/Population_parameters_Severe.csv")
		pars=Par_population0[ip,13:65]
		initstates0=Par_population0[ip,67]
	}	
	states["FreeVirus"]=initstates0	
	pars["timeout"]=30
	pars["unnamed"]=1	
#  Model run ----------------------------------------------------------------------------
source("eventRMD.R")	
if (ModelPlatform=="R") {
		pars["timeout"]=60	
		source("Model/delaymodelfun.R")		
		pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]			
		try({out_plot_paper=dede(y=states,times=timepoints,func=modelfun,method="lsodes",parms=unlist(pars),nout=3,events=list(data=eventdata))})
	}
	else{		
		isWindows<-Sys.info()[["sysname"]]=="Windows"		
		if (isWindows) {print("for running in Windows ModelPlatform should be=  'R'");WindowsError}
		extension<-".so"
		dyn.load(paste0(modeldir,modelname,extension))		
		pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]	
		pars=unlist(pars)
		try({out_plot_paper <- dede(states, timepoints, "derivs", pars, dllname=modelname,initfunc="initmod",rtol=1e-3,atol=1e-6,nout=3,events=list(data=eventdata))});		
	}
	out_plot_paper=data.frame(out_plot_paper)
	out_plot_paper$ip=ip	
	out_plot_paper$initVirus=(out_plot_paper$FreeVirus)[out_plot_paper$time==20]
	out_plot_paper$maxIL6=max(out_plot_paper$IL6)
	out_plot_paper$damage=max(4E8-out_plot_paper$ResistantEpithelium-out_plot_paper$EpitheliumLung)/4E8*100
	out_plot_paper$damage0=(4E8-out_plot_paper$ResistantEpithelium-out_plot_paper$EpitheliumLung)/4E8*100
	out_plot_paper$maxFreeVirus=max(out_plot_paper$FreeVirus)
	Il6KP=out_plot_paper$EpiKilledbyIL6/(out_plot_paper$EpiKilledbyIFN+out_plot_paper$EpiKilledbyTc+out_plot_paper$EpiKilledbyVirus+out_plot_paper$EpiKilledbyIL6)
	VKp=out_plot_paper$EpiKilledbyVirus/(out_plot_paper$EpiKilledbyIFN+out_plot_paper$EpiKilledbyTc+out_plot_paper$EpiKilledbyVirus+out_plot_paper$EpiKilledbyIL6)
	out_plot_paper$il6killedpercent=max(Il6KP[Il6KP!="Inf" & Il6KP!="NaN"])
	out_plot_paper$Vkilledpercent=max(VKp[VKp!="Inf" & VKp!="NaN"])
	Bioout_plot_paper=out_plot_paper
	if (ip==1 | caserun=="Optimum") {		
		out_plot1=Bioout_plot_paper		
	}
	eventtimes=Par_population0[ip,"minevent"]
	daybase=min(eventtimes)
	reme1=daybase%%2
	if (reme1==1){daybase=daybase+1}else{daybase=daybase}
	out_plot_paper$eventtimes=daybase
	out_plot_paper[,"time"]=out_plot_paper[,"time"]-daybase #one day befor first injection
	out_plot_paper=out_plot_paper[out_plot_paper[,"time"]>=0,]
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
		#----------- check what is the damage at different days -----------------------
		damageD0=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*0][1]
		damageD1=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*1][1]
		damageD10=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*10][1]
		damageD11=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*11][1]
		damageD12=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*12][1]
		damageD15=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*15][1]
		damageD28=out_plot1_paper_ip$Damageper[out_plot1_paper_ip$time>=24*28][1]
		#--------------------------------------------------------------
ipun=ip
NSvec=c("ipr","ip","Recoveredt","minevent","damageD0","damageD1","damageD10","damageD11","damageD12","damageD15","damageD28",names(pars),names(states))	
Svec0= c(ip,ipun,cross5pert0,min(eventtimes),damageD0,damageD1,damageD10,damageD11,damageD12,damageD15,damageD28,pars,states)
Svec=unlist(Svec0)
names(Svec)=NSvec
if (ip<=500) 			{datares1=data.frame(rbind(datares1,Svec))}
if (ip>500 & ip<=1000)  {datares2=data.frame(rbind(datares2,Svec))}
if (ip>1000 & ip<=1500) {datares3=data.frame(rbind(datares3,Svec))}
if (ip>1500) 			{datares4=data.frame(rbind(datares4,Svec))}	

	cross5per1=rbind(datares1,datares2,datares3,datares4)
	colnames(cross5per1)=NSvec
	cross5per1pre_RMD=data.frame(cross5per1)

}

