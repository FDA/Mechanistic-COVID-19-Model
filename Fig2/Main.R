# Codes are written by: Mohammadreza Samieegohar
# Figurpe 2
#---------------------------------------------------
proc_start<-proc.time()
options(warn=-1)
# load libraries -----------------------------------
library(optparse)
library(deSolve)
library(ggplot2)
library(gridExtra)
library(grid)
# Specify command line arguments-------------------
isWindows<-Sys.info()[["sysname"]]=="Windows"
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--severity"), default="mild",type="character", help="severity group,options:'mild' or 'severe'")
parser<-add_option(parser, c("-c", "--case"), default="Optimum", help="options: 'Population' or 'Optimum'")
parser<-add_option(parser, c("-m", "--ModelPlatform"), default="R",type="character", help="if 'R' used, then the R version of the model will be used, which is compatible with Windows and Linux, otherwise a *so format is used that is compatible with Linux.")
parser<-add_option(parser, c("-p", "--popNum"), default=2000,type="integer", help="Virtual population size")
args<-parse_args(parser)
args$drug=args$severity
caserun=args$case
ModelPlatform=args$ModelPlatform
popNum=as.numeric(args$popNum)
if (args$drug=="severe") {caseMS="severReady"}
if (args$drug=="mild") {caseMS="mildReady"}
iT=as.numeric(args$ERR)
print(sprintf("Severity is: -----%s-----",args$drug))
print(sprintf("Run is:    -----%s-----",args$case))
# Initial states and parameters and model names ---------------------------
mymodel<-list()
modelname<-"delaymymod"
modeldir<-"Model/"
nnnfile<-paste0(modeldir,"nnn.txt")
source(paste0(modeldir,"delaypars.R"))
source(paste0(modeldir,"delaystates.R"))
timepoints=seq(0,577,1)
if (caserun=="Optimum") {popNum=0}
# Loading clinical data----------------------------------------------------
VIRUS=read.csv(paste0("Data_Covid/Viral_Titer/",caseMS,".csv"))
if (caseMS=="mildReady") {
	VIRUS=VIRUS[-1:-2,]
	VIRUS=VIRUS[VIRUS$Rx<23,]
	addday=7
}
if (caseMS=="severReady") {
	VIRUS=VIRUS[VIRUS$Rx<23,]
	addday=5.5
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
IFNAalpha0<-read.csv(paste0("Data_Covid/IFNAalpha/","IFNA_Fig1A_withmaxmin",".csv"),sep = ",") 
IFA_sever=read.csv("Data_Covid/IFNAalpha/IFA_sever.csv")

if (caseMS=="severReady") {	
	NSIgM_data<-read.csv(paste0("Data_Covid/IgMIgG/","SIgM_Ready",".csv"),sep = ",") 
	NSIgM_data=NSIgM_data[NSIgM_data$day<23,] 		
	NSIgM_data$Time=(NSIgM_data$day+addday)*24    
	NSIgG_data<-read.csv(paste0("Data_Covid/IgMIgG/","SIgG_Ready",".csv"),sep = ",") 
	NSIgG_data=NSIgG_data[NSIgG_data$day<23,] 		
	NSIgG_data$Time=(NSIgG_data$day+addday)*24    	
}
if (caseMS=="mildReady") {	
	NSIgM_data<-read.csv(paste0("Data_Covid/IgMIgG/","NSIgM_Ready",".csv"),sep = ",") 
	NSIgM_data=NSIgM_data[NSIgM_data$day<23,] 		 
	NSIgM_data$Time=(NSIgM_data$day+addday)*24    
	NSIgG_data<-read.csv(paste0("Data_Covid/IgMIgG/","NSIgG_Ready",".csv"),sep = ",") 
	NSIgG_data=NSIgG_data[NSIgG_data$day<23,] 		
	NSIgG_data$Time=(NSIgG_data$day+addday)*24    
}
Bioout_plot_paperall=datares1=datares2=datares3=datares4=c()
# Virtual Population Calculator Loop-----------------------------------------------------
for (ip in 1:(popNum+1)) {
	#for (ip in 2) {
	if(ip == 250 | ip == 500 |ip == 1000 |ip == 1000 ) {print(ip)}
	pars["timeout"]=30
	pars["unnamed"]=1
	pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]		
# Input parameters ----------------------------------------------------------------------
	if (caseMS=="mildReady") {
		Par_population0=read.csv("Inputs_parameters/Population_parameters_Mild.csv")
		lastcoluse=which(colnames(Par_population0)=="starttime")
		pars=Par_population0[ip,13:lastcoluse]
		initstates0=Par_population0[ip,67]
#		initstates0=read.csv("Inputs_parameters/Mild_virusdiss.csv")
		#Optimal parameters mild----
		if (ip==1) {
			Par_population0=read.csv("Inputs_parameters/Optimal_pars_mild.csv")
			pars=Par_population0[,"x"]
			names(pars)=Par_population0[,"X"]
			initstates00=read.csv("Inputs_parameters/Optimal_initstates_mild.csv")	
			initstates0=(initstates00[2,2])
		}		
	}	
	if (caseMS=="severReady") {
		Par_population0=read.csv("Inputs_parameters/Population_parameters_Severe.csv")
		lastcoluse=which(colnames(Par_population0)=="starttime")
		pars=Par_population0[ip,13:lastcoluse]
		initstates0=Par_population0[ip,67]
#		initstates0=read.csv("Inputs_parameters/Severe_virusdiss.csv")
		#Optimal parameters severe----
		if (ip==1) {
			Par_population0=read.csv("Inputs_parameters/Optimal_pars_severe.csv")
			pars=Par_population0[,"x"]
			names(pars)=Par_population0[,"X"]
			initstates00=read.csv("Inputs_parameters/Optimal_initstates_severe.csv")	
			initstates0=(initstates00[2,2])
		}
	}	
#	states["FreeVirus"]=initstates0[ip,2]
	states["FreeVirus"]=initstates0
	
	pars["timeout"]=30
	pars["unnamed"]=1
#  Model run ----------------------------------------------------------------------------
	if (ModelPlatform=="R") {
		source("Model/delaymodelfun.R")		
		pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]			
		try({out_plot_paper=dede(y=states,times=timepoints,func=modelfun,method="lsoda",parms=unlist(pars),atol = 1e-6,rtol=1e-3)})
	
	}
	else{	
		isWindows<-Sys.info()[["sysname"]]=="Windows"		
		if (isWindows) {print("for running in Windows ModelPlatform should be=  'R'");WindowsError}
		extension<-ifelse(isWindows, ".dll", ".so")
		dyn.load(paste0(modeldir,modelname,extension))		
		pars["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]	
		pars=unlist(pars)
		try({out_plot_paper <- dede(states, timepoints, "derivs", pars, dllname=modelname,initfunc="initmod",rtol=1e-3,atol=1e-6,nout=3)});		
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
	
	if (ip<=500) 			{datares1=data.frame(rbind(datares1,Bioout_plot_paper))}
	if (ip>500 & ip<=1000)  {datares2=data.frame(rbind(datares2,Bioout_plot_paper))}
	if (ip>1000 & ip<=1500) {datares3=data.frame(rbind(datares3,Bioout_plot_paper))}
	if (ip>1500) 			{datares4=data.frame(rbind(datares4,Bioout_plot_paper))}	
	Bioout_plot_paperall=rbind(datares1,datares2,datares3,datares4)
	
	
#	Bioout_plot_paperall=rbind(Bioout_plot_paperall,Bioout_plot_paper)
	Bioout_plot_paperall=data.frame(Bioout_plot_paperall)
	if (ip==1 | caserun=="Optimum") {		
		out_plot1=Bioout_plot_paper		
	}
}

# Calculating the confidence interval for uncertainly bands ---------------------------------
BIL6q=BIgMq=BIgGq=BTnaiveq=BIFNA1q=BFreeVirusq=BREcell=BEpitheliumLungq=c()
Bioout_plot_paperall$FreeVirus=log10(Bioout_plot_paperall$FreeVirus)
stepQ=Bioout_plot_paperall[2,"time"]-Bioout_plot_paperall[1,"time"]
for (ti1 in seq(0,max(Bioout_plot_paperall[,"time"]),stepQ)) { 	
	tout_plot1=Bioout_plot_paperall[Bioout_plot_paperall[,1]==ti1,]
	tq=(lapply(tout_plot1,quantile,probs=c(0.05,.95),na.rm=TRUE))	
	BIL6q=rbind(BIL6q,c(ti1,tq$IL6))
	BFreeVirusq=rbind(BFreeVirusq,c(ti1,tq$FreeVirus))
	BIFNA1q=rbind(BIFNA1q,c(ti1,tq$IFNA1))
	BTnaiveq=rbind(BTnaiveq,c(ti1,tq$Tnaive))
	BIgGq=rbind(BIgGq,c(ti1,tq$IgG))
	BIgMq=rbind(BIgMq,c(ti1,tq$IgM))
	BEpitheliumLungq=rbind(BEpitheliumLungq,c(ti1,tq$EpitheliumLung))
}

colnames(BEpitheliumLungq)=c("time","qmin","qmax")
colnames(BIL6q)=c("time","qmin","qmax")
colnames(BFreeVirusq)=c("time","qmin","qmax")
colnames(BIFNA1q)=c("time","qmin","qmax")
colnames(BTnaiveq)=c("time","qmin","qmax")
colnames(BIgGq)=c("time","qmin","qmax")
colnames(BIgMq)=c("time","qmin","qmax")

BEpitheliumLungq1=data.frame(BEpitheliumLungq)
BIL6q1=data.frame(BIL6q)
BFreeVirusq1=data.frame(BFreeVirusq)
BIFNA1q1=data.frame(BIFNA1q)
BTnaiveq1=data.frame(BTnaiveq)
BIgGq1=data.frame(BIgGq)
BIgMq1=data.frame(BIgMq)

#Plotting -------------------------------------------------------
outdir2="Fig2/"
if (caserun=="Optimum") {source("Plot_Optimum.R")} 
if (caserun=="Population") {source("Plot_Population.R")}
