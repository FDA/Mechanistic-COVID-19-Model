# All codes are written by: Mohammadreza Samieegohar
# Figure 4b
#---------------------------------------------------
proc_start<-proc.time()
options(warn=-1)
library(optparse)
library(ggplot2)
library(gridExtra)
library(deSolve)
isWindows<-Sys.info()[["sysname"]]=="Windows"
parser<-OptionParser()
parser<-add_option(parser, c("-c", "--case"), default="Optimum",type="character", help="Population or Optimum")
parser<-add_option(parser, c("-p", "--popNum"), default="2000", help="population number")
args<-parse_args(parser)
popNum=args$popNum
print(sessionInfo())
case=args$case
github="Fig4b/"
caseN="HAE_Ave"
parsall=c("1apars")
ftn=0
for (ft in 1) {
	parnums=parsall[1]
	tmp<-read.table(paste0("Inputs_parameters/hergmod_states.txt"), col.names=c("param","value"))
	states<-setNames(tmp$value, tmp$param)
	tmp<-read.table(paste0("Inputs_parameters/hergmod_pars.txt"), col.names=c("param","value"))
	pars<-setNames(tmp$value, tmp$param)
#	pp<-read.table("Inputs_parameters/hergmod_drug_param_bounds.txt",header=T,as.is=T)
#	pnames<-pp$Parameter
#	high_bounds<-pp$High
#	low_bounds<-pp$Low
	pmax<-10000
#	encode_pars<-function(pars) pmax*log10(pars/low_bounds)/log10(high_bounds/low_bounds)
#	decode_pars<-function(ind) low_bounds*(high_bounds/low_bounds)^(ind/pmax)
	Verp6=cbind(c(8*3600,24*3600,48*3600),c(1.21,.5,.61)/10^6/.665)       
	Calu3=cbind(c(8*3600,24*3600,48*3600),c(2.87,2.17,2.00)/10^6/2.7)     
	HAE1=cbind(c(8*3600,24*3600,48*3600),c(18.3,15.3,2.45))               
	HAE2=cbind(c(8*3600,24*3600,48*3600),c(6.58,5.78,0.73))               
	HAE_Ave0=cbind(c(8*3600,24*3600,48*3600),c(12.44,10.54,1.59))         
	HAE_Ave=cbind(c(8*3600,24*3600,48*3600),c(12.44,10.54,1.59)/1000/9119)
	if (caseN=="Verp6")   {expData=Verp6}
	if (caseN=="Calu3")   {expData=Calu3}
	if (caseN=="HAE_Ave") {expData=HAE_Ave}
	incubation_time=c(48,72,48)
	names(incubation_time)=c("VeroE6","Calu3","HAE_Ave")
	colnames(expData)=c("time","NTP")
	ftime=seq(0,50*3600,50)
	caseN0=caseN
	if (caseN0=="Verp6") {C0=1E5; V_percell=0.66E-12}  
	if (caseN0=="Calu3") {C0=2E5; V_percell=2.7E-12} 
	if (caseN0=="HAE_Ave") {C0=1E6; V_percell=9.119E-12} 
	states["cell_number"]=C0
	hergmod1<-function(t, state, pars) {
		with(as.list(c(state, pars)),{
					Kcell=0
					F_intracellular=1   
					medium_volume=10^-3 					
					dC_intracellularRDV_total=CL_percell*(C_extracellularRDV_total*F_extracellular)/V_percell - CL_percell*(C_intracellularRDV_total*F_intracellular)/V_percell - Km*(C_intracellularRDV_total*F_intracellular)-Knew*C_intracellularRDV_total									
					dC_extracellularRDV_total=CL_percell*cell_number*(C_intracellularRDV_total*F_intracellular)/medium_volume - CL_percell*cell_number*(C_extracellularRDV_total*F_extracellular)/medium_volume
					dcell_number = Kcell*cell_number
					dTP = Km*(C_intracellularRDV_total*F_intracellular) - Kc*TP					
					list(c(dC_intracellularRDV_total,dC_extracellularRDV_total,dcell_number,dTP))
					
				}) 
	}	
	objfun<-function(ind){		
		parsin=(ind)		
		names(parsin)=names(pars)
		try({out_pre <- ode(states, ftime,func=hergmod1, parms=parsin,
		method="lsoda",atol = 1e-6, rtol = 1e-6)});	
		deptime<-expData[,"time"]
		deptimeround=round(deptime)
		idxPeaktime<-out_pre[,"time"]%in%deptimeround
		idxPeaktimeT=which(idxPeaktime==TRUE)
		yPred<-out_pre[idxPeaktimeT,] 
		fval1=sum(((yPred[,"TP"]-expData[,"NTP"])/max(expData[,"NTP"]))^2)
		print(fval1)		
		print("--------")			
		tmaxx=out_pre[,"time"][out_pre[,"TP"]==max(out_pre[,"TP"])]		
		if (tmaxx<6*3600 | is.na(tmaxx[1])) {maxPE=10}else{maxPE=0}			
		fval=fval1+maxPE
		(yPred[,"TP"]-expData[,"NTP"])
		
	}		
#	initpar=pp$Initial
#	names(initpar)=pp$Parameter
	parsin=pars	
	ttt<-read.csv(paste0("Inputs_parameters/1apars.csv"),header=T)
	Allpar=ttt
	Allpar=Allpar[,-1]
	use_spike="yes"
	if (use_spike=="yes") {		
		ttto<-read.table(paste0("Inputs_parameters/pars_HAE_Ave.txt"),header=F,as.is=T)	
		initpar_opt=(ttto[,2])
		names(initpar_opt)=ttto[,1]
	}
	out_plot_all=c()

	parsin=initpar_opt
	try({out_plotpp <- ode(states, ftime,func=hergmod1,parms=parsin,
    method="lsoda",atol = 1e-6, rtol = 1e-6)})	
	out_plotpp=data.frame(out_plotpp)
	out_plot1_opt0=out_plotpp
	if (case=="Optimum"){PoPcase=10}else{PoPcase=popNum}	
	for (ip1 in c(seq(1,PoPcase,1))) {	
		
		out_plot=c()
		parsin=Allpar[ip1,]		
		try({out_plot <- ode(states, ftime,func=hergmod1,parms=parsin,
		method="lsoda",atol = 1e-6, rtol = 1e-6)})		
		out_plot=data.frame(out_plot)
		out_plot$idn=rep(ip1,length(out_plot$TP))
		out_plot_all=rbind(out_plot_all,out_plot)		
	}	
	C_intracellularRDV_total=C_extracellularRDV_total=cell_number=TP=c()
	for (ti1 in ftime) {
		tout_plot1=out_plot_all[out_plot_all[,1]==ti1,]
		tq=(lapply(tout_plot1,quantile,probs=c(0.025,0.975),na.rm = TRUE))				
		C_intracellularRDV_total=rbind(C_intracellularRDV_total,c(ti1,tq$C_intracellularRDV_total))
		C_extracellularRDV_total=rbind(C_extracellularRDV_total,c(ti1,tq$C_extracellularRDV_total))
		cell_number=rbind(cell_number,c(ti1,tq$cell_number))
		TP=rbind(TP,c(ti1,tq$TP))		
	}
	
	colnames(C_intracellularRDV_total)=c("time","qmin","qmax")
	colnames(C_extracellularRDV_total)=c("time","qmin","qmax")
	colnames(cell_number)=c("time","qmin","qmax")
	colnames(TP)=c("time","qmin","qmax")	
	C_intracellularRDV_total=data.frame(C_intracellularRDV_total)
	C_extracellularRDV_total=data.frame(C_extracellularRDV_total)
	cell_number=data.frame(cell_number)
	TP=data.frame(TP)
	
	plotlist = list()
	expData=data.frame(expData)
	HAE1=data.frame(HAE1)
	HAE2=data.frame(HAE2)
	expData$mind=HAE2$X2/1000/9119
	expData$maxd=HAE1$X2/1000/9119	
	out_plot1_opt0$TP=out_plot1_opt0$TP/10^-6
	out_plot1_opt0$time=out_plot1_opt0$time/3600	
	TP$qmin=TP$qmin/10^-6
	TP$qmax=TP$qmax/10^-6
	TP$time=TP$time/3600	
	expData$time=expData$time/3600
	expData$NTP=expData$NTP/10^-6
	expData$mind=expData$mind/10^-6
	expData$maxd=expData$maxd/10^-6
	if (case!="Optimum"){		
		p4<-ggplot( )
		p4<-p4+geom_point(data=expData, aes(x=time, y=NTP), size=2, alpha=1, color="black")+
		geom_errorbar(data=expData,aes(x=time,ymin=mind, ymax=maxd),width=2,
		position=position_dodge(0.05))		
		p4<-p4+geom_line(data=out_plot1_opt0, aes(x=time, y=TP), size=0.85, alpha=.85,linetype = "solid",color = "black")
		p4<-p4+geom_ribbon(data=TP,aes(x=time,ymin=qmin,ymax=qmax),alpha=0.3)
		p4<-p4+ylab(paste0("TP"))
		p4<-p4+ylab(paste0("TP (uM)"))+xlab(paste0("Time post infection (hours)"))+ theme_bw()
		ggsave(paste0(github,"/","Population.png"), p4, width=3.38, height=3.38)
		print("************************************************************")
		print("--------------*******************************---------------")
		print("------------------------*********---------------------------")
		print("---------------------------***------------------------------")
		print(" The Population version of figure 4B is successfully plotted")
		
	}else{
		p4<-ggplot( )
		p4<-p4+geom_point(data=expData, aes(x=time, y=NTP), size=2, alpha=1, color="black")+
		geom_errorbar(data=expData,aes(x=time,ymin=mind, ymax=maxd),width=2,
		position=position_dodge(0.05))
		p4<-p4+geom_line(data=out_plot1_opt0, aes(x=time, y=TP), size=0.85, alpha=.85,linetype = "solid",color = "black")	
		p4<-p4+ylab(paste0("TP"))
		p4<-p4+ylab(paste0("TP (uM)"))+xlab(paste0("Time post infection (hours)"))+ theme_bw()
		ggsave(paste0(github,"/","Optimum.png"), p4, width=3.38, height=3.38)
		print("************************************************************")
		print("--------------*******************************---------------")
		print("------------------------*********---------------------------")
		print("---------------------------***------------------------------")
		print(" The Optimum version of figure 4B is successfully plotted ")
	}
}
