# All codes are written by: Mohammadreza Samieegohar
# Figure 4a
#---------------------------------------------------
proc_start<-proc.time()
options(warn=-1)
library(optparse)
library(ggplot2)
library(gridExtra)
library(deSolve)
isWindows<-Sys.info()[["sysname"]]=="Windows"
parser<-OptionParser()
parser<-add_option(parser, c("-p", "--populationSize"), default="2000", help="population size")
args<-parse_args(parser)
populationSize=args$populationSize
print(sessionInfo())
data_RMD=read.csv("data/all.csv")
data_RMD=data_RMD[data_RMD["time"]!=2,] 
data_RMD=data_RMD[data_RMD["time"]<4.11,]
data_RMD["time"]=data_RMD["time"]*3600
concvec=unique(data_RMD$dose)
tmp<-read.table(paste0("input_parameters/hergmod_states.txt"), col.names=c("param","value"))
states<-setNames(tmp$value, tmp$param)
hergmod1<-function(t, state, pars) {
	with(as.list(c(state, pars)),{				
				inf_stop=1
				if (t>(3600*2)) {inf_stop=0}										
				Rinter=-K14*CR+K41*CRP
				Rout=-KC*CR						
				dCR=Rinter+Rout+KCon/(3600*2)*inf_stop
				dCRP=-Rinter			
				dKCon=0
				list(c(dCR,dCRP,dKCon))				
			}) 
}
alldrugsweeps<-list()
drugsweeps_plot=c()
opt="yes"
drugsweeps_plotop=drugsweepsop1=c()
ttt<-read.csv(paste0("input_parameters/9apars.csv"),header=T)
Allpar=ttt
Allpar=Allpar[,-1]
for (ip1 in c(seq(1,populationSize,1))) {
	parsin=Allpar[ip1,]
	for(concj in concvec){		
		expData=data_RMD[data_RMD$dose==concj,][,c("time","M")]
		colnames(expData)=c("time","NTP")
		ftime=seq(0,max(expData["time"]),100)
		states["KCon"]=concj*10^6
		try({drugsweeps <- ode(states, ftime,func=hergmod1, parms=parsin,
		method="lsoda",atol = 1e-6, rtol = 1e-6)});
		if(length(drugsweeps)==0)
		return(c())
		drugsweeps=cbind(drugsweeps,CRv=drugsweeps[,"CR"]/as.numeric(parsin["V1"]))
		alldrugsweeps[[concj]]<-drugsweeps
		drugsweeps1=cbind(drugsweeps,concj,ip1)
		drugsweeps_plot=rbind(drugsweeps_plot,drugsweeps1)		
		use_spike="yes"
		if (use_spike=="yes") {
			ttt<-read.table(paste0("input_parameters/pars_Verp6.txt"),header=F,as.is=T)
			initparop=(ttt[,2]) 
			names(initparop)=ttt[,1]			
		}
		if (opt=="yes") {	
			try({drugsweepsop <- ode(states, ftime,func=hergmod1, parms=initparop,							method="lsoda",atol = 1e-6, rtol = 1e-6)});
			if(length(drugsweepsop)==0)
				return(c())
			drugsweepsop=cbind(drugsweepsop,CRv=drugsweepsop[,"CR"]/initparop["V1"])		
			drugsweepsop1=cbind(drugsweepsop,concj)
			drugsweeps_plotop=rbind(drugsweeps_plotop,drugsweepsop1)			
		}
	}
	opt="no"
}
drugsweeps_plot=data.frame(drugsweeps_plot)
out_plot_all=drugsweeps_plot
C_intracellularRDV_total=C_extracellularRDV_total=cell_number=TP=c()
for (ci1 in concvec) {
	ftimemax=max(out_plot_all[out_plot_all[,"concj"]==ci1,"time"])
	ftime1=seq(0,ftimemax,100)
	for (ti1 in ftime1) {
		tout_plot1=out_plot_all[out_plot_all[,"time"]==ti1 & out_plot_all[,"concj"]==ci1,]
		tq=(lapply(tout_plot1,quantile,probs=c(0.025,0.975),na.rm = TRUE))
		
		TP=rbind(TP,c(ti1,ci1,tq$CRv))
	}
}

colnames(TP)=c("time","conc","LL","HL")
data_RMD=data.frame((data_RMD))
data_RMD$dose=as.numeric(data_RMD$dose)

TP=as.data.frame(TP)
drugsweeps_plotop=as.data.frame(drugsweeps_plotop)
p_paper=ggplot(data=TP,group=conc,color=as.factor(conc)) +
		geom_ribbon(data=TP,aes(x=time/3600,ymin=LL, ymax=HL,fill=as.factor(conc)),alpha=.3)+
		geom_errorbar(data=data_RMD,aes(x=time/3600,ymin=L, ymax=H),width=.1,position=position_dodge(.1))+
		geom_line(data=drugsweeps_plotop,aes(x=time/3600, y=CRv,group=concj),linetype = "solid") +
		geom_point(data=data_RMD,aes(x=time/3600, y=M, group=dose,shape=as.factor(dose)),size=1.5)+	
		scale_color_manual(values=c("black", "black", "black", "black", "black","black"))+
		scale_color_manual(values=c("black", "black", "black", "black", "black","black"))+
		scale_fill_manual(values=c("gray", "gray", "gray", "gray", "gray","gray"))+
		scale_shape_manual(values=c(15,16,17,0,1,2))+ 
		scale_colour_discrete("dose(mg)")+
		guides(color = FALSE, size = FALSE,fill = FALSE, shape=guide_legend(title="Dose (mg)"))+
		labs(x ="Time (hour)", y = "Remdisivir in plasma (ng/ml)")+ theme_bw()+ theme(legend.position = c(0.8, 0.6))
ggsave(sprintf("results/Figure4A.png"),p_paper,width=3.38, height=3.38)
print("************************************************************")
print("--------------*******************************---------------")
print("------------------------*********---------------------------")
print("---------------------------***------------------------------")
print(" The figure 4A is successfully plotted ")
