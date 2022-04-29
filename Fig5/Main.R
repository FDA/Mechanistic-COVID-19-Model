# Codes are written by: Mohammadreza Samieegohar
# Figure 5. Calibrating (a,b) and validating (c,d) the model for the primary endpoint (time to recovery) used in the remdesivir trial
#---------------------------------------------------
proc_start<-proc.time()
options(warn=-1)
# Load libraries -----------------------------------
library(optparse)
library(deSolve)
library(ggplot2)
library(gridExtra)
library(grid)
# Specify command line arguments-------------------
isWindows<-Sys.info()[["sysname"]]=="Windows"
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--severity"), default="mild",type="character", help="severity group,options:'mild' or 'severe'")
parser<-add_option(parser, c("-m", "--ModelPlatform"), default="R",type="character", help="if 'R' used, then the R version of the model will be used, which is compatible with Windows and Linux, otherwise a *so format is used that is compatible with Linux.")
parser<-add_option(parser, c("-p", "--popNum"), default=2000,type="integer", help="Virtual population size")
args<-parse_args(parser)
args$drug=args$severity
caserun="Population"
ModelPlatform=args$ModelPlatform
popNum=as.numeric(args$popNum)
Runcase=args$severity
print(sprintf("Severity is: -----%s-----",args$severity))
# Output folder ------------------------------------------------------------
outdir2="Fig5/"
system(paste0("mkdir -p ",outdir2))

 #Calculation of the Placebo group (Calibration Section) -------------------
print("Calibration is started----------------------------------")
source("Calibration.R")


print("Calibration is done and Validation is started-----------")

# Calculation of the Remdesivire group (Validation Section) ----------------
source("Validation.R")
print("Validation is done--------------------------------------")

# Subpopulation Sampling----------------------------------------------------
All_adm_t=c()
Allcases=c("mildReady","severReady")
if (Runcase== "mild")   {caserun=1}
if (Runcase== "severe") {caserun=2}
for (caseMS in Allcases[caserun]) {
	if (caseMS=="severReady") {casename="Severe";print(casename)}
	if (caseMS=="mildReady") {casename="Mild";print(casename)}
	minP=MedP=maxP=minD=MedD=maxD=0
	MedianVP=MedianVD=c()
	vi=3
	Damageper2placebo_V=Damageper2RMD_V=c()
	number_p=100
	i1=1;seedn=1
	for (seedn in 1) { 
		allparnum=c()
		for (iM in seq(1,100)) {			
			if (caseMS=="B" | caseMS=="mildReady"){
				severity="mildgroup" 
			}else{
				severity="severegroup" 				
			}			
			if (severity=="severegroup" | caseMS=="severReady")	 {	
					Sever_Control=read.csv("Data_Covid/Recovery/Sever_Control.csv")
					Sever_RMD=read.csv("Data_Covid/Recovery/Sever_RMD.csv")
					allparnumR=read.csv("Inputs_parameters/ID_Severe.csv")					
					cross5per1=cross5per1pre
					cross5per1$Roundedday=round(cross5per1$Recoveredt/24)
					cross5per1$NRoundedday=(cross5per1$Recoveredt/24)					
					cross5perR=cross5per1pre_RMD	
					cross5perR$Roundedday=round(cross5perR$Recoveredt/24)
					cross5perR$NRoundedday=(cross5perR$Recoveredt/24)			
				}						
				if (severity=="mildgroup" | caseMS=="mildReady")	 {				
					Sever_Control=read.csv("Data_Covid/Recovery/Mild_Control.csv")
					Sever_RMD=read.csv("Data_Covid/Recovery/Mild_RMD.csv")
					allparnumR=read.csv("Inputs_parameters/ID_Mild.csv")
					cross5per1=cross5per1pre
					cross5per1$Roundedday=round(cross5per1$Recoveredt/24)
					cross5per1$NRoundedday=(cross5per1$Recoveredt/24)					
					cross5perR=cross5per1pre_RMD
					cross5perR$Roundedday=round(cross5perR$Recoveredt/24)
					cross5perR$NRoundedday=(cross5perR$Recoveredt/24)
				}
			allparnumR1=allparnumR[,-1]
			parnum=allparnumR1[iM,]
			parnum=as.numeric(unlist(parnum))			
			cross5per2=cross5per1[cross5per1$ip%in%parnum,]
			cross5perR=cross5perR[cross5perR$ip%in%parnum,]				
			Damageper1=c()
				for (ita in seq(0,27,1)) {
				Damageper0=sum(round(cross5per2$Recoveredt/24)<=ita)/(number_p)*100
				Damageper1=rbind(Damageper1,c((ita)*24,Damageper0))
				}
			Damageper2=data.frame(Damageper1)
			colnames(Damageper2)=c("time","Recovered")
			Damageper2placebo=Damageper2	
			Damageper1R=c()
				for (ita in seq(0,27,1)) {
				Damageper0R=sum(round(cross5perR$Recoveredt/24)<=ita)/(number_p)*100
				Damageper1R=rbind(Damageper1R,c(ita*24,Damageper0R))
				}
			Damageper2R=data.frame(Damageper1R)
			colnames(Damageper2R)=c("time","Recovered")
			Damageper2RMD=Damageper2R			
			Damageper2placebo_V=rbind(Damageper2placebo_V,cbind(Damageper2placebo,iM,caseMS))
			Damageper2RMD_V=rbind(Damageper2RMD_V,cbind(Damageper2RMD,iM,caseMS))						
			RC_C=Sever_Control;RC_D=Sever_RMD;
			RC_C[,2:4]=RC_C[,2:4]*100 
			RC_D[,2:4]=RC_D[,2:4]*100 			
		}

# Ploting Calibration section Figure a,b ------------------------------------------
		MPall1=MDall1=c()
		for (ti1 in unique(Damageper2placebo_V$time)) {
			Pdf=Damageper2placebo_V$Recovered[Damageper2placebo_V$time==ti1]
			MPall=median(Pdf)
			PQw=quantile(Pdf, probs = c(0, 1))
			minPall=PQw[1]
			maxPall=PQw[2]
			MPall1=rbind(MPall1,c(ti1,MPall,minPall,maxPall))
			Ddf=Damageper2RMD_V$Recovered[Damageper2RMD_V$time==ti1]
			MDall=median(Ddf)
			DQw=quantile(Ddf, probs = c(0, 1))	
			minDall=DQw[1]
			maxDall=DQw[2]			
			MDall1=rbind(MDall1,c(ti1,MDall,minDall,maxDall))
		}		
		colnames(MPall1)=c("time","Recovered","min","max")
		colnames(MDall1)=c("time","Recovered","min","max")
		MDall1=data.frame(MDall1)
		MPall1=data.frame(MPall1)
		
# Ploting Validation section Figure c,d ------------------------------------------
numsiz=16
labsiz=18
		Rec_JustCtrl<-ggplot()+
				geom_line(data=MPall1, aes(x=time/24, y=Recovered ), size=.6, alpha=1,linetype = "solid",col="black")+
				geom_ribbon(data=MPall1,aes(x=time/24,ymin=min,ymax=max),alpha=0.3)+
				geom_errorbar(data=RC_C, aes(x=day, ymin=(min),ymax=(max)),width=.8,
						position=position_dodge(2),color="black",fill="black")+
				geom_point(data=RC_C, aes(x=day, y=opt), size=1.5, alpha=1,color="black")+
				xlab(paste0("Time(day)"))+ylab("Recovered%")+ylim(0,100)+ scale_x_continuous(breaks = seq(0, 27, by = 9),limits = c(0, 27.8))+
				theme(axis.text=element_text(size=20),axis.title=element_text(size=labsiz,face="bold"))+ theme_bw()
		ggsave(paste0(outdir2,casename,"_Placebo_Fig5.png"), Rec_JustCtrl, width=3.38, height=3.38)
				
		Recnoctr<-ggplot()+				
				geom_line(data=MDall1, aes(x=time/24, y=Recovered ), size=.6, alpha=1,linetype = "solid",col="black")+
				geom_ribbon(data=MDall1,aes(x=time/24,ymin=min,ymax=max),alpha=0.3)+
				geom_errorbar(data=RC_D, aes(x=day, ymin=(min),ymax=(max)),color="black",fill="black",width=.8,
						position=position_dodge(2),color="black",fill="black")+
				geom_point(data=RC_D, aes(x=day, y=opt), size=1.5, alpha=1,color="black")+
				xlab(paste0("Time(day)"))+ylab("Recovered%")+ylim(0,100)+xlim(0,27)+ scale_x_continuous(breaks = seq(0, 27, by = 9),limits = c(0, 27.8))+
				theme(axis.text=element_text(size=20),axis.title=element_text(size=labsiz,face="bold"))+ theme_bw()		
		ggsave(paste0(outdir2,casename,"_RMD_Fig5.png"), Recnoctr, width=3.38, height=3.38)
	}			
}	
if (casename=="Severe") {
	print("************************************************************")
	print("--------------*******************************---------------")
	print("------------------------*********---------------------------")
	print("---------------------------***------------------------------")
	print("----------------------------*-------------------------------")
	print("    figures 5c and 5d (severe) are successfully plotted     ")
}

if (casename=="Mild") {
	print("************************************************************")
	print("--------------*******************************---------------")
	print("------------------------*********---------------------------")
	print("---------------------------***------------------------------")
	print("----------------------------*-------------------------------")
	print("    figures 5a and 5b (mild) are successfully plotted       ")
}