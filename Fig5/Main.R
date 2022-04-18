# All codes are written by: Mohammadreza Samieegohar
# Figure 5
#---------------------------------------------------
proc_start<-proc.time()
options(warn=-1)
library(optparse)
parser<-OptionParser()
parser<-add_option(parser, c("-c", "--Runcase"), default="mild",type="character", help="severe or mild")
args<-parse_args(parser)
Runcase=args$Runcase
library(ggplot2)
library(grid)
numsiz=16
labsiz=18
iT=00
DamP=0
outdir2="Fig5/"
system(paste0("mkdir -p ",outdir2))
byplot=""
EC=""
All_adm_t=c()
Allcases=c("mildReady","severReady","B","C","D","E")
allC1=c()
if (Runcase== "mild")   {caserun=1}
if (Runcase== "severe") {caserun=2}
for (caseMS in Allcases[caserun]) {
	if (caseMS=="severReady") {casename="Severe";print(casename)}
	if (caseMS=="mildReady") {casename="Mild";print(casename)}
	minP=MedP=maxP=minD=MedD=maxD=0
	MedianVP=MedianVD=c()
	vi=3
	Damageper2placebo_V=Damageper2RMD_V=c()
	LW1=1
	LW2=0
	i1=1;seedn=1
	for (seedn in 1) { 
		for (iM in seq(1,100)) {
			
			if (caseMS=="B" | caseMS=="mildReady"){
				severity="mildgroup" 
			}else{
				severity="severegroup" 				
			}			
			if (i1>0) {
				if (caseMS=="mildReady") { 
					Sever_Control=read.csv("Data_Covid/Recovery/Mild_Control.csv")
					Sever_RMD=read.csv("Data_Covid/Recovery/Mild_RMD.csv")
					d_diff=diff(Sever_Control[,3])
					d_diffmin=diff(Sever_Control[,2])
					d_diffmax=diff(Sever_Control[,4])
					number_p=100
					DamP=0
					LW1=1
					LW3=29
					LW2=0
				}
				if ( caseMS=="severReady") { 
					Sever_Control=read.csv("Data_Covid/Recovery/Sever_Control.csv")
					Sever_RMD=read.csv("Data_Covid/Recovery/Sever_RMD.csv")
					d_diff=diff(Sever_Control[,3])
					d_diffmin=diff(Sever_Control[,2])
					d_diffmax=diff(Sever_Control[,4])
					DamP=0
					number_p=100
					LW2=0
					LW1=1
					LW3=29
				}
				if ( caseMS=="B") { 
					Sever_Control=read.csv(sprintf("Data_Covid/Recovery/BCDE/%s_Control.csv","B"))
					Sever_RMD=read.csv(sprintf("Data_Covid/Recovery/BCDE/%s_RMD.csv","B"))
					d_diff=diff(Sever_Control[,3])
					d_diffmin=diff(Sever_Control[,2])
					d_diffmax=diff(Sever_Control[,4])
					DamP=0
					number_p=100
					LW1=1.00
					LW2=0
					LW3=29
				}
				if ( caseMS=="C") { 
					Sever_Control=read.csv(sprintf("Data_Covid/Recovery/BCDE/%s_Control.csv","C"))
					Sever_RMD=read.csv(sprintf("Data_Covid/Recovery/BCDE/%s_RMD.csv","C"))
					d_diff=diff(Sever_Control[,3])
					d_diffmin=diff(Sever_Control[,2])
					d_diffmax=diff(Sever_Control[,4])					
					LW1=1.00
					LW2=0
					number_p=100
					LW3=28
					DamP=0
				}
				if ( caseMS=="D") { 
					Sever_Control=read.csv(sprintf("Data_Covid/Recovery/BCDE/%s_Control.csv","D"))
					Sever_RMD=read.csv(sprintf("Data_Covid/Recovery/BCDE/%s_RMD.csv","D"))
					d_diff=diff(Sever_Control[,3])
					d_diffmin=diff(Sever_Control[,2])
					d_diffmax=diff(Sever_Control[,4])					
					LW1=1
					LW2=0
					number_p=100
					LW3=29
					DamP=15
				}
				if ( caseMS=="E") { 
					Sever_Control=read.csv(sprintf("Data_Covid/Recovery/BCDE/%s_Control.csv","E"))
					Sever_RMD=read.csv(sprintf("Data_Covid/Recovery/BCDE/%s_RMD.csv","E"))
					d_diff=diff(Sever_Control[,3])
					d_diffmin=diff(Sever_Control[,2])
					d_diffmax=diff(Sever_Control[,4])
					LW1=1
					LW2=0
					number_p=100
					LW3=29
					DamP=15
				}				
				if (severity=="severegroup")	 {					
					cross5per1=read.csv("Inputs_parameters/_i1is1_cross5per1_severReady.csv")
					cross5per1$Roundedday=round(cross5per1$Recoveredt/24)
					cross5per1$NRoundedday=(cross5per1$Recoveredt/24)					
					cross5perR=read.csv("Inputs_parameters/_i1is2_cross5per1_severReady.csv")	
					cross5perR$Roundedday=round(cross5perR$Recoveredt/24)
					cross5perR$NRoundedday=(cross5perR$Recoveredt/24)					
				}				
				if (severity=="mildgroup")	 {
					cross5per1=read.csv("Inputs_parameters/_i1is1_cross5per1_mildReady.csv")
					cross5per1$Roundedday=round(cross5per1$Recoveredt/24)
					cross5per1$NRoundedday=(cross5per1$Recoveredt/24)					
					cross5perR=read.csv("Inputs_parameters/_i1is2_cross5per1_mildReady.csv")	
					cross5perR$Roundedday=round(cross5perR$Recoveredt/24)
					cross5perR$NRoundedday=(cross5perR$Recoveredt/24)
				}

				wantip1=c();version="Nsimple";Er_round=0;round_correct=0
				for (id1 in 1:27) {
					whip=which(cross5per1$Roundedday==id1)
					if(iM==1) {
						sample_diff=sample(d_diffmin[id1],1)
						
					}
					if(iM==2) {
						sample_diff=sample(d_diffmax[id1],1)
					}										
					if(iM>2 ) {
						if( iM<=100/2+1) {
							sample_diff=sample(seq(d_diffmin[id1],d_diff[id1],(d_diff[id1]-d_diffmin[id1])/1000),1)
						}
						if(iM>=100/2-1) {
							sample_diff=sample(seq(d_diff[id1],d_diffmax[id1],(d_diffmax[id1]-d_diff[id1])/1000),1)
						}
					}
					if(sample_diff<0 | is.na(sample_diff)) {sample_diff=0}
					howmany_want=round(LW1*sample_diff*number_p+Er_round)+round_correct 			
					if (length(whip)==0 | is.na(howmany_want)) {next}
					needP=howmany_want-length(whip)
					if (needP==0) {rep_ip=cross5per1$ip[whip]}
					if (needP>0) {rep_ip=c(cross5per1$ip[whip],rep(cross5per1$ip[whip[1]],needP))}
					if (needP<0) {whip1=sample(whip,howmany_want);rep_ip=cross5per1$ip[whip1]}
					
					wantip=rep_ip
					wantip1=c(wantip1,wantip)
				}
			}
			if (length(wantip1)>number_p) {wantip1=wantip1[1:number_p] }
			ucov_PN=number_p-length(wantip1)			
			cross5per1=cross5per1[!is.na(cross5per1$NRoundedday),]
			deadman=0
			if (caseMS=="severReady") {deadman=0}
			Uc2=number_p-length(wantip1)-deadman			
			ucov_P2=cross5per1$ip[cross5per1$NRoundedday>LW3 & cross5per1$minevent>LW2 & cross5per1$damageD11>=DamP]			
			if (length(ucov_P2)==0) {ucov_P2=ucov_P; print("length ucov_P2 is zero");}
			ucov_Pt2=sample(ucov_P2,Uc2,replace=TRUE) 
			Uc1=0
			if (Uc1>0) {
				ucov_Pt1=sample(ucov_P,Uc1,replace=TRUE)
				parnum=c(wantip1,ucov_Pt1,ucov_Pt2)
			}else{
				parnum=c(wantip1,ucov_Pt2)
			}
			
			if (caseMS!="severReady" & caseMS!="mildReady") {
				cross5per2=cross5per1[cross5per1$ip%in%parnum,]
				cross5perR=cross5perR[cross5perR$ip%in%parnum,]				
			}else{				
				cross5per2=cross5per1[cross5per1$ip%in%parnum,]
				cross5perR=cross5perR[cross5perR$ip%in%parnum,]				
			}			
			Damageper1=c()
			for (ita in seq(0,27,1)) {
				Damageper0=sum(round(cross5per2$Recoveredt/24)<=ita)/(number_p)*100
				Damageper1=rbind(Damageper1,c((ita)*24,Damageper0))
			}
			Damageper2=data.frame(Damageper1)
			colnames(Damageper2)=c("time","Recovered")
			Damageper2placebo=Damageper2
			if ( caseMS=="E" ) {
			}			
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
			maxdd=1E8			
			Pgroup=cross5per2$Recoveredt[cross5per2$Recoveredt/24<maxdd]
			Dgroup=cross5perR$Recoveredt[cross5perR$Recoveredt/24<maxdd]			
			MedianVP[iM]=median((Pgroup/24))
			MedianVD[iM]=median((Dgroup/24))			
			minP=min(MedianVP,na.rm = TRUE)
			MedP=median(MedianVP,na.rm = TRUE)
			maxP=max(MedianVP,na.rm = TRUE)
			minD=min(MedianVD,na.rm = TRUE)
			MedD=median(MedianVD,na.rm = TRUE)
			maxD=max(MedianVD,na.rm = TRUE)
		}	
		allC=c(caseMS,number_p,number_p,MedP,paste0("(",minP,",",maxP,")"),MedD,paste0("(",minD,",",maxD,")"))		
		allC1=rbind(allC1,allC)		
		colN=c("Group","# Patient(Placebo)","# Patient(Remdesivir)","Median time to recovery(placebo)","95% CI (placebo)","Median time to recovery(Remdesivir)","95% CI (Remdesivir)")
		colnames(allC1)=colN	
		clinicalT=matrix(0,6,7)
		clinicalT[1,]=c("mild",50,55,"5","(4,7)","5","(4,6)")
		clinicalT[2,]=c("severe",471,486,"18","(15,20)","11","(10,14)")
		clinicalT[3,]=c("B",60,67,"6","(4,7)","5","(4,6)")
		clinicalT[4,]=c("C",199,222,"9","(7,10)","7","(6,8)")
		clinicalT[5,]=c("D",99,98,"20","(14-26)","15","(10,27)")
		clinicalT[6,]=c("E",147,125,"28","(24-NE)","29","(24,NE)")
		colnames(clinicalT)=colN
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
		Rec_JustCtrl<-ggplot()+
				geom_line(data=MPall1, aes(x=time/24, y=Recovered ), size=.6, alpha=1,linetype = "solid",col="black")+
				geom_ribbon(data=MPall1,aes(x=time/24,ymin=min,ymax=max),alpha=0.3)+
				geom_errorbar(data=RC_C, aes(x=day, ymin=(min),ymax=(max)),width=.8,
						position=position_dodge(2),color="black",fill="black")+
				geom_point(data=RC_C, aes(x=day, y=opt), size=1.5, alpha=1,color="black")+
				xlab(paste0("Time(day)"))+ylab("Recovered%")+ylim(0,100)+ scale_x_continuous(breaks = seq(0, 27, by = 9),limits = c(0, 27.8))+
				theme(axis.text=element_text(size=20),axis.title=element_text(size=labsiz,face="bold"))+ theme_bw()
		ggsave(paste0(outdir2,byplot,EC,casename,"_Placebo_Fig5.png"), Rec_JustCtrl, width=3.38, height=3.38)
				
		Recnoctr<-ggplot()+				
				geom_line(data=MDall1, aes(x=time/24, y=Recovered ), size=.6, alpha=1,linetype = "solid",col="black")+
				geom_ribbon(data=MDall1,aes(x=time/24,ymin=min,ymax=max),alpha=0.3)+
				geom_errorbar(data=RC_D, aes(x=day, ymin=(min),ymax=(max)),color="black",fill="black",width=.8,
						position=position_dodge(2),color="black",fill="black")+
				geom_point(data=RC_D, aes(x=day, y=opt), size=1.5, alpha=1,color="black")+
				xlab(paste0("Time(day)"))+ylab("Recovered%")+ylim(0,100)+xlim(0,27)+ scale_x_continuous(breaks = seq(0, 27, by = 9),limits = c(0, 27.8))+
				theme(axis.text=element_text(size=20),axis.title=element_text(size=labsiz,face="bold"))+ theme_bw()		
		ggsave(paste0(outdir2,byplot,EC,casename,"_RMD_Fig5.png"), Recnoctr, width=3.38, height=3.38)
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
