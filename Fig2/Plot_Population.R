#Ploting sizes------------------------------------------------------------
point_size=3
line_size=1.5
numsiz=16
labsiz=18
#Decrease T lymphocytes ----------------------------------------------------
out_plot1$yPred_TNAIVE_perc=(10^5-out_plot1$Tnaive)/10^5*100
BTnaiveq2=BTnaiveq1
BTnaiveq2$qmin=(10^5-BTnaiveq1$qmax)/10^5*100
BTnaiveq2$qmax=(10^5-BTnaiveq1$qmin)/10^5*100
ErrbarTn=ErrbarTn[ErrbarTn$time<=576,]
TNAIVE_percent1plot=TNAIVE_percent1plot[TNAIVE_percent1plot$time<=576,]

p<-ggplot()
p<-p+geom_point(data=TNAIVE_percent1plot, aes( x=time,y=Tnaivepercent), size=point_size, alpha=1, color="black")
p<-p+geom_line(data=out_plot1, aes(x=time, y=yPred_TNAIVE_perc), size=line_size, alpha=1,linetype = "solid")
p<-p+geom_errorbar(data=ErrbarTn,aes(x=time,ymin=ymin, ymax=ymax),color="black",linetype = "solid",width=20,
		position=position_dodge(20))
p<-p+geom_ribbon(data=BTnaiveq2, aes(x=time, ymin=qmax,ymax=qmin),alpha=0.3)
p<-p+ylab(paste0("Decrease % of \n T lymphocytes"))+ylim(0,100)+xlim(0,576)+scale_x_continuous(breaks=seq(0,576,8*24))+ theme_bw()
p<-p+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())

#IL6 ----------------------------------------------------
BIL6q1=data.frame(BIL6q1)
BIL6q1$qmin=BIL6q1$qmin/75
BIL6q1$qmax=BIL6q1$qmax/75
BIL6q1$time=BIL6q1$time
out_plot1$IL6=out_plot1$IL6/75
IL6$y1min[IL6$y1min<0]=0
IL6=IL6[IL6$time<=576,]

p4<-ggplot()
p4<-p4+geom_point(data=IL6, aes(x=time, y=y1mean), size=point_size, alpha=1, color="black")
p4<-p4+geom_line(data=out_plot1, aes(x=time, y=IL6), size=line_size, alpha=1,linetype = "solid")
p4<-p4+geom_ribbon(data=BIL6q1, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)
p4<-p4+geom_errorbar(data=IL6,aes(x=time,ymin=y1min, ymax=y1max),width=20,
		position=position_dodge(20),col="black",na.rm=TRUE)
p4<-p4+ylab(paste0("IL6 \n (pg/ml)"))+ylim(0,175)+xlim(0,576)+scale_x_continuous(breaks=seq(0,576,8*24))+ theme_bw()
p4<-p4+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
VIRUS$FreeVirus1=log10(VIRUS$yP)
out_plot1$FreeVirus=log10(out_plot1$FreeVirus)

#Virus Titer ----------------------------------------------------
VIRUS=VIRUS[VIRUS$time<=576,]
p1<-ggplot()
p1<-p1+geom_line(data=out_plot1, aes(x=time, y=FreeVirus), size=line_size, alpha=1,linetype = "solid")
p1<-p1+geom_point(data=VIRUS, aes(x=time, y=FreeVirus1), size=point_size, alpha=1, color="black")
p1<-p1+geom_ribbon(data=BFreeVirusq1, aes(x=time, ymin=(qmin),ymax=(qmax)),alpha=0.3)
p1<-p1+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title=element_text(size=15,face="bold"))+ theme_bw()
p1<-p1+ylab("Virus Titer \n (log[copy num/ml])")+ylim(0,9)+xlim(0,576)+scale_x_continuous(breaks=seq(0,576,8*24))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())

#IFN ----------------------------------------------------
IFNAalpha=data.frame(cbind(unique(IFNAalpha0$Day),unique(IFNAalpha0$INFA_mean_pgperml),unique(IFNAalpha0$min),unique(IFNAalpha0$max))) 
colnames(IFNAalpha)=c("Day","INFA_mean_pgperml","Imin","Imax") 
IFNAalpha=IFNAalpha[IFNAalpha$Day!=8,]
IFNAalpha$INFA_mean_pgperml=1*IFNAalpha$INFA_mean_pgperml    
IFNAalpha$Imin=1*IFNAalpha$Imin
IFNAalpha$Imax=1*IFNAalpha$Imax
IFA_sever=read.csv("Data_Covid/IFNAalpha/IFA_sever.csv")
IFA_sever=data.frame(IFA_sever)
IFA_sever$IFN_S=IFA_sever$IFN_S*1*.1
IFNAalpha=IFNAalpha[IFNAalpha$Day<23,] 		
IFNAalpha$Time=(IFNAalpha$Day+addday)*24
IFA_sever=IFA_sever[IFA_sever$Day<23,] 		
IFA_sever$Time=(IFA_sever$Day+addday)*24

#if (IFNAalpha$Imax[2]>300/10) {IFNAalpha$Imax[2]=300/10}
p5<-ggplot()
if (caseMS=="mildReady") {
	IFNAalpha=IFNAalpha[IFNAalpha$Time<=576,]
	p5<-p5+geom_point(data=IFNAalpha, aes(x=Time, y=INFA_mean_pgperml), size=point_size, alpha=1, color="black")
	p5<-p5+geom_errorbar(data=IFNAalpha,aes(x=Time,ymin=Imin, ymax=Imax),color="black",linetype = "solid",na.rm=TRUE, width=20,
			position=position_dodge(20))
}
if (caseMS=="severReady") {
	IFA_sever=IFA_sever[IFA_sever$Time<=576,]
	
	p5<-p5+geom_point(data=IFA_sever, aes(x=Time, y=IFN_S), size=point_size, alpha=1, color="black")	
}
lim=c(-1,500/10)

p5<-p5+geom_line(data=out_plot1, aes(x=time, y=IFNA1/12), size=line_size, alpha=1,linetype = "solid")
p5<-p5+geom_ribbon(data=BIFNA1q1, aes(x=time, ymin=qmin/12,ymax=qmax/12),alpha=0.3)
p5<-p5+ylab(paste0("IFN Alpha \n (pg/ml)"))+ylim(0,500/10)+xlim(0,576)+scale_x_continuous(breaks=seq(0,576,8*24))+scale_y_continuous(expand = c(0, 0), limits = lim) +
		theme_bw()
p5<-p5+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title.y=element_text(size=labsiz,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())

#IgM ----------------------------------------------------
NSIgM_data=NSIgM_data[NSIgM_data$Time<=576,]
p9<-ggplot()
p9<-p9+geom_point(data=NSIgM_data, aes(x=Time, y=mean), size=point_size, alpha=1, color="black")
p9<-p9+geom_errorbar(data=NSIgM_data,aes(x=Time,ymin=min, ymax=max),color="black",linetype = "solid",na.rm=TRUE,width=20,
		position=position_dodge(20))
p9<-p9+geom_line(data=out_plot1, aes(x=time, y=IgM), size=line_size, alpha=1,linetype = "solid")
p9<-p9+geom_ribbon(data=BIgMq1, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)+ylim(0,690)+xlim(0,576)+scale_x_continuous(breaks=seq(0,576,8*24))
p9<-p9+ylab(paste0("Titer of IgM \n "))+ theme_bw()
p9<-p9+xlab(paste0("Time post infection(hr)"))+theme(axis.text=element_text(size=numsiz),axis.title=element_text(size=20,face="bold"))+theme(axis.title.x = element_blank(), axis.text.x = element_blank())

#IgG ----------------------------------------------------
NSIgG_data=NSIgG_data[NSIgG_data$Time<=576,]
p10<-ggplot()
p10<-p10+geom_point(data=NSIgG_data, aes(x=Time, y=mean), size=point_size, alpha=1, color="black")
p10<-p10+geom_line(data=out_plot1, aes(x=time, y=IgG), size=line_size, alpha=1,linetype = "solid")
p10<-p10+geom_errorbar(data=NSIgG_data,aes(x=Time,ymin=min, ymax=max),color="black",linetype = "solid",na.rm=TRUE,width=20,
		position=position_dodge(20))
p10<-p10+geom_ribbon(data=BIgGq1, aes(x=time, ymin=qmin,ymax=qmax),alpha=0.3)
p10<-p10+ylab(paste0("Titer of IgG \n "))+ylim(0,69)+scale_x_continuous(breaks=seq(0,576,8*24))+ theme_bw()

p10<-p10+xlab(paste0("Time post infection(Day)"))+theme(axis.text=element_text(size=numsiz),axis.title=element_text(size=20,face="bold"))#+theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p10<-p10+scale_x_continuous( breaks= c(0,8*24,16*24,24*24),labels=c("0","8","16","24"))

#Saving plots --------------------------------------------------------------------------
pall=grid.arrange(arrangeGrob(p1,p5,p4,p,p9,p10, ncol=1))
biologicalfig=caserun
system(paste0("mkdir -p ",outdir2,biologicalfig))		
ggsave(paste0(outdir2,biologicalfig,"/",args$drug,"Fig2.png"),pall,height=16, width=8)
print("************************************************************")
print("--------------*******************************---------------")
print("------------------------*********---------------------------")
print("---------------------------***------------------------------")
print(" The population version of figure 2 is successfully plotted ")