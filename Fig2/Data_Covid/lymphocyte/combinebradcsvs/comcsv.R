datafile="C:/Users/Mohammad.Samieegohar/Desktop/Dr_Li/Brad_paper/combinebradcsvs/phase1/"
datadf1<-read.csv(paste0(datafile,"AstraZeneca_herg.csv"))
#d3=unique(datadf1$drug)

datadf2<-read.csv(paste0(datafile,"BristolMeyersSquibb_herg.csv"))
datadf3<-read.csv(paste0(datafile,"BSYSGmbH_herg.csv"))
datadf4<-read.csv(paste0(datafile,"CRLChantestRamp_herg.csv"))
datadf5<-read.csv(paste0(datafile,"Eurofins_herg.csv"))
datadf6<-read.csv(paste0(datafile,"Merck_herg.csv"))
datadf7<-read.csv(paste0(datafile,"Metrion_herg.csv"))
datadf8<-read.csv(paste0(datafile,"NanionMunichAmbientPatchliner_herg.csv"))
datadf9<-read.csv(paste0(datafile,"NanionMunichAmbientSyncropatch_herg.csv"))
datadf10<-read.csv(paste0(datafile,"NanionUSAmbientSyncropatch_herg.csv"))
datadf11<-read.csv(paste0(datafile,"NMI_herg.csv"))
datadf12<-read.csv(paste0(datafile,"Roche_herg.csv"))
datadf13<-read.csv(paste0(datafile,"Sophion_herg.csv"))


datafile2="C:/Users/Mohammad.Samieegohar/Desktop/Dr_Li/Brad_paper/combinebradcsvs/phase2/"
datadf14<-read.csv(paste0(datafile2,"BristolMeyersSquibb_herg.csv"))
datadf15<-read.csv(paste0(datafile2,"BSYSGmbH_herg.csv"))
datadf16<-read.csv(paste0(datafile2,"CRLChantestRamp_herg.csv"))
datadf17<-read.csv(paste0(datafile2,"Eurofins_herg.csv"))
datadf18<-read.csv(paste0(datafile2,"Merck_herg.csv"))
datadf19<-read.csv(paste0(datafile2,"Metrion_herg.csv"))
datadf20<-read.csv(paste0(datafile2,"NanionPatchlinerAmbient_herg.csv"))
datadf21<-read.csv(paste0(datafile2,"NanionSynchropatchAmbient_herg.csv"))
datadf22<-read.csv(paste0(datafile2,"Sophion_herg.csv"))

cc0=c("AstraZeneca","BristolMeyersSquibb","BSYSGmbH","CRLChantestRamp")



datadf1[,6]=1
datadf2[,6]=2
datadf3[,6]=3
datadf4[,6]=4
datadf5[,6]=5
datadf6[,6]=6
datadf7[,6]=7
datadf8[,6]=8
datadf9[,6]=9
datadf10[,6]=10
datadf11[,6]=11
datadf12[,6]=12
datadf13[,6]=13
datadf14[,6]=14
datadf15[,6]=15
datadf16[,6]=4
datadf17[,6]=17
datadf18[,6]=18
datadf19[,6]=7
datadf20[,6]=20
datadf21[,6]=21
datadf22[,6]=22


datadfall=rbind(datadf1,datadf2,datadf3,datadf4,datadf5,datadf6,datadf7,datadf8,datadf9,datadf10,datadf11,datadf12,datadf13,datadf14,datadf15,datadf16,datadf17,datadf18,datadf19,datadf20,datadf21,datadf22)


datadfall[,3]=datadfall[,6]
datadfall[,6]=datadfall[,2]
datadfall[,2]=datadfall[,5]
datadfall[,5]=datadfall[,6]
datadfall[,6]=datadfall[,4]
datadfall[,4]=datadfall[,5]
datadfall[,5]=datadfall[,6]


datadfall[,6]=NULL
datadfalum=datadfall
datadfalum[,4]=datadfalum[,4]/1000

colnames(datadfalum)=c("Compound","Channel","Experiment","Dose","Response")
colnames(datadfall)=c("Compound","Channel","Experiment","Dose","Response")

#m=unique.array(datadfall$Compound)
nams=levels(unlist(datadfall$Compound))
m=1
numexp=list()
for (i in nams){
my_data= subset(datadfall, Compound == "i")
#d1=levels(my_data$Experiment)
expe=unique(unlist(my_data$Experiment))
numexp=length(expe)
print(i)}





write.csv(datadfalum,"allumsz.csv")
write.csv(datadfall,"alls.csv")


    
    
  

