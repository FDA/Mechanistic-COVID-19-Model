modelfun <- function(Time, State, Pars){
currenttime<-unclass(as.POSIXct(strptime(date(),"%c")))[1] 
if(currenttime - Pars["starttime"]>=Pars["timeout"])
 stop("timeout!"); 
  with(as.list(c(State, Pars)), {
if(Time<12){
lag12_TcLymphnode<-0
}else{
lag12_TcLymphnode<- lagvalue(Time - 12,9)}
if(Time<24){
lag24_InfectedAPC<-0
}else{
lag24_InfectedAPC<- lagvalue(Time - 24,5)}
if(Time<48){
lag48_FreeInfluenza<-0
}else{
lag48_FreeInfluenza<- lagvalue(Time - 48,2)}
ReactionFlux1=K14*Epitheliumwithreplicatedinfluenza
ReactionFlux2=K16*Epitheliumwithreplicatedinfluenza
ReactionFlux3=K13*Epitheliumwithinfluenzaparticle
ReactionFlux4=K11*FreeInfluenza*EpitheliumLung
ReactionFlux5=K18*FreeInfluenza
ReactionFlux6=K48*FreeInfluenza*A
ReactionFlux7=K29*IL6
ReactionFlux8=K15*InfectedAPC+K49*Epitheliumwithreplicatedinfluenza
ReactionFlux14=K34*Epitheliumwithinfluenzaparticle*0.15*lag12_TcLymphnode
ReactionFlux15=K34*Epitheliumwithreplicatedinfluenza*0.15*lag12_TcLymphnode
ReactionFlux18=K27*ResistantEpithelium
ReactionFlux188=K288*IL6*ResistantEpithelium
ReactionFlux19=K24*0.15*lag12_TcLymphnode+K47*Epitheliumwithreplicatedinfluenza+K57*InfectedAPC
ReactionFlux20=K26*IFNA1*EpitheliumLung
ReactionFlux200=K288*IL6*EpitheliumLung
ReactionFlux21=K25*IFNA1
ReactionFlux22=K28*IFNA1*Epitheliumwithreplicatedinfluenza
ReactionFlux222=K288*IL6*Epitheliumwithreplicatedinfluenza
ReactionFlux23=K28*IFNA1*Epitheliumwithinfluenzaparticle
ReactionFlux233=K288*IL6*Epitheliumwithinfluenzaparticle
ReactionFlux25=K37*lag24_InfectedAPC
ReactionFlux26=K32*ActivatedAPCLymphnode*Tnaive
ReactionFlux266=K322*(1E5-Tnaive)
ReactionFlux2666=K3222*FreeInfluenza*Tnaive
ReactionFlux27=K42*ActivatedAPCLymphnode
ReactionFlux29=K53*InfectedAPC
ReactionFlux32=0
ReactionFlux244=K45*TcLymphnode
ReactionFlux2444=K3222*FreeInfluenza*TcLymphnode
ReactionFlux66=K59*(1E5-A)
ReactionFlux88=K245*lag48_FreeInfluenza
dEpitheliumwithreplicatedinfluenza=1/unnamed*(-ReactionFlux1+ReactionFlux3-ReactionFlux15-ReactionFlux22-ReactionFlux222)
dFreeInfluenza=1/unnamed*(ReactionFlux2-ReactionFlux5)
dEpitheliumwithinfluenzaparticle=1/unnamed*(-ReactionFlux3+ReactionFlux4-ReactionFlux14-ReactionFlux23-ReactionFlux233)
dEpitheliumLung=1/unnamed*(-ReactionFlux4+ReactionFlux18-ReactionFlux20-ReactionFlux200)
dInfectedAPC=1/unnamed*(ReactionFlux6-ReactionFlux29-ReactionFlux32)
dIL6=1/unnamed*(-ReactionFlux7+ReactionFlux8+ReactionFlux88)
dResistantEpithelium=1/unnamed*(-ReactionFlux18+ReactionFlux20-ReactionFlux188)
dIFNA1=1/unnamed*(ReactionFlux19-ReactionFlux21)
dTcLymphnode=1/unnamed*(-ReactionFlux244+ReactionFlux26-ReactionFlux2444)
dActivatedAPCLymphnode=1/unnamed*(ReactionFlux25-ReactionFlux27)
dA=ReactionFlux66-ReactionFlux6
dEpiKilledbyIFN=ReactionFlux22+ReactionFlux23
dEpiKilledbyIL6=ReactionFlux222+ReactionFlux233+ReactionFlux200+ReactionFlux188
dEpiKilledbyTc=ReactionFlux15+ReactionFlux14
dTnaive=ReactionFlux266-ReactionFlux2666
    return(list(c(dEpitheliumwithreplicatedinfluenza,dFreeInfluenza,dEpitheliumwithinfluenzaparticle,dEpitheliumLung,dInfectedAPC,dIL6,dResistantEpithelium,dIFNA1,dTcLymphnode,dActivatedAPCLymphnode,dA,dEpiKilledbyIFN,dEpiKilledbyIL6,dEpiKilledbyTc,dTnaive)))
})}