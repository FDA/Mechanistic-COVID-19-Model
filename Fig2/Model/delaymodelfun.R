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
lag48_FreeVirus<-0
}else{
lag48_FreeVirus<- lagvalue(Time - 48,2)}

# Reactions----------------------------------
## Parameters are explained in "Table S3. The Full Mechanistic Model Parameters" 
ReactionFlux1=K14*EpitheliumwithreplicatedVirus
ReactionFlux2=K16*EpitheliumwithreplicatedVirus
ReactionFlux3=K13*EpitheliumwithVirusparticle
ReactionFlux4=K11*FreeVirus*EpitheliumLung
ReactionFlux5=K18*FreeVirus
ReactionFlux6=K48*FreeVirus*A
ReactionFlux7=K29*IL6
ReactionFlux8=K15*InfectedAPC+K49*EpitheliumwithreplicatedVirus
ReactionFlux14=K34*EpitheliumwithVirusparticle*K40*lag12_TcLymphnode
ReactionFlux15=K34*EpitheliumwithreplicatedVirus*K40*lag12_TcLymphnode
ReactionFlux18=K27*ResistantEpithelium
ReactionFlux188=K288*IL6*ResistantEpithelium
ReactionFlux19=(K24*K40*lag12_TcLymphnode+K47*EpitheliumwithreplicatedVirus+K57*InfectedAPC)*1
ReactionFlux20=K26*IFNA1*EpitheliumLung
ReactionFlux200=K288*IL6*EpitheliumLung
ReactionFlux21=K25*IFNA1
ReactionFlux22=K28*IFNA1*EpitheliumwithreplicatedVirus
ReactionFlux222=K288*IL6*EpitheliumwithreplicatedVirus
ReactionFlux23=K28*IFNA1*EpitheliumwithVirusparticle
ReactionFlux233=K288*IL6*EpitheliumwithVirusparticle
ReactionFlux25=K37*lag24_InfectedAPC
ReactionFlux26=K32*ActivatedAPCLymphnode*Tnaive
ReactionFlux266=K322*(1E5-Tnaive)
ReactionFlux2666=K3222*FreeVirus*Tnaive
ReactionFlux27=K42*ActivatedAPCLymphnode
ReactionFlux29=K53*InfectedAPC

ReactionFlux244=K45*TcLymphnode
ReactionFlux2444=K3222*FreeVirus*TcLymphnode
ReactionFlux66=K59*(1E5-A)
ReactionFlux88=K245*lag48_FreeVirus
ReactionFlux40=K400*8*2*TcLymphnode-K401*IgM
ReactionFlux50=K500*8*2*TcLymphnode-K501*IgG
ReactionFlux55=K555*FreeVirus*IgM+K556*FreeVirus*IgG
ReactionFlux9=Krenew*(4e8-EpitheliumLung-ResistantEpithelium-EpitheliumwithreplicatedVirus-EpitheliumwithVirusparticle)

# The differential equations----------------------------------
## Virus Life Cycle Module:-----------------------------------
# The annotated equations are coresponses to the supplementary documant equations.
dEpitheliumwithreplicatedVirus=1/unnamed*(-ReactionFlux1+ReactionFlux3-ReactionFlux15-ReactionFlux22-ReactionFlux222) #Equation (2)  
dFreeVirus=1/unnamed*(ReactionFlux2-ReactionFlux5-ReactionFlux55) #Equation (1)   
dEpitheliumwithVirusparticle=1/unnamed*(-ReactionFlux3+ReactionFlux4-ReactionFlux14-ReactionFlux23-ReactionFlux233) #Equation (3)  

## Immune systems Module:---------------------------------------------
dEpitheliumLung=1/unnamed*(-ReactionFlux4+ReactionFlux18-ReactionFlux20-ReactionFlux200+ReactionFlux9) #Equation (6)  
dInfectedAPC=1/unnamed*(ReactionFlux6-ReactionFlux29) #Equation (9)  
dIL6=1/unnamed*(-ReactionFlux7+ReactionFlux8+ReactionFlux88) #Equation (5)  
dResistantEpithelium=1/unnamed*(-ReactionFlux18+ReactionFlux20-ReactionFlux188) #Equation (7)  
dIFNA1=1/unnamed*(ReactionFlux19-ReactionFlux21) #Equation (4)  
dTcLymphnode=1/unnamed*(-ReactionFlux244+ReactionFlux26-ReactionFlux2444) #Equation (12)  
dActivatedAPCLymphnode=1/unnamed*(ReactionFlux25-ReactionFlux27) #Equation (10)  
dA=ReactionFlux66-ReactionFlux6 #Equation (8)  
dEpiKilledbyIFN=ReactionFlux22+ReactionFlux23 #Death rate of epithelial cells induced by IFN. (Calculated for checking)
dEpiKilledbyIL6=ReactionFlux222+ReactionFlux233+ReactionFlux200+ReactionFlux188 #Death rate of epithelial cells induced by IL6. (Calculated for checking)
dEpiKilledbyTc=ReactionFlux15+ReactionFlux14 #Death rate of epithelial cells induced by TCell. (Calculated for checking)
dEpiKilledbyVirus=ReactionFlux4 #Death rate of epithelial cells induced by Virus. (Calculated for checking)
dTnaive=ReactionFlux266-ReactionFlux2666 #Equation (11)  
dIgM=ReactionFlux40 #Equation (13)  
dIgG=ReactionFlux50 #Equation (14)  
    return(list(c(dEpitheliumwithreplicatedVirus,dFreeVirus,dEpitheliumwithVirusparticle,dEpitheliumLung,dInfectedAPC,dIL6,dResistantEpithelium,dIFNA1,dTcLymphnode,dActivatedAPCLymphnode,dA,dEpiKilledbyIFN,dEpiKilledbyIL6,dEpiKilledbyTc,dEpiKilledbyVirus,dTnaive,dIgG,dIgM)))
})}