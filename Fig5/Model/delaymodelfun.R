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
if (CR<1e-14) {CR=0}
if (CRP<1e-14) {CRP=0}
if (KCon<1e-14) {KCon=0}
if (C_intracellularRDV_total<1e-14) {C_intracellularRDV_total=0}
if (C_extracellularRDV_total<1e-14) {C_extracellularRDV_total=0}
if (cell_number<1e-14) {cell_number=0}
if (TP<1e-14) {TP=0}
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
ReactionFlux19=K24*K40*lag12_TcLymphnode+K47*EpitheliumwithreplicatedVirus+K57*InfectedAPC
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
ReactionFlux32=0
ReactionFlux244=K45*TcLymphnode
ReactionFlux2444=K3222*FreeVirus*TcLymphnode
ReactionFlux66=K59*(1E5-A)
ReactionFlux88=K245*lag48_FreeVirus
ReactionFlux40=K400*8*2*TcLymphnode-K401*IgM
ReactionFlux50=K500*8*2*TcLymphnode-K501*IgG
ReactionFlux55=K555*FreeVirus*IgM+K556*FreeVirus*IgG
ReactionFlux1311=-Kp14*CR+Kp41*CRP
ReactionFlux1312=-KpC*CR
ReactionFlux1313=KCon/2
ReactionFlux1314=(1/(1+(TP/K1314)^Kn1))
ReactionFlux1315=Km*(C_intracellularRDV_total*F_intracellular)
ReactionFlux1316=Kc*TP
ReactionFlux1317=Kcell*cell_number
ReactionFlux1318=Knew*C_intracellularRDV_total
ReactionFlux1319=CL_percell*(C_extracellularRDV_total*F_extracellular)/V_percell
ReactionFlux1320=CL_percell*(C_intracellularRDV_total*F_intracellular)/V_percell
ReactionFlux1321=((-Kp14*CR+Kp41*CRP)-KpC*CR+KCon/2)/V1
ReactionFlux9=Krenew*(4e8-EpitheliumLung-ResistantEpithelium-EpitheliumwithreplicatedVirus-EpitheliumwithVirusparticle)

dEpitheliumwithreplicatedVirus=1/unnamed*(-ReactionFlux1+ReactionFlux1314*ReactionFlux3-ReactionFlux15-ReactionFlux22-ReactionFlux222)
dFreeVirus=1/unnamed*(ReactionFlux2-ReactionFlux5-ReactionFlux55)
dEpitheliumwithVirusparticle=1/unnamed*(-ReactionFlux1314*ReactionFlux3+ReactionFlux4-ReactionFlux14-ReactionFlux23-ReactionFlux233)
dEpitheliumLung=1/unnamed*(-ReactionFlux4+ReactionFlux18-ReactionFlux20-ReactionFlux200+ReactionFlux9)
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
dEpiKilledbyVirus=ReactionFlux4
dTnaive=ReactionFlux266-ReactionFlux2666
dIgM=ReactionFlux40
dIgG=ReactionFlux50
dCR=ReactionFlux1311+ReactionFlux1312+ReactionFlux1313
dCRP=-(ReactionFlux1311)
dKCon=0
dC_intracellularRDV_total=ReactionFlux1319-ReactionFlux1320-ReactionFlux1315-ReactionFlux1318
dC_extracellularRDV_total=ReactionFlux1321*1000/1000000000/603
dcell_number=ReactionFlux1317
dTP=ReactionFlux1315-ReactionFlux1316
    return(list(c(dEpitheliumwithreplicatedVirus,dFreeVirus,dEpitheliumwithVirusparticle,dEpitheliumLung,dInfectedAPC,dIL6,dResistantEpithelium,dIFNA1,dTcLymphnode,dActivatedAPCLymphnode,dA,dEpiKilledbyIFN,dEpiKilledbyIL6,dEpiKilledbyTc,dEpiKilledbyVirus,dTnaive,dIgG,dIgM,dCR,dCRP,dKCon,dC_intracellularRDV_total,dC_extracellularRDV_total,dcell_number,dTP)))
})}