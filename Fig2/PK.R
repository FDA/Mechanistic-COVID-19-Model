tttmt<-read.csv(paste0("Inputs_parameters/","1apars.csv"),header=T)
tttmt=tttmt[,-1]
tttmt2<-read.csv(paste0("Inputs_parameters/","4apars.csv"),header=T)
tttmt2=tttmt2[,-1]
tttmt_all=cbind(tttmt,tttmt2)
if (ip==0){
}
set.seed(seedn)
r=rnorm(2*parnum_init, mean = 9, sd = 2.5)
r1=r[r>=6 & r<=12]

set.seed(seedn)
Day_After_sypd=sample(r1,parnum_init)
if (ip==0) {Day_After_syp=9}else{Day_After_syp=Day_After_sypd[ip]}
t3=dt3=c()
st_day=adddaysim+Day_After_syp
j=24 
infusion_time=2
num_day=9
if (i1==1) {ini_doseRDV=0;doseRDV=0;}
if (i1==0) {ini_doseRDV=0;doseRDV=0;}
if (i1==2) {ini_doseRDV=200;doseRDV=100;}
for(i in st_day:(st_day+num_day)) {
	t2=c(i*j,(i*j+infusion_time))
	t3=c(t3,t2)
	if (i==st_day) {dt2=c(1,-1)*ini_doseRDV}
	if (i>st_day) {	dt2=c(1,-1)*doseRDV}
	dt3=c(dt3,dt2)
}

eventtimes<-round(as.vector(t3))
reme=eventtimes[1]%%2
if (reme==1){eventtimes=eventtimes+1}else{eventtimes=eventtimes} 
eventdose<-as.vector(dt3)*10^6  
nanidx<-is.nan(eventdose)
eventdata<-data.frame(var="KCon",time=eventtimes[!nanidx],value=eventdose[!nanidx],method="add")
if (ip==0) {
	pars["Kp14"]=0.000272*3600 
	pars["Kp41"]=0.00014*3600
	pars["KpC"]=0.002*3600
	pars["K1314"]=IC50_scale*0.1/(9119*1000) 
	pars["Km"]=0.0002541*3600
	pars["F_intracellular"]=1
	pars["Kc"]=2.446e-05*3600
	pars["Kcell"]=0*3600
	pars["Knew"]=0.006602*3600
	pars["CL_percell"]=2.417e-06
	pars["F_extracellular"]=0.12
	pars["V_percell"]=1.17/4E8/1000
	pars["V1"]=6697 
}else{
	timeconv=3600
	INTip=sample(seq(1,2000), parnum_init, replace = TRUE)
	INTipi=INTip[ip]
	pars["Kp14"]=tttmt_all[INTipi,"K14"]*timeconv 
	pars["Kp41"]=tttmt_all[INTipi,"K41"]*timeconv
	pars["KpC"]=tttmt_all[INTipi,"KC"]*timeconv
	pars["K1314"]=IC50_scale*0.1/(9119*1000) 
	pars["Km"]=tttmt_all[INTipi,"Km"]*timeconv
	pars["F_intracellular"]=1
	pars["Kc"]=tttmt_all[INTipi,"Kc"]*timeconv
	pars["Kcell"]=0*timeconv
	pars["Knew"]=tttmt_all[INTipi,"Knew"]*timeconv
	pars["CL_percell"]=tttmt_all[INTipi,"CL_percell"]
	pars["F_extracellular"]=0.12
	pars["V_percell"]=1.17/4E8/1000
	pars["V1"]=tttmt_all[INTipi,"V1"]	
}
initstates["KCon"]=doseRDV*10^6*0
initstates["cell_number"]=10^6
xxday=30 
lasttim=max(c(((st_day+num_day)*24+xxday*24),650))
timepoints=seq(0,lasttim,2)
if (ip==0) {timepoints=seq(0,650,1)}