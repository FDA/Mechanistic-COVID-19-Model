#invitro part
ss="twoerror_aNewrun_Notscaled_10000_cov10_jumpis10per_initpar20percent_6hrpick_sigmaper2_1o5WC3edPointFinal"
mcmcpar=sprintf("/scratch/mohammadreza.samieegohar/Covid/MCMC_Newmodel_Knew/twoErro_step1_Newmodel_MoreNewmodelKnew_rev03_NoFint_KCellfix_justHAE_Pmax/results/%s/results/",ss)
tttmt<-read.csv(paste0(mcmcpar,"1apars.csv"),header=T)
tttmt=tttmt[,-1]
#initpar=encode_pars(ttt[,2])

ss2="new_testtwoerror_cov10_jumpis10per_initpar20percentPN_6hrpick_sigmaper1_dt100_BRAD1_LIKeHoo_JFenable"
mcmcpar2=sprintf("/scratch/mohammadreza.samieegohar/Covid/PlasmaPK/MCMC/results/%s/results/",ss2)
tttmt2<-read.csv(paste0(mcmcpar2,"4apars.csv"),header=T)
tttmt2=tttmt2[,-1]
tttmt_all=cbind(tttmt,tttmt2)
if (ip==0){
#write.csv(tttmt_all,sprintf("suplimantary_table_q95/%s_tttmt_all.csv",caseMS))
}
#set.seed(seedn)
#Day_After_syp=sample(seq(6,12), 1, replace = TRUE)
#admin_days=read.csv("Day_after.csv")
#admin_days=read.csv("dis17.csv")
# random day-----------------------------------S
set.seed(seedn)
r=rnorm(2*parnum_init, mean = 9, sd = 2.5)
r1=r[r>=6 & r<=12]
set.seed(seedn)
Day_After_sypd=sample(r1,parnum_init)
#----------------------------------------------E

#print("Day_After_syp")
#ip=0
#print("dddd")

#if (ip==0) {Day_After_syp=9}else{Day_After_syp=admin_days$x[ip]}
if (ip==0) {Day_After_syp=9}else{Day_After_syp=Day_After_sypd[ip]}

#print("Day_After_syp")
#print(Day_After_syp)
#xxxx
#uuuu
t3=dt3=c()
#if (caseMS=="mildReady") {
#	st_day=addday+Day_After_syp
#	
#}
#
#if (caseMS=="severReady") {
#	st_day=addday+Day_After_syp
#}
st_day=adddaysim+Day_After_syp
#print(addday)
#print(Day_After_syp)
#print("**************")

#st_day=3 
j=24 #each day
#ini_doseRDV=200 #mg
#doseRDV=100 #mg
infusion_time=2 #hr
num_day=9
#num_day=30

if (i1==1) {ini_doseRDV=0;doseRDV=0;}
if (i1==0) {ini_doseRDV=0;doseRDV=0;}

#if (i1==2) {ini_doseRDV=0;doseRDV=0;}

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
if (reme==1){eventtimes=eventtimes+1}else{eventtimes=eventtimes} #ineed to have aeven number for save calculation time
eventdose<-as.vector(dt3)*10^6  #kcn unit is ng consistat with pk model, in pk also we chang eit from mg to ng 
#print(eventtimes)
#print(st_day)



nanidx<-is.nan(eventdose)
eventdata<-data.frame(var="KCon",time=eventtimes[!nanidx],value=eventdose[!nanidx],method="add")

if (ip==0) {
	pars["Kp14"]=0.000272*3600 #is should change to 1/h
	pars["Kp41"]=0.00014*3600
	pars["KpC"]=0.002*3600
	pars["K1314"]=IC50_scale*0.1/(9119*1000) #unit mol/l  ~1E-8
	pars["Km"]=0.0002541*3600
	pars["F_intracellular"]=1
	pars["Kc"]=2.446e-05*3600
	pars["Kcell"]=0*3600
	pars["Knew"]=0.006602*3600
	pars["CL_percell"]=2.417e-06
	pars["F_extracellular"]=0.12
	pars["V_percell"]=1.17/4E8/1000
	pars["V1"]=6697 #ml
}else{
	timeconv=3600
#	if (ip>2000) {ip=sample(seq(1,2000), 1, replace = TRUE)}
	INTip=sample(seq(1,2000), parnum_init, replace = TRUE)
	INTipi=INTip[ip]
#	print("readyyyyyyy")
#	print(ip)
#	print(INTipi)
#	print("--------------------")
	
	pars["Kp14"]=tttmt_all[INTipi,"K14"]*timeconv #is should change to 1/h
	pars["Kp41"]=tttmt_all[INTipi,"K41"]*timeconv
	pars["KpC"]=tttmt_all[INTipi,"KC"]*timeconv
	pars["K1314"]=IC50_scale*0.1/(9119*1000) #unit mol/l  ~1E-8
	pars["Km"]=tttmt_all[INTipi,"Km"]*timeconv
	pars["F_intracellular"]=1
	pars["Kc"]=tttmt_all[INTipi,"Kc"]*timeconv
	pars["Kcell"]=0*timeconv
	pars["Knew"]=tttmt_all[INTipi,"Knew"]*timeconv
	pars["CL_percell"]=tttmt_all[INTipi,"CL_percell"]
	pars["F_extracellular"]=0.12
	pars["V_percell"]=1.17/4E8/1000
	pars["V1"]=tttmt_all[INTipi,"V1"]
	#pars["Kn1"]=.61
}
#write.csv(pars,sprintf("suplimantary_table_q95/%s_all.csv",caseMS))

initstates["KCon"]=doseRDV*10^6*0
initstates["cell_number"]=10^6
xxday=30 #it is for having atleast runtime equal to optm and for 9dosing we are cut them 
#timepoints=seq(0,((st_day+num_day)*24+xxday*24),2)
lasttim=max(c(((st_day+num_day)*24+xxday*24),650))
timepoints=seq(0,lasttim,2)
#timepoints=c(seq(0,lasttim,10),eventtimes)
#timepoints=sort(timepoints)
#timepoints=unique(timepoints)
if (ip==0) {timepoints=seq(0,650,1)}
#print(max(timepoints))
#if (max(eventtimes)>max(timepoints)) {print("time not right"); sstop}
#if (i1==1) {parnum=0}else{parnum=50}
