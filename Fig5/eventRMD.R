Day_After_syp=Par_population0[ip,ncol(Par_population0)]
incDay=Par_population0[ip,(ncol(Par_population0)-1)]
adddaysim=incDay
i1=2
if (i1==2) {ini_doseRDV=200;doseRDV=100;}
if (ip==0) {Day_After_syp=9}else{Day_After_syp=Day_After_syp}
st_day=adddaysim+Day_After_syp
j=24 
infusion_time=2 
num_day=9
t3=dt3=c()
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
