library(ggplot2)
library(gridExtra)
plotlist = list()
mymodel<-list()
modelname<-"delaymymod"
modeldir<-"Model/"
nnnfile<-paste0(modeldir,"nnn.txt")
source("funs/model_init.R")
source(paste0(modeldir,"delaypars.R"))
source(paste0(modeldir,"delaystates.R"))
inipath=""
if (caserun!="Optimum") {
	if (i1==0) {parnum_init=2*popNum}
	number_p=popNum	
}else{
	if (i1==0) {parnum_init=10}
	number_p=10	
}
parnumq="parnumq"
cutday="daybase"
xuncS=.39 
xuncM=.07
ML=.7;MM=1.15;SL=.7;SM=1.15
onepat=0
parnum_init_2=parnum_init*0.05*onepat 
#subp="nnn_80pop_500step_100pmax_Dmage0o6_ty3_EComIFN"
outdir2=paste0("Figs/","","/","","","","/","","/","/")
system(paste0("mkdir -p ",outdir2))
BCinf_pars0=read.csv(paste0('Inputs_parameters/',"viruscytokines_pars.txt"),sep = "",header = F,stringsAsFactors =FALSE)
set.seed(seedn)
vec=BCinf_pars0$V2
par_rand1_1=par_rand2=c()

iT2=42
BCinf_pars02=read.csv(paste0('Inputs_parameters/',"viruscytokines_pars.txt"),sep = "",header = F,stringsAsFactors =FALSE)
set.seed(seedn)
vec2=BCinf_pars02$V2
parnum_init_1=parnum_init-parnum_init_2
for (i in 1:length(vec)) {
	if (caseMS=="severReady") {par_rand=runif(parnum_init_1,SL*vec[i],SM*vec[i]);par_rand_2=runif(parnum_init_2,SL*vec2[i],SM*vec2[i])}
	if (caseMS=="mildReady") {par_rand=runif(parnum_init_1,ML*vec[i],MM*vec[i]);par_rand_2=runif(parnum_init_2,ML*vec2[i],MM*vec2[i])}
	par_rand1_1=cbind(par_rand1_1,par_rand)
	par_rand2=cbind(par_rand2,par_rand_2)
	
}
if (length(par_rand2)!=0){
	par_rand1=rbind(par_rand1_1,par_rand2)
}else{par_rand1=(par_rand1_1)}
colnames(par_rand1)=BCinf_pars0$V1
