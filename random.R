isWindows<-Sys.info()[["sysname"]]=="Windows"
#--- load libraries
#----------------------
library(cmaes)
library(parallel)
library(ggplot2)
library(gridExtra)
plotlist = list()
mymodel<-list()
modelname<-"delaymymod"
modeldir<-"Data/"
nnnfile<-paste0(modeldir,"nnn.txt")
source("funs/model_init.R")
#--- load ODE model, states, and parameters
extension<-ifelse(isWindows, ".dll", ".so")
dyn.load(paste0(modeldir,modelname,extension))
source(paste0(modeldir,"delaypars.R"))
#---------Best case in Influenza
source(paste0(modeldir,"delaystates.R"))
inipath=""

if (i1==0) {parnum_init=4000}
number_p=2000
#parnum=5000
parnumq="parnumq"
cutday="daybase"
xuncS=.39 #already is define in side the code 
xuncM=.07 #1-xuncM is the NA cases that never recoverd

#sample(seq(0.65,1,0.02),1)

ML=.7;MM=1.15;SL=.7;SM=1.15
#SL=sample(seq(0.64,1,0.02),1)
#SM=sample(seq(1,1.2,0.02),1)

#ML=.5;MM=1.5;SL=.5;SM=1.5
#parnum_init_2=parnum_init*.01
onepat=0
#if (seedn>100 & seedn<200) {onepat=1}
#if (seedn>300 & seedn<400) {onepat=1}

parnum_init_2=parnum_init*0.05*onepat # 

subp="nnn_80pop_500step_100pmax_Dmage0o6_ty3_EComIFN"

outdiropt=sprintf("/scratch/mohammadreza.samieegohar/Covid/Rev06_Cofitting_DrLiBased/model09012020_10_rev03_nnn_adj_epthCellRecov_readme5_T153_improve_rev01/%s/T%s/results/",subp,iT)
outfolfig=sprintf("Virus2_0o5x_q90_1testT%s_twouc_xunSM%s_%s_%s_%s_L%s_%s_number_p%s_%s_parnum_init%s_parnum_init_2%s",iT,xuncS,xuncM,SL,SM,ML,MM,number_p,cutday,parnum_init,parnum_init_2)
#outdir2=paste0("New_AfterReview_2000out4000_yes550/",subp,"/",args$seed,"_",outfolfig,"/",IC50_scale,"/","/")
outdir2=paste0("BothCases_norestriction_reallyFinal_01/",subp,"/",args$seed,"_",outfolfig,"/",IC50_scale,"/","/")

system(paste0("mkdir -p ",outdir2))
pars_Changed_path=outdiropt
BCinf_pars0=read.csv(paste0(pars_Changed_path,"viruscytokines_pars.txt"),sep = "",header = F,stringsAsFactors =FALSE)
set.seed(seedn)
vec=BCinf_pars0$V2
par_rand1_1=par_rand2=c()
#--------------------------------------
iT2=42
subp="nnn_200pop_3000step_1000pmax_EConE3also_Topt115nnn_2xalsomils_try2"

outdiropt2=sprintf("/scratch/mohammadreza.samieegohar/Covid/Rev06_Cofitting_DrLiBased/model09012020_10_rev03_nnn_adj_epthCellRecov_readme5_T153_improve/%s/T%s/results/",subp,iT2)
pars_Changed_path2=outdiropt2
BCinf_pars02=read.csv(paste0(pars_Changed_path2,"viruscytokines_pars.txt"),sep = "",header = F,stringsAsFactors =FALSE)
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
print(BCinf_pars0)

