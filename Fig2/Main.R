# All codes are written by: Mohammadreza Samieegohar
# Figure 2
#---------------------------------------------------
proc_start<-proc.time()
options(warn=1)
library(optparse)
library(deSolve)
isWindows<-Sys.info()[["sysname"]]=="Windows"
parser<-OptionParser()
parser<-add_option(parser, c("-d", "--severity"), default="mild",type="character", help="severity group,options:'mild' or 'severe'")
parser<-add_option(parser, c("-c", "--case"), default="Optimum", help="options: 'Population' or 'Optimum'")
parser<-add_option(parser, c("-m", "--ModelPlatform"), default="R",type="character", help="if 'R' used, then the R version of the model will be used, which is compatible with Windows and Linux, otherwise a *so format is used that is compatible with Linux.")
parser<-add_option(parser, c("-p", "--popNum"), default=2000,type="integer", help="population number")
parser<-add_option(parser, c("-s", "--seed"), default=69, type="integer", help="Random seed")
args<-parse_args(parser)
args$drug=args$severity
caserun=args$case
ModelPlatform=args$ModelPlatform
popNum=as.numeric(args$popNum)
if (args$drug=="severe") {caseMS="severReady"}
if (args$drug=="mild") {caseMS="mildReady"}
seedn=as.numeric(args$seed)
args$ERR=1
iT=as.numeric(args$ERR)
print(sprintf("Severity is: -----%s-----",args$drug))
print(sprintf("Run is:    -----%s-----",args$case))
for (IC50_scale in c(1)) {		
	for (EC in 1 ) {
		byplot=""
		popis="false"
		for (i1 in 0:1) {
			source("Core_Code.R")
		}
	}
}

