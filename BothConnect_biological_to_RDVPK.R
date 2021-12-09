

proc_start<-proc.time()

options(warn=1)

library(optparse)

isWindows<-Sys.info()[["sysname"]]=="Windows"

#--- specify command line arguments
parser<-OptionParser()

parser<-add_option(parser, c("-d", "--drug"), default="sever",type="character", help="Drug name [required]")
parser<-add_option(parser, c("-o", "--opioids"), default="viruscytokines",type="character", help="opioids name [required]")
parser<-add_option(parser, c("-i", "--sample"), default="0", help="Bootstrap sample number, range, and/or list of samples to be fitted [default 0 fits to original data]")
parser<-add_option(parser, c("-e", "--ERR"), default="104", help="xxxrap sample number, range, and/or list of samples to be fitted [default 0 fits to original data]")
parser<-add_option(parser, c("-s", "--seed"), default=214, type="integer", help="Random seed [default 100]")
parser<-add_option(parser, c("-c", "--cores"), default=1, type="integer", help="Number of cores to use during fitting [default 1]")
parser<-add_option(parser, c("-f", "--forking"), default=FALSE, action="store_true", help="Flag to turn on forking for parallelization (not supported in Windows)")
parser<-add_option(parser, c("-l", "--lambda"), default=103,type="integer", help="Population size to use for CMA-ES [default 4+floor(3*log(N))]")
parser<-add_option(parser, c("-m", "--maxiter"), default=3000,type="integer", help="Maximum number of generations for CMA-ES [default 100*N^2]")
parser<-add_option(parser, c("-t", "--tol"), default=.01,type="double", help="Stopping tolerance for CMA-ES [default 0.5e-12]")

args<-parse_args(parser)

#--- load libraries
library(cmaes)
library(parallel)
library(deSolve)
if (args$drug=="sever") {caseMS="severReady"}
if (args$drug=="mild") {caseMS="mildReady"}
seedn=as.numeric(args$seed)
print(as.numeric(args$ERR))
print(as.numeric(args$seed))
iT=as.numeric(args$ERR)
print("caseIMS is:")
print(caseMS)
for (IC50_scale in c(1)) {		
	for (EC in 1 ) {
		byplot=""
		popis="false"
		for (i1 in 0:2) {
					source("Bothindependent_paperplot_Covid_fitting.R")

		}

	}
}
