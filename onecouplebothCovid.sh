#!/bin/sh
#$ -cwd
#$ -pe thread 1
#$ -j y
#$ -P CDERID0047
#$ -N Imp_Both_RDVBIO
#$ -R y
#$ -l h_rt=5:00:00
#$ -l h_vmem=4G
#$ -t 1-2        #num. drugs(28) x WORKERPERDRUG(1000)
#$ -o /dev/null
ERR=$(((SGE_TASK_ID-1)))
CASENAMES=(mild sever)
DRUG=${CASENAMES[ERR]}

source /projects/mikem/applications/R-4.0.2/set_env.sh



Rscript BothConnect_biological_to_RDVPK.R -d "$DRUG" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID".txt

