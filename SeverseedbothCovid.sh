#!/bin/sh
#$ -cwd
#$ -pe thread 1
#$ -j y
#$ -P CDERID0047
#$ -N Seed
#$ -l s_rt=24:00:00
#$ -R y
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00
#$ -t 1-100
#$ -o NULL
CASENAMES=(sever)
WORKERPERDRUG=100
ERR=$(((SGE_TASK_ID)))
IDX=$(((SGE_TASK_ID-1)%$WORKERPERDRUG+1))
IDX2=$(((SGE_TASK_ID-IDX)/$WORKERPERDRUG))
DRUG=${CASENAMES[IDX2]}
source /projects/mikem/applications/R-4.0.2/set_env.sh

Rscript BothConnect_biological_to_RDVPK.R -s "$IDX" -d "$DRUG" >& logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID".txt

