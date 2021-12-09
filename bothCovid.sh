#!/bin/sh
#$ -cwd
#$ -pe thread 1
#$ -j y
#$ -P CDERID0047
#$ -N Both_RDVBIO
#$ -R y
#$ -l h_rt=5:00:00
#$ -l h_vmem=4G
#$ -t 1-80        #num. drugs(28) x WORKERPERDRUG(1000)
#$ -o /dev/null
CASENAMES=(mild sever)
WORKERPERDRUG=40
ERR=$(((SGE_TASK_ID)))
IDX=$(((SGE_TASK_ID-1)%$WORKERPERDRUG+1))
IDX2=$(((SGE_TASK_ID-IDX)/$WORKERPERDRUG))
DRUG=${CASENAMES[IDX2]}

Rscript BothConnect_biological_to_RDVPK.R -e "$IDX" -d "$DRUG" >& /scratch/mohammadreza.samieegohar/new/Connect_biological_to_RDVPK_population_incubDayinMCMC_epthCellRecov1_readme4_Rev01/logfiles/"$JOB_NAME".o"$JOB_ID"."$SGE_TASK_ID".txt

