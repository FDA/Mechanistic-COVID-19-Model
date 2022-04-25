#!/bin/sh
#$ -cwd
#$ -pe thread 1
#$ -j y
#$ -N hERG_boot
#$ -l s_rt=240:00:00
#$ -R y
#$ -l h_rt=240:00:00
#$ -t 1-32000
#$ -o logfiles

echo "Running task $SGE_TASK_ID of job $JOB_NAME (ID of $JOB_ID) on $HOSTNAME"

source /projects/mikem/applications/centos7/python3/set-run-env.sh

time python3 helloworld.py    
