#!/bin/bash
#$ -S /bin/bash
#$ -cwd -N BRE_spray -V -e /dev/null -o /dev/null
#
# Start the job as qsub -t 1-24 sge-fds-array.sh
#
# Put your Job commands here.
#------------------------------------------------
$FDS BRE_Spray_$SGE_TASK_ID.fds
