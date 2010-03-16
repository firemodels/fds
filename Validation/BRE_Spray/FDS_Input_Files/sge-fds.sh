#!/bin/bash
#$ -S /bin/bash
#$ -cwd -N BRE_Spray -V -e /dev/null -o /dev/null
#
# Put your Job commands here.
#------------------------------------------------
$FDS $in
