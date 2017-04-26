#!/bin/bash
if [ $# -lt 1 ]; then
  echo "Usage: copycase.sh casename [HOST]"
  echo ""
  echo "The script copycase.sh creates a tar file for the FDS case with CHID"
  echo "casename (excluding the restart files).  The parameter HOST is optional."
  echo "If HOST is present it will copy the tar file to that host (for example"
  echo "blaze.el.nist.gov). If the parameter HOST is not present it will copy"
  echo "the tar file to the directory $HOME/TARCASES"
  exit
fi
CHID=$1
HOST=$2
TARDIR=$HOME/TARCASES
TARFILE=$TARDIR/${CHID}.tar
if [ ! -d $TARDIR ]; then
  mkdir $TARDIR
fi
tar -cvf $TARFILE --exclude='./*.restart' ./${CHID}*
if [ "$HOST" == "" ]; then
echo
echo FDS case: $CHID 
echo copied to: $TARFILE
else
scp -q $TARFILE $HOST\:.
echo FDS case: $CHID 
echo copied to host $HOST at ${CHID}.tar
fi
