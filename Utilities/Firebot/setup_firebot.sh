#!/bin/bash

CURDIR=`pwd`

if [ ! -d ~/firebot ] ; then
  cd 
  echo ~\firebot does not exist - creating
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Utilities/Firebot firebot
  echo ~\firebot created.
fi

# create a clean FDS repository

if [ ! -d ~/FDS-SMVclean ] ; then
  cd 
  echo ~\FDS-SMVclean does not exist - creating
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk FDS-SMVclean
  echo ~\FDS-SMVclean created.
fi

cd $CURDIR
