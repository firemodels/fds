#!/bin/bash

CURDIR=`pwd`

# create smokebot directory

if [ ! -d ~/smokebot ] ; then
  cd 
  echo ~\smokebot does not exist - creating
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk/Utilities/Firebot smokebot
  echo ~\smokebot created.
fi

# create a clean cfast repository 

if [ ! -d ~/cfastclean ] ; then
  cd 
  echo ~\cfastclean does not exist - creating
  svn co http://cfast.googlecode.com/svn/trunk/cfast/trunk cfastclean
  echo ~\cfastclean created.
fi

# create a clean FDS-SMV repository 

if [ ! -d ~/FDS-SMVclean ] ; then
  cd 
  echo ~\FDS-SMVclean does not exist - creating
  svn co http://fds-smv.googlecode.com/svn/trunk/FDS/trunk FDS-SMVclean
  echo ~\FDS-SMVclean created.
fi

cd $CURDIR
