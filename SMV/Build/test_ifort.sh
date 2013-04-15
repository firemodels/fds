#!/bin/bash
if [ "$IFORT_COMPILER" == "" ] ; then
  echo "*** fatal error: The shell variable IFORT_COMPILER is not defined"
  echo "    It should point to the Intel compiler location"
  echo "    For example: /opt/intel/composerxe"
  echo "    Define IFORT_COMPILER in one of your shell startup scripts"
  echo "    (bash --> .bashrc : IFORT_COMPILER=/opt/intel/composerxe"
  echo "    (csh/tcsh -> .cshrc :  setenv IFORT_COMPILER /opt/intel/composerxe"
  exit
fi
