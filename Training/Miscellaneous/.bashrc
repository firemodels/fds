# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi
ulimit -s unlimited

# User specific environment and startup programs

export IFORT_COMPILER=/opt/intel/composerxe
export IFORT_COMPILER_LIB=$IFORT_COMPILER/lib

. /usr/local/Modules/3.2.10/init/bash

module load null modules torque-maui

MPIDIST=/shared/openmpi_64ib

export FDSNETWORK=infiniband

alias qfds.sh="~/FDS-SMV/Utilities/Scripts/qfds.sh"
LD_LIBRARY_PATH=$MPIDIST/lib:$LD_LIBRARY_PATH
PATH=$MPIDIST/bin:$PATH
export LD_LIBRARY_PATH PATH

