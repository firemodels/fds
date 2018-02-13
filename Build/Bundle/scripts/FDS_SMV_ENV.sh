#/bin/bash 
  
# ---- FDS and smokeview versions

export fds_version=FDS6.6.0
export smv_version=SMV6.6.2
  
#  ---- MPI versions
  
export linux_mpi_version=INTEL 
export osx_mpi_version=3.0.0 

#  ---- Linux/OSX repo location
  
export linux_svn_root=FireModels_fork 

# ---------------------------------------
# shouldn't have to change anything below here

#  ---- Guide location

export GUIDE_DIR=$HOME/.bundle/FDS_Guides

#  ---- openmpi location

export OPENMPI_DIR=$HOME/.bundle/OPENMPI

#  ---- bundle location

export BUNDLE_DIR=$HOME/.bundle/BUNDLE
  
#  ---- bot locations
  
export firebotrepo=/home2/smokevis2/firebot/FireModels_central 
export firebothome=/home2/smokevis2/firebot 

export smokebotrepo=/home2/smokevis2/smokebot/FireModels_central 
export smokebothome=/home2/smokevis2/smokebot 
  
# ---- Linux login info
  
export linux_hostname=`hostname`
export linux_username=`whoami`
export linux_logon=$linux_username@$linux_hostname
  
# ---- OSX login info
  
export osx_hostname=`hostname`
export osx_username=`whoami`
export osx_logon=$osx_username@$osx_hostname
