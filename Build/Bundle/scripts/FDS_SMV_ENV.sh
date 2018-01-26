#/bin/bash 
  
# ---- FDS and smokeview version ---- 
  
export fds_version=6.6.0test 
export smv_version=6.6.1test 
  
#  ---- repo locations ---- 
  
  
#  *** Linux/OSX 
  
export linux_svn_root=FireModels_fork 
export compiler_dir=fire-notes/INSTALL/LINUX/INTEL_17u4 
export misc_dir=fire-notes/INSTALL/LIBS/LINUX/LIB64 
  
#  ---- MPI library locations ---- 
  
export linux_mpi_version=INTEL 
export osx_mpi_version=3.0.0 
  
#  ---- bot locations ---- 
  
export firebotrepo=/home2/smokevis2/firebot/FireModels_central 
export firebothome=/home2/smokevis2/firebot 
export smokebotrepo=/home2/smokevis2/smokebot/FireModels_central 
export smokebothome=/home2/smokevis2/smokebot 
  
#  ---- hostnames ---- 
  
  
#  *** linux 
  
export linux_hostname=blaze.el.nist.gov 
export linux_username=`whoami`
export linux_logon=$linux_username@$linux_hostname
  
#  *** osx 
  
export osx_hostname=ignis.el.nist.gov 
export osx_username=`whoami`
export osx_logon=$osx_username@$osx_hostname
