#/bin/bash 

#   
if [ "$fds_version" == "" ]; then
  echo "***warning: fds_version is not defined"
fi
if [ "$smv_version" == "" ]; then
  echo "***warning: smv_version is not defined"
fi
  
if [ "$linux_mpi_version" == "" ]; then
  echo "***warning: linux_mpi_version is not defined"
fi
if [ "$osx_mpi_version" == "" ]; then
  echo "***warning: osx_mpi_version is not defined"
fi
  
if [ "$linux_svn_root" == "" ]; then
  echo "***warning: linux_svn_root is not defined"
fi

if [ "$GUIDE_DIR" == "" ]; then
  echo "***warning: GUIDE_DIR is not defined"
fi
  
if [ "$firebotrepo" == "" ]; then
  echo "***warning: firebotrepo is not defined"
fi
if [ "$firebothome" == "" ]; then
  echo "***warning: firebothome is not defined"
fi
if [ "$smokebotrepo" == "" ]; then
  echo "***warning: smokebotrepo is not defined"
fi
if [ "$smokebothome" == "" ]; then
  echo "***warning: smokebothome is not defined"
fi

if [ "$linux_hostname" == "" ]; then
  echo "***warning: linux_hostname is not defined"
fi
if [ "$linux_username" == "" ]; then
  echo "***warning: linux_username is not defined"
fi
if [ "$linux_logon" == "" ]; then
  echo "***warning: linux_logon is not defined"
fi
  
if [ "$osx_hostname" == "" ]; then
  echo "***warning: osx_hostname is not defined"
fi
if [ "$osx_username" == "" ]; then
  echo "***warning: osx_username is not defined"
fi
if [ "$osx_logon" == "" ]; then
  echo "***warning: osx_logon is not defined"
fi
