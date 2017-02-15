#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: make_installer.sh -o ostype -i FDS_TAR.tar.gz -d installdir INSTALLER.sh"
  echo ""
  echo "Creates an FDS/Smokeview installer sh script. "
  echo ""
  echo "  -o ostype - OSX or LINUX"
  echo "  -i FDS.tar.gz - compressed tar file containing FDS distribution"
  echo "  -d installdir - default install directory"
  echo "   INSTALLER.sh - bash shell script containing self-extracting FDS installer"
  echo
  exit
fi

INSTALLDIR=
FDS_TAR=
ostype=
INSTALLER=

while getopts 'd:i:o:' OPTION
do
case $OPTION in
  d)
  INSTALLDIR="$OPTARG"
  ;;
  i)
  FDS_TAR="$OPTARG"
  ;;
  o)
  ostype="$OPTARG"
  ;;
esac
done 
shift $(($OPTIND-1))

INSTALLER=$1

if [ "$ostype" == "" ]
then
echo "*** fatal error: OS type (OSX or LINUX) not specified"
exit 0
fi

if [ "$FDS_TAR" == "" ]
then
echo "*** fatal error: FDS distribution file not specified"
exit 0
fi

if [ "$INSTALLDIR" == "" ]
then
echo "*** fatal error: default install directory not specified"
exit 0
fi

if [ "$INSTALLER" == "" ]
then
echo "*** fatal error: installer not specified"
exit 0
fi

BASHRC2=.bashrc
PLATFORM=linux
LDLIBPATH=LD_LIBRARY_PATH
if [ "$ostype" == "OSX" ]; then
  LDLIBPATH=DYLD_LIBRARY_PATH
  BASHRC2=.bash_profile
  PLATFORM=osx
fi
OPENMPIFILE=openmpi_1.8.4_${PLATFORM}_64.tar.gz

size2=64

ostype2=$ostype
if [ "$ostype" == "LINUX" ]
then
ostype2=Linux
fi

cat << EOF > $INSTALLER
#!/bin/bash

OVERRIDE=\$1
echo ""
echo "Installing $size2 bit $ostype2 FDS $FDSVERSION and Smokeview $SMVVERSION"
echo ""
echo "Options:"
echo "  1) Press <Enter> to begin installation"
echo "  2) Type \"extract\" to copy the installation files to $FDS_TAR"

BAK=_\`date +%Y%m%d_%H%M%S\`

#--- make a backup of a file

BACKUP_FILE()
{
  INFILE=\$1
  if [ -e \$INFILE ]
  then
  echo
  echo Backing up \$INFILE to \$INFILE\$BAK
  cp \$INFILE \$INFILE\$BAK
fi
}

#--- convert a path to it absolute equivalent

function ABSPATH() {
  pushd . > /dev/null;
  if [ -d "\$1" ];
  then
    cd "\$1";
    dirs -l +0;
  else
    cd "\`dirname \"\$1\"\`";
    cur_dir=\`dirs -l +0\`;
    if [ "\$cur_dir" == "/" ]; then
      echo "\$cur_dir\`basename \"\$1\"\`";
    else
      echo "\$cur_dir/\`basename \"\$1\"\`";
    fi;
  fi;
  popd > /dev/null;
}

#--- make a directory, checking if the user has permission to create it

MKDIR()
{
  DIR=\$1
  CHECK=\$2
  if [ ! -d \$DIR ]
  then
    echo "Creating directory \$DIR"
    mkdir -p \$DIR>&/dev/null
  else
    if [ "\$CHECK" == "1" ] 
    then
      while true; do
          echo "The directory, \$DIR, already exists."
          if [ "\$OVERRIDE" == "y" ]
            then
              yn="y"
          else
              read -p "Do you wish to overwrite it? (yes/no) " yn
          fi
          case \$yn in
              [Yy]* ) break;;
              [Nn]* ) echo "Installation cancelled";exit;;
              * ) echo "Please answer yes or no.";;
          esac
      done
      rm -rf \$DIR>&/dev/null
      mkdir -p \$DIR>&/dev/null
    fi
  fi
  if [ ! -d \$DIR ]
  then
    echo "Creation of \$DIR failed.  Likely cause,"
    echo "\`whoami\` does not have permission to create \$DIR."
    echo "FDS installation aborted."
    exit 0
  else
    echo The installation directory, "\$DIR, has been created."
  fi
  touch \$DIR/temp.\$\$>&/dev/null
  if ! [ -e \$DIR/temp.\$\$ ]
  then
    echo "\`whoami\` does not have permission to write to \$DIR"
    echo "FDS installation aborted."
    exit 0
  fi
  rm \$DIR/temp.\$\$
}

#--- record the name of this script and the name of the directory 
#    it will run in

THISSCRIPT=\`ABSPATH \$0\`
THISDIR=\`pwd\`

#--- record temporary startup file names

BASHFDS=/tmp/bashrc_fds.\$\$
BASHUNINSTALL=/tmp/uninstall_fds.\$\$

#--- Find the beginning of the included FDS tar file so that it 
#    can be subsequently un-tar'd
 
SKIP=\`awk '/^__TARFILE_FOLLOWS__/ { print NR + 1; exit 0; }' \$0\`

#--- extract tar.gz file from this script if 'extract' specified

if [ "\$OVERRIDE" == "y" ] 
then
  option=""
else
  read  option
fi

if [ "\$option" == "extract" ]
then
  name=\$0
  THAT=$FDS_TAR
  if [ -e \$THAT ]
  then
    while true; do
      echo "The file, \$THAT, already exists."
      read -p "Do you wish to overwrite it? (yes/no) " yn
      case \$yn in
        [Yy]* ) break;;
        [Nn]* ) echo "Extraction cancelled";exit;;
        * ) echo "Please answer yes or no.";;
      esac
    done
  fi
  echo Extracting the file embedded in this installer to \$THAT
  tail -n +\$SKIP \$THISSCRIPT > \$THAT
  exit 0
fi

OSSIZE=\`getconf LONG_BIT\`
if [ "\$OSSIZE" != "64" ] ; then
  if [ "\$OSSIZE" == "32" ] ; then
    echo "***Fatal error: FDS and Smokeview require a 64 bit operating system."
    echo "   The size of the operating system found is \$OSSIZE."
    exit 0
  fi
  echo "***Warning: FDS and Smokeview require a 64 bit operating system."
  echo "   The size of the operating system found is \$OSSIZE."
  echo "   Proceed with caution."
fi

#--- get FDS root directory

echo ""
echo "FDS install options"
EOF

if [ "$ostype" == "OSX" ]
then
cat << EOF >> $INSTALLER
    echo "  Press 1 to install in /Applications/$INSTALLDIR [default]"
    echo "  Press 2 to install in \$HOME/$INSTALLDIR"
EOF
  else
cat << EOF >> $INSTALLER
    echo "  Press 1 to install in \$HOME/$INSTALLDIR [default]"
    echo "  Press 2 to install in /opt/$INSTALLDIR"
    echo "  Press 3 to install in /usr/local/bin/$INSTALLDIR"
EOF
  fi
cat << EOF >> $INSTALLER
echo "  Enter a directory path to install elsewhere"

if [ "\$OVERRIDE" == "y" ] 
then
  answer="1"
else
  read answer
fi
EOF

if [ "$ostype" == "OSX" ]
then
cat << EOF >> $INSTALLER
  if [[ "\$answer" == "1" || "\$answer" == "" ]]; then
    eval FDS_root=/Applications/$INSTALLDIR
  elif [[ "\$answer" == "2" ]]; then
    eval FDS_root=\$HOME/$INSTALLDIR
  else
    eval FDS_root=\$answer
  fi
EOF
else
cat << EOF >> $INSTALLER
  if [[ "\$answer" == "1" || "\$answer" == "" ]]; then
    eval FDS_root=\$HOME/$INSTALLDIR
  elif [ "\$answer" == "2" ]; then
    FDS_root=/opt/$INSTALLDIR
  elif [ "\$answer" == "3" ]; then
    FDS_root=/usr/local/bin/$INSTALLDIR
  else
    eval FDS_root=\$answer
  fi
EOF
fi

#--- specify MPI location

cat << EOF >> $INSTALLER
valid_answer=
while true; do
  OPTION=0
  OPTION1=
  OPTION2=
  OPTION3=
  OPTION4=
  echo ""
  echo "OpenMPI options"

  OPTION=\$(echo \$OPTION + 1 | bc)
  OPTION2=\$OPTION
  echo "  Press \$OPTION2 to use \$FDS_root/bin/openmpi_64 [default]"

  OPTION=\$(echo \$OPTION + 1 | bc)
  OPTION1=\$OPTION
  echo "  Press \$OPTION1 to install OpenMPI manually"
  echo "     See \$FDS_root/bin/README.html for details"

  mpipath=
  mpipatheth=
  mpiused=
  if [ -d /shared/openmpi_64 ] ; then
     mpipath=\$MPIDIST_ETH
     mpipatheth=/shared/openmpi_64
     OPTION=\$(echo \$OPTION + 1 | bc)
     OPTION3=\$OPTION
     echo "  Press \$OPTION3 to use /shared/openmpi_64"
  fi
  mpipathib=
  if [ -d /shared/openmpi_64ib ] ; then
     mpipathib=/shared/openmpi_64ib
     mpipath=\$MPIDIST_IB
     OPTION=\$(echo \$OPTION + 1 | bc)
     OPTION4=\$OPTION
     echo "  Press \$OPTION4 to use /shared/openmpi_64ib"
  fi

  if [ "\$OVERRIDE" == "y" ]
  then
    answer="1"
  else
    read answer
  fi
  if [[ "\$answer" == "\$OPTION2" || "\$answer" == "" ]]; then
     answer=\$OPTION2
     eval MPIDIST_FDS=\$FDS_root/bin/openmpi_64
     eval MPIDIST_FDSROOT=\$FDS_root/bin
     mpiused=\$FDS_root/bin/openmpi_64
     valid_answer=1
  else
    eval MPIDIST_FDS=
  fi
  eval MPIDIST_FDS=\$FDS_root/bin/openmpi_64
  if [[ "\$answer" == "\$OPTION3" ]]; then
     mpipath2=\\\$MPIDIST_ETH
     mpiused=\$mpipatheth
     valid_answer=1
  fi
  if [[ "\$answer" == "\$OPTION4" ]]; then
     mpipath2=\\\$MPIDIST_IB
     mpiused=\$mpipathib
     valid_answer=1
  fi
  if [[ "\$valid_answer" == "" ]]; then
    echo ""
    echo "An invalid option was selected"
  else
    break;
  fi
done

mpipathfds=
if [ "\$MPIDIST_FDS" != "" ]; then
   mpipathfds=\$MPIDIST_FDS
   mpipath=\$MPIDIST_FDS
   if [[ "\$answer" == "\$OPTION2" ]]; then
     mpipath2=\\\$MPIDIST_FDS
   fi
fi

#--- do we want to proceed

while true; do
   echo ""
   echo "Installation directory: \$FDS_root"
   if [ "\$mpiused" == "" ] ; then
     echo "     OpenMPI directory: to be specified later" 
   else
     echo "     OpenMPI directory: \$mpiused"
   fi
   if [ "\$OVERRIDE" == "y" ] ; then
     yn="y"
   else
     read -p "Do you wish to proceed with the installation? (yes/no) " yn
   fi
   case \$yn in
      [Yy]* ) break;;
      [Nn]* ) echo "Installation cancelled";exit;;
      * ) echo "Please answer yes or no.";;
   esac
done
 
#--- make the FDS root directory

echo ""
echo "Installation beginning"
 
MKDIR \$FDS_root 1
MKDIR \$FDS_root/Uninstall 1

#--- copy installation files into the FDS_root directory

echo
echo "Copying FDS installation files to"  \$FDS_root
cd \$FDS_root
tail -n +\$SKIP \$THISSCRIPT | tar -xz
if [ "\$MPIDIST_FDSROOT" != "" ]; then
  echo unpacking OpenMPI distribution to \$MPIDIST_FDSROOT
  cd \$MPIDIST_FDSROOT
  tar xvf $OPENMPIFILE >& /dev/null
fi
echo "Copy complete."

#--- create uninstall file

cat << BASH > \$BASHUNINSTALL
#/bin/bash
FDSDIR=\$FDS_root
UNINSTALL=
BASH
if [ "$ostype" == "OSX" ] ; then
cat << BASH >> \$BASHUNINSTALL
BASHRC=~/.bash_profile
BASH
else
cat << BASH >> \$BASHUNINSTALL
BASHRC=~/.bashrc
BASH
fi
cat << BASH >> \$BASHUNINSTALL
while true; do
  read -p "Do you wish to remove \\\$FDSDIR, .bashrc_fds and FDS entries from \\\$BASHRC? (yes/no) " yn
  case \\\$yn in
      [Yy]* ) 
        UNINSTALL=1
        break;;
      [Nn]* ) 
        break;;
      * ) 
        echo "Please answer yes or no.";;
  esac
done
if [[ "\\\$UNINSTALL" == "1" ]]; then
  if [[ -d \\\$FDSDIR ]]; then
    echo removing \\\$FDSDIR
    rm -r \\\$FDSDIR
  else
    echo "***warning: The directory \\\$FDSDIR does not exist."
  fi
  if [[ -e \\\$BASHRC ]]; then
    BAK=_\\\`date +%Y%m%d_%H%M%S\\\`
    echo removing FDS entries from \\\$BASHRC
    grep -v bashrc_fds \\\$BASHRC | grep -v "#FDS" | grep -v INTEL_SHARED_LIB | grep -v MPIDIST_ETH | grep -v MPIDIST_FDS | grep -v MPIDIST_IB > ~/.bashrc_new
    mv \\\$BASHRC ~/.bashrc\\\$BAK
    mv ~/.bashrc_new \\\$BASHRC
  else
    echo "***warning: the file \\\$BASHRC does not exist."
  fi
  if [ -e ~/.bashrc_fds ]; then
    echo removing ~/.bashrc_fds
    rm ~/.bashrc_fds
  else
    echo "***warning: the file ~/.bashrc_fds does not exist."
  fi
  echo "Uninstall of FDS and Smokeview complete."
else
  echo "Uninstall of FDS and Smokeview cancelled."
fi
BASH

chmod +x \$BASHUNINSTALL
mv \$BASHUNINSTALL \$FDS_root/Uninstall/uninstall_fds.sh

#--- create BASH startup file

cat << BASH > \$BASHFDS
#/bin/bash

# OpenMPI location
ARG1=\\\$1

# Intel shared library location (default fds_install_dir/bin/INTELLIBS16)
ARG2=\\\$2

# FDS location

export FDSBINDIR=\$FDS_root/bin

# OpenMPI location

export MPIDIST=\\\$ARG1

# Intel shared library location

INTEL_SHARELIB=\\\$FDSBINDIR/INTELLIBS16
if [ "\\\$ARG2" != "" ]; then
  INTEL_SHARELIB=\\\$ARG2
fi

# unalias application names used by FDS

unalias fds >& /dev/null
unalias smokeview >& /dev/null
unalias smokezip >& /dev/null
unalias smokediff >& /dev/null
unalias fds6 >& /dev/null
unalias smokeview6 >& /dev/null
unalias smokezip6 >& /dev/null
unalias smokediff6 >& /dev/null

if [[ "\\\$MPIDIST" != "" && ! -d \\\$MPIDIST ]]; then
  echo "*** Warning: the MPI distribution, \\\$MPIDIST, does not exist"
  echo "      check the 'source ./bashrc_fds' line in your $BASHRC2 startup file"
  MPIDIST=
fi

# FDS network type

FDSNETWORK=
if [[ "\\\$MPIDIST" == *ib ]] ; then
  FDSNETWORK=infiniband
fi
export FDSNETWORK

# Update $LDLIBPATH

BASH
if [ "$ostype" == "LINUX" ] ; then
cat << BASH >> \$BASHFDS
$LDLIBPATH=\\\$FDSBINDIR/LIB64:\\\$INTEL_SHARELIB:\\\$$LDLIBPATH
BASH
fi
cat << BASH >> \$BASHFDS
if [ "\\\$MPIDIST" != "" ]; then
  $LDLIBPATH=\\\$MPIDIST/lib:\\\$$LDLIBPATH
fi

# Update PATH

if [ "\\\$MPIDIST" == "" ]; then
  PATH=\\\$FDSBINDIR:\\\$PATH
else
  PATH=\\\$MPIDIST/bin:\\\$FDSBINDIR:\\\$PATH
fi
export $LDLIBPATH PATH

# Set number of OMP threads

export OMP_NUM_THREADS=4
BASH

#--- create .bash_fds startup file

BACKUP_FILE ~/.bashrc_fds

if [ -e ~/.bashrc_fds ] ; then
  echo Updating .bashrc_fds
else
  echo Creating .bashrc_fds
fi
cp \$BASHFDS ~/.bashrc_fds
rm \$BASHFDS

#--- update .bash_profile
EOF
if [ "$ostype" == "OSX" ]; then
cat << EOF >> $INSTALLER
  BACKUP_FILE ~/.bash_profile

  BASHSTARTUP=/tmp/.bash_profile_temp_\$\$
  cd \$THISDIR
  echo "Updating .bash_profile"
  grep -v bashrc_fds ~/.bash_profile | grep -v INTEL_SHARED_LIB | grep -v "#FDS" | grep -v MPIDIST_ETH | grep -v MPIDIST_FDS | grep -v MPIDIST_IB > \$BASHSTARTUP
  echo "#FDS --------------------------------------------" >> \$BASHSTARTUP
  echo "export MPIDIST_FDS=\$mpipathfds"                >> \$BASHSTARTUP
  echo "export MPIDIST_ETH=\$mpipatheth"                >> \$BASHSTARTUP
  echo "export MPIDIST_IB=\$mpipathib"                  >> \$BASHSTARTUP
  echo "source ~/.bashrc_fds \$mpipath2"                >> \$BASHSTARTUP
  echo "#FDS --------------------------------------------" >> \$BASHSTARTUP
  cp \$BASHSTARTUP ~/.bash_profile
  rm \$BASHSTARTUP
EOF
else
cat << EOF >> $INSTALLER
#--- update .bashrc
  BACKUP_FILE ~/.bashrc

  BASHSTARTUP=/tmp/.bashrc_temp_\$\$
  cd \$THISDIR
  echo "Updating .bashrc"
  grep -v bashrc_fds ~/.bashrc | grep -v INTEL_SHARED_LIB | grep -v "#FDS" | grep -v MPIDIST_ETH | grep -v MPIDIST_FDS | grep -v MPIDIST_IB > \$BASHSTARTUP
  echo "#FDS -----------------------------------" >> \$BASHSTARTUP
  echo "export MPIDIST_FDS=\$mpipathfds"          >> \$BASHSTARTUP
  echo "export MPIDIST_ETH=\$mpipatheth"          >> \$BASHSTARTUP
  echo "export MPIDIST_IB=\$mpipathib"            >> \$BASHSTARTUP
if [ "\\\$IFORT_COMPILER_LIB" != "" ]; then
  echo "INTEL_SHARED_LIB=\\\$IFORT_COMPILER_LIB/intel64" >> \$BASHSTARTUP
else
  if [ "\\\$IFORT_COMPILER" != "" ]; then
    echo "INTEL_SHARED_LIB=\\\$IFORT_COMPILER/lib/intel64" >> \$BASHSTARTUP
  fi
fi
  echo "source ~/.bashrc_fds \$mpipath2"          >> \$BASHSTARTUP
  echo "#FDS -----------------------------------" >> \$BASHSTARTUP
  cp \$BASHSTARTUP ~/.bashrc
  rm \$BASHSTARTUP
EOF
fi

cat << EOF >> $INSTALLER
echo ""
echo "*** Log out and log back in so changes will take effect."
echo ""
echo "Installation complete."
exit 0


__TARFILE_FOLLOWS__
EOF
chmod +x $INSTALLER
cat $FDS_TAR >> $INSTALLER
echo "Installer created."
