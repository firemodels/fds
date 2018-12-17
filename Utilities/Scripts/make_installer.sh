#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: make_installer.sh -i FDS_TAR.tar.gz -d installdir INSTALLER.sh"
  echo ""
  echo "Creates an FDS/Smokeview installer sh script. "
  echo ""
  echo "  -i FDS.tar.gz - compressed tar file containing FDS distribution"
  echo "  -d installdir - default install directory"
  echo "   INSTALLER.sh - bash shell script containing self-extracting Installer"
  echo
  exit
fi

FDSEDITION=FDS6
FDSMODULE=$FDSEDITION

FDSVARS=${FDSEDITION}VARS.sh
INSTALLDIR=
FDS_TAR=
INSTALLER=
ostype="LINUX"
ostype2="Linux"
if [ "`uname`" == "Darwin" ] ; then
  ostype="OSX"
  ostype2="OSX"
fi

while getopts 'd:i:' OPTION
do
case $OPTION in
  d)
  INSTALLDIR="$OPTARG"
  ;;
  i)
  FDS_TAR="$OPTARG"
  ;;
esac
done 
shift $(($OPTIND-1))

INSTALLER=$1

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
if [ "$MPI_VERSION" != "INTEL" ]; then
OPENMPIFILE=openmpi_${MPI_VERSION}_${PLATFORM}_64.tar.gz
fi

cat << EOF > $INSTALLER
#!/bin/bash

OVERRIDE=\$1
INSTALL_LOG=/tmp/fds_install.log
echo "" > \$INSTALL_LOG
echo ""
echo "Installing FDS $FDSVERSION and Smokeview $SMVVERSION for $ostype2"
echo ""
echo "Options:"
echo "  1) Press <Enter> to begin installation [default]"
echo "  2) Type \"extract\" to copy the installation files to:"
echo "     $FDS_TAR"

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
          echo \$yn >> \$INSTALL_LOG
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
  fi
  touch \$DIR/temp.\$\$>&/dev/null
  if ! [ -e \$DIR/temp.\$\$ ]
  then
    echo ""
    echo "***error: \`whoami\` does not have permission to overwrite \$DIR"
    echo ""
    ls -ld \$DIR
    echo ""
    echo "Either: "
    echo "  1. change to a user that has permission, "
    echo "  2. remove \$DIR or,"
    echo "  3. change the owner/permissions of \$DIR" 
    echo "     to allow acceess to \`whoami\`"
    echo "FDS installation aborted."
    exit 0
  fi
  echo The installation directory, "\$DIR, has been created."
  rm \$DIR/temp.\$\$
}

#--- record the name of this script and the name of the directory 
#    it will run in

THISSCRIPT=\`ABSPATH \$0\`
THISDIR=\`pwd\`

#--- record temporary startup file names

BASHRCFDS=/tmp/bashrc_fds.\$\$
FDSMODULEtmp=/tmp/fds_module.\$\$
STARTUPtmp=/tmp/readme.\$\$

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
echo \$option >> \$INSTALL_LOG

if [ "\$option" == "extract" ]
then
  name=\$0
  THAT=$FDS_TAR
  if [ -e \$THAT ]
  then
    while true; do
      echo "The file, \$THAT, already exists."
      read -p "Do you wish to overwrite it? (yes/no) " yn
      echo \$yn >> \$INSTALL_LOG
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
echo \$answer >> \$INSTALL_LOG
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
eval MPIDIST_FDS=\$FDS_root/bin/openmpi_64
mpiused=\$FDS_root/bin/openmpi_64
eval MPIDIST_FDSROOT=\$FDS_root/bin
eval MPIDIST_FDS=\$FDS_root/bin/openmpi_64

#--- do we want to proceed

while true; do
   echo ""
   echo "Installation directory: \$FDS_root"
EOF
if [ "$OPENMPIFILE" != "" ]; then
cat << EOF >> $INSTALLER
   echo "     OpenMPI directory: \$mpiused"
EOF
fi
cat << EOF >> $INSTALLER
   if [ "\$OVERRIDE" == "y" ] ; then
     yn="y"
   else
     read -p "Proceed? (yes/no) " yn
   fi
   echo \$yn >> \$INSTALL_LOG
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

#--- copy installation files into the FDS_root directory

echo
echo "Copying FDS installation files to"  \$FDS_root
cd \$FDS_root
tail -n +\$SKIP \$THISSCRIPT | tar -xz
EOF

if [ "$OPENMPIFILE" != "" ]; then
cat << EOF >> $INSTALLER
if [ "\$MPIDIST_FDSROOT" != "" ]; then
  echo unpacking OpenMPI distribution to \$MPIDIST_FDSROOT
  cd \$MPIDIST_FDSROOT
  tar xvf $OPENMPIFILE >& /dev/null
fi
EOF
fi

cat << EOF >> $INSTALLER

echo "Copy complete."

#--- create fds module

MKDIR \$FDS_root/bin/modules

cat << MODULE > \$FDSMODULEtmp
#%Module1.0#####################################################################
###
### FDS6 modulefile
###

proc ModulesHelp { } {
        puts stderr "\tAdds FDS bin location to your PATH environment variable"
}

module-whatis   "Loads fds paths and libraries."

conflict FDS6
conflict openmpi
conflict intel

prepend-path    PATH            \$FDS_root/bin
prepend-path    LD_LIBRARY_PATH \$FDS_root/bin/LIB64
MODULE
if [ "$ostype" == "LINUX" ] ; then
cat << MODULE >> \$FDSMODULEtmp
prepend-path    LD_LIBRARY_PATH /usr/lib64
MODULE
fi
if [ "$MPI_VERSION" != "INTEL" ] ; then
cat << MODULE >> \$FDSMODULEtmp
prepend-path    PATH            \$FDS_root/bin/openmpi_64/bin
setenv          OPAL_PREFIX     \$FDS_root/bin/openmpi_64
MODULE
fi

cp \$FDSMODULEtmp \$FDS_root/bin/modules/$FDSMODULE
rm \$FDSMODULEtmp

#--- create BASH startup file

cat << BASH > \$BASHRCFDS
#/bin/bash
FDSBINDIR=\$FDS_root/bin
export PATH=\\\$FDSBINDIR:\\\$PATH
BASH

if [ "$MPI_VERSION" != "INTEL" ] ; then
cat << BASH >> \$BASHRCFDS
export PATH=\\\$FDSBINDIR/openmpi_64/bin:\\\$PATH
export OPAL_PREFIX=\\\$FDSBINDIR/openmpi_64
BASH
fi

if [ "$ostype" == "LINUX" ] ; then
OMP_COMMAND="grep -c processor /proc/cpuinfo"
else
OMP_COMMAND="system_profiler SPHardwareDataType"
fi

if [ "$ostype" == "LINUX" ] ; then
cat << BASH >> \$BASHRCFDS
export $LDLIBPATH=/usr/lib64:\\\$FDSBINDIR/LIB64:\\\$$LDLIBPATH
BASH
fi

cat << BASH >> \$BASHRCFDS
#  set OMP_NUM_THREADS to max of 4 and "Total Number of Cores" 
#  obtained from running:
#  \$OMP_COMMAND
export OMP_NUM_THREADS=4
BASH

#--- create startup and readme files

mv \$BASHRCFDS \$FDS_root/bin/$FDSVARS
chmod +x \$FDS_root/bin/$FDSVARS

#--- create startup readme file

cat << STARTUP > \$STARTUPtmp
<h3>Defining Environment Variables Used by FDS</h3>
Options:
<ul>
<li>Add following lines to one of your startup files
(usually \$HOME/.bashrc).<br>
<pre>
STARTUP

cat \$FDS_root/bin/$FDSVARS | grep -v bash >> \$STARTUPtmp

cat << STARTUP >> \$STARTUPtmp
</pre>
<li>or add:
<pre>
source \$FDS_root/bin/$FDSVARS
</pre>
<li>or if you are using modules, add:
<pre>
export MODULEPATH=\$FDS_root/bin/modules:\\\$MODULEPATH
module load $FDSMODULE
</pre>

<li>Log out and log back in so changes will take effect.

<li>To uninstall fds, erase the directory:<br>
\$FDS_root 
<p>and remove changes you made to your startup files.

STARTUP

mv \$STARTUPtmp \$FDS_root/Documentation/README_startup.html

EOF

cat << EOF >> $INSTALLER
echo ""
echo "-----------------------------------------------"
echo "*** To complete the installation add the following lines to your startup file"
echo "   (usually \$HOME/.bashrc)."
echo ""
cat \$FDS_root/bin/$FDSVARS | grep -v bash
echo ""
echo "or if you are using modules, add:"
echo ""
echo "export MODULEPATH=\$FDS_root/bin/modules:\\\$MODULEPATH"
echo "module load $FDSMODULE"
echo ""
echo "*** Log out and log back in so the changes will take effect."
echo ""
echo "*** To uninstall fds, erase the directory: "
echo "\$FDS_root"
echo "and remove any changes made to your startup file."
echo ""
echo "*** Installation complete."
exit 0


__TARFILE_FOLLOWS__
EOF
chmod +x $INSTALLER
cat $FDS_TAR >> $INSTALLER
echo "Installer created."
