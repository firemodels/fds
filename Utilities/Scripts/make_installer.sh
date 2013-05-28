#!/bin/bash
EXPECTED_ARGS=5

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: make_installer.sh ostype ossize FDS_TAR.tar.gz INSTALLER.sh"
  echo ""
  echo "Creates an FDS/Smokeview installer sh script. "
  echo ""
  echo "  ostype - OSX or LINUX"
  echo "  ossize - ia32, intel64"
  echo "  FDS.tar.gz - compressed tar file containing FDS distribution"
  echo "  INSTALLER.sh - .sh script containing self-extracting installer"
  echo "  installdir - default install directory"
  echo
  exit
fi

ostype=$1
ossize=$2
FDS_TAR=$3
INSTALLER=$4
INSTALLDIR=$5

LDLIBPATH=LD_LIBRARY_PATH
if [ "$ostype" == "OSX" ]
then
LDLIBPATH=DYLD_LIBRARY_PATH
fi

cat << EOF > $INSTALLER
#!/bin/bash

echo ""
echo "FDS $FDSVERSION and Smokeview $SMVVERSION installer for $ostype $ossize"
echo ""
echo "Options:"
echo "  1) Press <Enter> to begin installation"
echo "  2) Type \"extract\" to copy the installation files to $FDS_TAR"

BAK=_\`date +%Y%m%d_%H%M%S\`

#--- output message

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

#--- create a script that points to an installed executable

MAKESHORTCUT()
{
  FDS_base=\$1
  EXE=\$2
  OUTFILE=\$3

cat << SHORTCUTFILE > \$TEMPFILE
#!/bin/bash 
\$FDS_base/bin/\$EXE \\\$@
SHORTCUTFILE

cp \$TEMPFILE \$OUTFILE
chmod +x \$OUTFILE
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
        read -p "Do you wish to overwrite it? (yes/no) " yn
        case \$yn in
            [Yy]* ) break;;
            [Nn]* ) echo "Installation cancelled";exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
    fi
  fi
  if [ ! -d \$DIR ]
  then
    echo "\`whoami\` does not have permission to create \$DIR."
    echo "FDS installation aborted."
    exit 0
  else
    echo The directory, "\$DIR, has been created."
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
# it will run in

THISSCRIPT=\`ABSPATH \$0\`
THISDIR=\`pwd\`

#--- record temporary start up file names

CSHFDS=/tmp/cshrc_fds.\$\$
BASHFDS=/tmp/bashrc_fds.\$\$

#--- Find the beginning of the included FDS tar file so that it can be
# subsequently un-tar'd
 
SKIP=\`awk '/^__TARFILE_FOLLOWS__/ { print NR + 1; exit 0; }' \$0\`

#--- extract tar.gz file from this script if 'extract' specified

read  option
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

#--- get FDS root directory

echo ""
echo "Where would you like to install FDS?"
  if [ "$ostype" == "OSX" ]; then
    echo "  Press 1 to install in /Applications/$INSTALLDIR"
  else
    echo "  Press 1 to install in \$HOME/$INSTALLDIR"
    echo "  Press 2 to install in /opt/$INSTALLDIR"
    echo "  Press 3 to install in /usr/local/bin/$INSTALLDIR"
  fi
echo "  Enter directory path to install elsewhere"
read answer

if [ "$ostype" == "OSX" ]; then
  if [[ "\$answer" == "1" || "\$answer" == "" ]]; then
    eval FDS_root=/Applications/$INSTALLDIR
  else
    eval FDS_root=\$answer
  fi
else
  if [[ "\$answer" == "1" || "\$answer" == "" ]]; then
    eval FDS_root=\$HOME/$INSTALLDIR
  elif [ "\$answer" == "2" ]; then
    FDS_root=/opt/$INSTALLDIR
  elif [ "\$answer" == "3" ]; then
    FDS_root=/usr/local/bin/$INSTALLDIR
  else
    eval FDS_root=\$answer
  fi
fi

#--- do we want to proceed

while true; do
    echo ""
    echo "Installation directory: \$FDS_root"
    read -p "Do you wish to begin the installation? (yes/no) " yn
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

#--- make shortcuts

SHORTCUTDIR=\`ABSPATH \$FDS_root/../shortcuts\`
MKDIR \$SHORTCUTDIR 0

TEMPFILE=/tmp/exetemp_\$\$
echo Adding shortcuts to \$SHORTCUTDIR
MAKESHORTCUT \$FDS_root fds \$SHORTCUTDIR/fds6
MAKESHORTCUT \$FDS_root smokeview \$SHORTCUTDIR/smokeview6
MAKESHORTCUT \$FDS_root smokediff \$SHORTCUTDIR/smokediff6
MAKESHORTCUT \$FDS_root smokezip \$SHORTCUTDIR/smokezip6
rm -f \$TEMPFILE

#--- copy installation files into the FDS_root directory

echo
echo "Copying FDS installation files to"  \$FDS_root
cd \$FDS_root
tail -n +\$SKIP \$THISSCRIPT | tar -xz
echo "Copy complete."

#--- create CSH startup file

cat << CSHRC > \$CSHFDS
#/bin/csh -f
set platform=\\\$1

# unalias application names used by FDS

unalias fds >& /dev/null
unalias smokeview >& /dev/null
unalias smokezip >& /dev/null
unalias smokediff >& /dev/null
unalias fds6 >& /dev/null
unalias smokeview6 >& /dev/null
unalias smokezip6 >& /dev/null
unalias smokediff6 >& /dev/null

# define FDS bin directory location

setenv FDSBINDIR \`pwd\`/bin

# environment for 64 bit Infiniband

if ( "\\\$platform" == "intel64ib" ) then
setenv MPIDIST /shared/openmpi_64ib
endif

# environment for 64 bit gigabit ethernet

if ( "\\\$platform" == "intel64" ) then
setenv MPIDIST /shared/openmpi_64
endif

# environment for 32 bit gigabit ethernet

if ( "\\\$platform" == "ia32" ) then
setenv MPIDIST \\\$MPIDIST32
endif

# Update LD_LIBRARY_PATH and PATH variables

setenv $LDLIBPATH \\\$MPIDIST/lib:\\\$$LDLIBPATH
set path=(\\\$FDSBINDIR \\\$MPIDIST/bin ~/bin \\\$path)

# if compilers are present then pre-define environment for their use

if ( \\\$?IFORT_COMPILER ) then

if ( -e \\\$IFORT_COMPILER/bin/compilervars.csh ) then

if ( "\\\$platform" == "intel64ib" ) then
source \\\$IFORT_COMPILER/bin/compilervars.csh intel64
endif

if ( "\\\$platform" == "intel64" ) then
source \\\$IFORT_COMPILER/bin/compilervars.csh intel64
endif

if ( "\\\$platform" == "ia32" ) then
source \\\$IFORT_COMPILER/bin/compilervars.csh ia32
endif

endif
endif

CSHRC

#--- create BASH startup file

cat << BASH > \$BASHFDS
#/bin/bash

platform=\\\$1

# unalias application names used by FDS

unalias fds >& /dev/null
unalias smokeview >& /dev/null
unalias smokezip >& /dev/null
unalias smokediff >& /dev/null
unalias fds6 >& /dev/null
unalias smokeview6 >& /dev/null
unalias smokezip6 >& /dev/null
unalias smokediff6 >& /dev/null

# define FDS bin directory location

export FDSBINDIR=\`pwd\`/bin
SHORTCUTDIR=\$SHORTCUTDIR

# environment for MPI

case "\\\$platform" in
  "intel64ib" )
    export MPIDIST=/shared/openmpi_64ib
  ;;
  "intel64" )
    export MPIDIST=/shared/openmpi_64
  ;;
  "ia32" )
    export MPIDIST=/shared/openmpi_32
  ;;
esac

# Update LD_LIBRARY_PATH and PATH variables

export $LDLIBPATH=\\\$MPIDIST/lib:\\\$$LDLIBPATH
export PATH=\\\$FDSBINDIR:\\\$SHORTCUTDIR:\\\$MPIDIST/bin:\\\$PATH

# if compilers are present then pre-define environment for their use

if [ -e "\\\$IFORT_COMPILER/bin/compilervars.sh" ]
then
case "\\\$platform" in
  "intel64" | "ia32" )
    source \\\$IFORT_COMPILER/bin/compilervars.sh \\\$platform
  ;;
  "intel64ib" )
    source \\\$IFORT_COMPILER/bin/compilervars.sh intel64
  ;;
esac
fi

BASH

#--- create .cshrc_fds startup file

echo

BACKUP_FILE ~/.cshrc_fds

if [ -e ~/.cshrc_fds ] ; then
  echo Updating .cshrc_fds
else
  echo Creating .cshrc_fds
fi
cp \$CSHFDS ~/.cshrc_fds
rm \$CSHFDS

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

if [ "$ostype" == "OSX" ]; then
  BACKUP_FILE ~/.bash_profile

  BASHPROFILETEMP=/tmp/.bash_profile_temp_\$\$
  cd \$THISDIR
  echo "Updating .bash_profile"
  grep -v bashrc_fds ~/.bash_profile | grep -v "#FDS" > \$BASHPROFILETEMP
  echo "#FDS " >> \$BASHPROFILETEMP
  echo "#FDS Setting environment for FDS and Smokeview.  The original version" >> \$BASHPROFILETEMP
  echo "#FDS of .bash_profile is saved in ~/.bash_profile\$BAK" >> \$BASHPROFILETEMP
  echo source \~/.bashrc_fds $ossize >> \$BASHPROFILETEMP
  cp \$BASHPROFILETEMP ~/.bash_profile
  rm \$BASHPROFILETEMP
fi

#--- update .bashrc

if [ "$ostype" != "OSX" ]; then
  BACKUP_FILE ~/.bashrc

  BASHRCTEMP=/tmp/.bashrc_temp_\$\$
  cd \$THISDIR
  echo "Updating .bashrc"
  grep -v bashrc_fds ~/.bashrc | grep -v "#FDS" > \$BASHRCTEMP
  echo "#FDS " >> \$BASHRCTEMP
  echo "#FDS Setting environment for FDS and Smokeview.  The original version" >> \$BASHRCTEMP
  echo "#FDS of .bashrc is saved in ~/.bashrc\$BAK" >> \$BASHRCTEMP
  echo source \~/.bashrc_fds $ossize >> \$BASHRCTEMP
  cp \$BASHRCTEMP ~/.bashrc
  rm \$BASHRCTEMP
fi

#--- update .cshrc

BACKUP_FILE ~/.cshrc

CSHTEMP=/tmp/.cshrc_temp_\$\$
echo "Updating .cshrc"
grep -v cshrc_fds ~/.cshrc | grep -v "#FDS" > \$CSHTEMP
echo "#FDS Setting environment for FDS and Smokeview.  The original version" >> \$CSHTEMP
echo "#FDS of .cshrc is saved in ~/.cshrc\$BAK" >> \$CSHTEMP
echo source \~/.cshrc_fds $ossize >> \$CSHTEMP
cp \$CSHTEMP ~/.cshrc
rm \$CSHTEMP

echo ""
echo "Installation complete."
exit 0


__TARFILE_FOLLOWS__
EOF
chmod +x $INSTALLER
cat $FDS_TAR >> $INSTALLER
echo "Installer created."
