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
echo "  2) Type \"extract\" to copy the installation files to"
echo "     the file $FDS_TAR"

FDS_root=$INSTALLDIR
BAK=\`date +%Y%m%d_%H%M%S\`

# record the name of this script and the name of the directory 
# it will run in

THIS=\`pwd\`/\$0
THISDIR=\`pwd\`

CSHFDS=/tmp/cshrc_fds.\$\$
BASHFDS=/tmp/bashrc_fds.\$\$

# Find the beginning of the included FDS tar file so that it can be
# subsequently un-tar'd
 
SKIP=\`awk '/^__TARFILE_FOLLOWS__/ { print NR + 1; exit 0; }' \$0\`

# extract tar.gz file from this script if 'extract' specified

read  option
if [ "\$option" == "extract" ]
then
name=\$0
THAT=$FDS_TAR
if [ -e \$THAT ]
then
while true; do
    echo "The file, \$THAT, already exists."
    read -p "Do you wish to overwrite it? (yes/no)" yn
    case \$yn in
        [Yy]* ) break;;
        [Nn]* ) echo "Extraction cancelled";exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
fi
echo Extracting the file embedded in this installer to \$THAT
tail -n +\$SKIP \$THIS > \$THAT
exit 0
fi

# get FDS root directory

echo ""
echo "Where would you like to install FDS? (default: \$FDS_root)"
read answer
if [ "\$answer" != "" ]; then
FDS_root=\$answer
fi
 
# make the FDS root directory
 
if [ ! -d \$FDS_root ]
then
echo "Creating directory \$FDS_root"
mkdir -p \$FDS_root>&/dev/null
else
while true; do
    echo "The directory, \$FDS_root, already exists."
    read -p "Do you wish to overwrite it? (yes/no)" yn
    case \$yn in
        [Yy]* ) break;;
        [Nn]* ) echo "Installation cancelled";exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
fi
  
HOMEBINDIR=~/bin
if [ ! -d \$HOMEBINDIR ]
then
mkdir \$HOMEBINDIR
fi

echo Creating fds6, smokeview6, smokediff6 and smokezip6 scripts

cat << FDS > \$HOMEBINDIR/fds6
#!/bin/bash 
\$FDS_root/bin/fds \\\$@
FDS
chmod +x \$HOMEBINDIR/fds6

cat << SMV > \$HOMEBINDIR/smokeview6
#!/bin/bash 
\$FDS_root/bin/smokeview \\\$@
SMV
chmod +x \$HOMEBINDIR/smokeview6

cat << SMV > \$HOMEBINDIR/smokediff6
#!/bin/bash 
\$FDS_root/bin/smokediff \\\$@
SMV
chmod +x \$HOMEBINDIR/smokediff6

cat << SMV > \$HOMEBINDIR/smokezip6
#!/bin/bash 
\$FDS_root/bin/smokezip \\\$@
SMV
chmod +x \$HOMEBINDIR/smokezip6

# did we succeed?

if [ ! -d \$FDS_root ]
then
echo "\`whoami\` does not have permission to create \$FDS_root."
echo "FDS installation aborted."
exit 0
else
echo "\$FDS_root has been created."
fi

# can we write to the FDS root directory?

touch \$FDS_root/temp.\$\$>&/dev/null
if ! [ -e \$FDS_root/temp.\$\$ ]
then
echo "\`whoami\` does not have permission to write to \$FDS_root."
echo "FDS installation aborted."
exit 0
fi
rm \$FDS_root/temp.\$\$

# now copy installation files into the FDS_root directory

echo
echo "Copying FDS installation files to"  \$FDS_root
cd \$FDS_root
tail -n +\$SKIP \$THIS | tar -xz
echo "Copy complete."

cat << CSHRC > \$CSHFDS
#/bin/csh -f
set platform=\\\$1

# unalias application names used by FDS

unalias fds >& /dev/null
unalias smokeview >& /dev/null
unalias smokezip >& /dev/null
unalias smokediff >& /dev/null

# define FDS bin directory location

setenv FDSBINDIR \`pwd\`/bin

# define openmpi library locations:
#   32/64 bit gigabit ethernet

set MPIDIST32=/shared/openmpi_32
set MPIDIST64=/shared/openmpi_64

# environment for 64 bit gigabit ethernet

if ( "\\\$platform" == "intel64" ) then
setenv MPIDIST \\\$MPIDIST64
set FORTLIB=\\\$FDSBINDIR/LIB64
endif

# environment for 32 bit gigabit ethernet

if ( "\\\$platform" == "ia32" ) then
setenv MPIDIST \\\$MPIDIST32
set FORTLIB=\\\$FDSBINDIR/LIB32
endif

if ( "\\\$platform" == "intel64" ) then
setenv $LDLIBPATH \\\`pwd\\\`/bin/LIB64
endif

if ( "\\\$platform" == "ia32" ) then
setenv $LDLIBPATH \\\`pwd\\\`/bin/LIB32
endif

# Update LD_LIBRARY_PATH and PATH variables

setenv $LDLIBPATH \\\$MPIDIST/lib:\\\${FORTLIB}:\\\$$LDLIBPATH
set path=(\\\$FDSBINDIR \\\$MPIDIST/bin ~/bin \\\$path)

# if compilers are present then pre-define environment for their use

if ( \\\$?IFORT_COMPILER ) then

if ( -e \\\$IFORT_COMPILER/bin/ifortvars.csh ) then

if ( "\\\$platform" == "intel64" ) then
source \\\$IFORT_COMPILER/bin/ifortvars.csh intel64
endif

if ( "\\\$platform" == "ia32" ) then
source \\\$IFORT_COMPILER/bin/ifortvars.csh ia32
endif

endif
endif

CSHRC

cat << BASH > \$BASHFDS
#/bin/bash

platform=\\\$1

# unalias application names used by FDS

unalias fds >& /dev/null
unalias smokeview >& /dev/null
unalias smokezip >& /dev/null
unalias smokediff >& /dev/null

# define FDS bin directory location

export FDSBINDIR=\`pwd\`/bin

# define openmpi library locations:
#   32/64 bit gigabit ethernet

MPIDIST32=/shared/openmpi_32
MPIDIST64=/shared/openmpi_64

# environment for 64 bit gigabit ethernet

case "\\\$platform" in
  "intel64" )
    export MPIDIST=\\\$MPIDIST64
    FORTLIB=\\\$FDSBINDIR/LIB64
  ;;
  "ia32" )
    export MPIDIST=\\\$MPIDIST32
    FORTLIB=\\\$FDSBINDIR/LIB32
  ;;
esac

# Update LD_LIBRARY_PATH and PATH variables

export $LDLIBPATH=\\\$MPIDIST/lib:\\\$FORTLIB:\\\$$LDLIBPATH
export PATH=\\\$FDSBINDIR:\\\$MPIDIST/bin:~/bin:\\\$PATH

# if compilers are present then pre-define environment for their use

if [ -e "\\\$IFORT_COMPILER/bin/ifortvars.sh" ]
then
case "\\\$platform" in
  "intel64" | "ia32" )
    source \\\$IFORT_COMPILER/bin/ifortvars.sh \\\$platform
  ;;
esac
fi

BASH

echo
echo Creating .bashrc_fds and .cshrc_fds startup files.
if [ -e ~/.cshrc_fds ]
then
echo Backing up .cshrc_fds to ~/.cshrc_fds_\$BAK
cp ~/.cshrc_fds ~/.cshrc_fds_\$BAK
fi

if [ -e ~/.bashrc_fds ]
then
echo Backing up .bashrc_fds to ~/.bashrc_fds_\$BAK
cp ~/.bashrc_fds ~/.bashrc_fds_\$BAK
fi

cp \$BASHFDS ~/.bashrc_fds
rm \$BASHFDS

cp \$CSHFDS ~/.cshrc_fds
rm \$CSHFDS

if [ -e ~/.bash_profile ]
then
echo Backing up .bash_profile to ~/.bash_profile_\$BAK
cp ~/.bash_profile ~/.bash_profile_\$BAK
fi

BASHPROFILETEMP=/tmp/.bash_profile_temp_\$\$
cd \$THISDIR
echo "Updating .bash_profile"
grep -v bashrc_fds ~/.bash_profile | grep -v "#FDS" > \$BASHPROFILETEMP
echo "#FDS Setting environment for FDS and Smokeview.  The original version" >> \$BASHPROFILETEMP
echo "#FDS of .bash_profile is saves in ~/.bash_profile_\$BAK" >> \$BASHPROFILETEMP
echo source \~/.bashrc_fds $ossize >> \$BASHPROFILETEMP
cp \$BASHPROFILETEMP ~/.bash_profile

touch ~/.cshrc
nlines=\$(grep cshrc_fds ~/.cshrc | wc -l)
if [ \$nlines -eq 0 ]
then
echo Updating .cshrc
echo source \~/.cshrc_fds $ossize >> ~/.cshrc
fi
echo ""
echo "Installation complete."
exit 0


__TARFILE_FOLLOWS__
EOF
chmod +x $INSTALLER
cat $FDS_TAR >> $INSTALLER
echo "Installer created."
