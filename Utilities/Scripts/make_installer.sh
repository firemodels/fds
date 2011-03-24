#!/bin/bash
EXPECTED_ARGS=4

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: make_installer.sh platform FDS_TAR.tar.gz INSTALLER.sh"
  echo ""
  echo "Creates an FDS/Smokeview installer sh script. "
  echo ""
  echo "  ostype - OSX or LINUX"
  echo "  platform - default platform (ia32, intel64 or ib64)"
  echo "  FDS.tar.gz - compressed tar file containing FDS distribution"
  echo "  INSTALLER.sh - .sh script containing self-extracting installer"
  echo
  exit
fi

ostype=$1
platformin=$2
FDS_TAR=$3
INSTALLER=$4

LDLIBPATH=LD_LIBRARY_PATH
if [ "$ostype" == "OSX" ]
then
LDLIBPATH=DYLD_LIBRARY_PATH
fi

cat << EOF > $INSTALLER
#!/bin/bash
FDS_root=~/FDS/FDS6

# record the name of this script and the name of the directory 
# it will be run in

THIS=\`pwd\`/\$0
THISDIR=\`pwd\`

CSHFDS=/tmp/cshrc_fds.\$\$
BASHFDS=/tmp/bashrc_fds.\$\$

# Find the beginning of the included FDS tar file so that it can be
# subsequently un-tar'd
 
SKIP=\`awk '/^__TARFILE_FOLLOWS__/ { print NR + 1; exit 0; }' \$0\`

# extract tar.gz file from this script if option 1 is 'extract'

if [ \$# -eq 1 ]
then
option=\$1
if [ "\$option" == "extract" ]
then
name=\$0
THAT=\${name%.*}.tar.gz
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
echo Extracting the compressed tar file contained in this installer to \$THAT
tail -n +\$SKIP \$THIS > \$THAT
exit 0
fi
fi

#
# get FDS root directory
#
echo ""
echo "*** FDS 6 Installer ***"
echo ""
echo "Enter FDS installation directory (default: \$FDS_root)"
read answer
if [ "\$answer" != "" ]; then
FDS_root=\$answer
fi
 
# make the FDS root directory
 
if [ ! -d \$FDS_root ]
then
echo "creating directory \$FDS_root"
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

# Setting PATH and LD_LIBRARY_PATH environment
# variables for use by FDS
setenv FDSBINDIR \`pwd\`/bin
if ( "\\\$platform" == "intel64" ) then
setenv $LDLIBPATH \\\`pwd\\\`/bin/LIB64
endif

if ( "\\\$platform" == "ia32" ) then
setenv $LDLIBPATH \\\`pwd\\\`/bin/LIB32
endif

# add FDS bin to path

set path=(\$FDSBINDIR \\\$path)

CSHRC

cat << BASH > \$BASHFDS
#/bin/bash
platform=\\\$1

# define FDS bin directory location

export FDSBINDIR=\`pwd\`/bin

# define openmpi library locations:
#   32/64 bit gigabit ethernet
#   64 bit infiniband

MPIDIST32=/shared/openmpi_32
MPIDIST64=/shared/openmpi_64
MPIDIST64IB=/shared/openmpi_64ib

# environment for 64 bit gigabit ethernet

if [ "\\\$platform" == "intel64" ]
then
export MPIDIST=\\\$MPIDIST64
FORTLIB=\\\$FDSBINDIR/LIB64
fi

# environment for 32 bit gigabit ethernet

if [ "\\\$platform" == "ia32" ]
then
export MPIDIST=\\\$MPIDIST32
FORTLIB=\\\$FDSBINDIR/LIB32
fi

# environment for 64 bit infiniband

if [ "\\\$platform" == "ib64" ]
then
export MPIDIST=\\\$MPIDIST64IB
FORTLIB=\\\$FDSBINDIR/LIB64
fi

# Update LD_LIBRARY_PATH and PATH variables

export $LDLIBPATH=\\\$MPIDIST/lib:\\\$FORTLIB
export PATH=\\\$FDSBINDIR:\\\$MPIDIST/bin:\\\$PATH

BASH

echo
echo Creating .bashrc_fds and .cshrc_fds startup files.
if [ -e ~/.cshrc_fds ]
then
echo Backing up .cshrc_fds
cp ~/.cshrc_fds ~/.cshrc_fds.BAK_\$\$
fi

if [ -e ~/.bashrc_fds ]
then
echo Backing up .bashrc_fds
cp ~/.bashrc_fds ~/.bashrc_fds.BAK_\$\$
fi

cp \$BASHFDS ~/.bashrc_fds
rm \$BASHFDS

cp \$CSHFDS ~/.cshrc_fds
rm \$CSHFDS

cd \$THISDIR
nlines=\$(grep bashrc_fds ~/.bashrc | wc -l)
if [ \$nlines -eq 0 ]
then
echo Updating .bashrc
echo source .bashrc_fds $platformin
echo source .bashrc_fds $platformin >> ~/.bashrc
fi
nlines=\$(grep cshrc_fds ~/.cshrc | wc -l)
if [ \$nlines -eq 0 ]
then
echo Updating .cshrc
echo source .cshrc_fds $platformin >> ~/.cshrc
fi
echo ""
echo "Installation complete."
exit 0


__TARFILE_FOLLOWS__
EOF
chmod +x $INSTALLER
cat $FDS_TAR >> $INSTALLER
echo "Installer created."
