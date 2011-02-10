#!/bin/bash
EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: make_installer.sh FORTLIB FDS_TAR INSTALLER"
  echo ""
  echo "Creates an FDS/Smokeview installer sh script. "
  echo ""
  echo "  FDS_TAR - compressed tar file contining distribution"
  echo "  INSTALLER - .sh script containing self-extracting installer"
  echo
  exit
fi

FDS_TAR=$1
INSTALLER=$2

cat << EOF > $INSTALLER
#!/bin/bash
FDS_root=~/FDS/FDS6

# record the name of this script and the name of the directory 
# it will be run in

THIS=\`pwd\`/\$0
THISDIR=\`pwd\`

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
    read -p "Do you wish to overwrite contents? (yes/no)" yn
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

echo "Copying FDS installation files to"  \$FDS_root
cd \$FDS_root
tail -n +\$SKIP \$THIS | tar -xz
echo "Copy complete."

echo "Setting path and LD_LIBRARY_PATH environment variables"
echo "in .bashrc_fds and .cshrc_fds"

echo "" > ~/.cshrc_fds
echo "#/bin/csh -f" >> ~/.cshrc_fds
echo "set platform=\\\$1" >> ~/.cshrc_fds
echo "" >> ~/.cshrc_fds
echo "# Setting PATH and LD_LIBRARY_PATH environment" >> ~/.cshrc_fds
echo "# variables for use by FDS" >> ~/.cshrc_fds
echo "" >> ~/.cshrc_fds

echo "#/bin/bash" > ~/.bashrc_fds
echo "platform=\\\$1" >> ~/.bashrc_fds
echo "" >> ~/.bashrc_fds
echo "# Setting PATH and LD_LIBRARY_PATH environment" >> ~/.bashrc_fds
echo "# variables for use by FDS" >> ~/.bashrc_fds
echo "" >> ~/.bashrc_fds
echo "setenv FDSBINDIR \`pwd\`/bin" >> ~/.cshrc_fds
echo "export FDSBINDIR=\`pwd\`/bin" >> ~/.bashrc_fds
if [ "a\$LD_LIBRARY_PATH" = "a" ]
then
echo "if ( \"\\\$platform\" == \"intel64\" ) then" >> ~/.cshrc_fds
echo "setenv LD_LIBRARY_PATH \`pwd\`/bin/LIB64" >> ~/.cshrc_fds
echo "endif" >> ~/.cshrc_fds

echo "if ( \"\\\$platform\" == \"ia32\" ) then" >> ~/.cshrc_fds
echo "setenv LD_LIBRARY_PATH \`pwd\`/bin/LIB32" >> ~/.cshrc_fds
echo "endif" >> ~/.cshrc_fds

echo "if [ \"\\\$platform\" == \"intel64\" ]" >> ~/.bashrc_fds
echo "then" >> ~/.bashrc_fds
echo "export LD_LIBRARY_PATH=\`pwd\`/bin/LIB64" >> ~/.bashrc_fds
echo "fi" >> ~/.bashrc_fds

echo "if [ \"\\\$platform\" == \"ia32\" ]" >> ~/.bashrc_fds
echo "then" >> ~/.bashrc_fds
echo "export LD_LIBRARY_PATH=\`pwd\`/bin/LIB32" >> ~/.bashrc_fds
echo "fi" >> ~/.bashrc_fds
else
echo "setenv LD_LIBRARY_PATH \`pwd\`/bin/$FORTLIB:\\\$LD_LIBRARY_PATH" >> ~/.cshrc_fds

echo "if [ \"\\\$platform\" == \"intel64\" ]" >> ~/.bashrc_fds
echo "then" >> ~/.bashrc_fds
echo "export LD_LIBRARY_PATH=\`pwd\`/bin/LIB64:\\\$LD_LIBRARY_PATH" >> ~/.bashrc_fds
echo "fi" >> ~/.bashrc_fds

echo "if [ \"\\\$platform\" == \"ia32\" ]" >> ~/.bashrc_fds
echo "then" >> ~/.bashrc_fds
echo "export LD_LIBRARY_PATH=\`pwd\`/bin/LIB32:\\\$LD_LIBRARY_PATH" >> ~/.bashrc_fds
echo "fi" >> ~/.bashrc_fds
fi

# add FDS bin to path

echo "set path=(\`pwd\`/bin \\\$path)" >> ~/.cshrc_fds
echo "export PATH=\\\$PATH:\`pwd\`/bin" >> ~/.bashrc_fds

cd \$THISDIR
echo ""
echo "If you use the csh or tcsh shell,"
echo "add the following line to your .cshrc file:"
echo "source ~/.cshrc_fds"
echo ""
echo "If you use the bash or sh shell,"
echo "add the following line to your .bashrc file:"
echo "source ~/.bashrc_fds"
echo ""
echo "Installation complete."
exit 0


__TARFILE_FOLLOWS__
EOF
chmod +x $INSTALLER
cat $FDS_TAR >> $INSTALLER
echo "Installer created."
