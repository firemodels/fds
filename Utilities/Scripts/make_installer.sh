#!/bin/bash
EXPECTED_ARGS=3

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: make_installer.sh FORTLIB FDS_TAR INSTALLER"
  echo ""
  echo "Creates an FDS/Smokeview installer sh script. "
  echo ""
  echo "  FORTLIB - directory containing run-time libraries"
  echo "  FDS_TAR - compressed tar file contining distribution"
  echo "  INSTALLER - .sh script containing self-extracting installer"
  echo
  exit
fi

FORTLIB=$1
FDS_TAR=$2
INSTALLER=$3

cat << EOF > $INSTALLER
#!/bin/bash
FDS_root=~/FDS/FDS6

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

# Find the beginning of the included FDS tar file
# so it can be subsequently un-tar'd
 
SKIP=\`awk '/^__TARFILE_FOLLOWS__/ { print NR + 1; exit 0; }' \$0\`

# record the name of this script and the name of the directory it will be run in

THIS=\`pwd\`/\$0
THISDIR=\`pwd\`
echo "Preliminaries complete."
while true; do
    read -p "Do you wish to proceed with the installation? (yes/no)" yn
    case \$yn in
        [Yy]* ) break;;
        [Nn]* ) echo "Installation cancelled";exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

# now copy installation files into the FDS_root directory

echo "Copying FDS installation files to"  \$FDS_root
cd \$FDS_root
tail -n +\$SKIP \$THIS | tar -xzv
echo "Copy complete."

echo "Setting path and LD_LIBRARY_PATH environment variables."

echo "" > ~/.cshrc_fds
echo "# Setting PATH and LD_LIBRARY_PATH environment" >> ~/.cshrc_fds
echo "# variables for use by FDS" >> ~/.cshrc_fds
echo "" >> ~/.cshrc_fds
echo "" > ~/.bashrc_fds
echo "# Setting PATH and LD_LIBRARY_PATH environment" >> ~/.bashrc_fds
echo "# variables for use by FDS" >> ~/.bashrc_fds
echo "" >> ~/.bashrc_fds
if [ "a\$LD_LIBRARY_PATH" = "a" ]
then
echo "setenv LD_LIBRARY_PATH \`pwd\`/bin/$FORTLIB" >> ~/.cshrc_fds
echo "export LD_LIBRARY_PATH=\`pwd\`/bin/$FORTLIB" >> ~/.bashrc_fds
else
echo "setenv LD_LIBRARY_PATH \`pwd\`/bin/$FORTLIB:\\\$LD_LIBRARY_PATH" >> ~/.cshrc_fds
echo "export LD_LIBRARY_PATH=\`pwd\`/bin/$FORTLIB:\\\$LD_LIBRARY_PATH" >> ~/.bashrc_fds
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
