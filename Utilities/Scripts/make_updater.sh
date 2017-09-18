#!/bin/bash
EXPECTED_ARGS=5

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: make_updater.sh ostype SMV_TAR.tar.gz INSTALLER.sh"
  echo ""
  echo "Creates a Smokeview updater sh script. "
  echo ""
  echo "  ostype - OSX or LINUX"
  echo "  SMV.tar.gz - compressed tar file containing Smokeview update"
  echo "  INSTALLER.sh - .sh script containing self-extracting updater"
  echo "  updatedir - default update directory"
  echo
  exit
fi

ostype=$1
revision=$2
SMV_TAR=$3
INSTALLER=$4
INSTALLDIR=$5

ossize=64

cat << EOF > $INSTALLER
#!/bin/bash

OVERRIDE=\$1
echo ""
echo "Updating Smokeview to $revision for $ossize bit $ostype"
echo ""
echo "Options:"
echo "  1) Press <Enter> to begin update"
echo "  2) Type \"extract\" to copy the update files to $SMV_TAR"

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

THISSCRIPT=\`ABSPATH \$0\`
THISDIR=\`pwd\`

#--- Find the beginning of the included SMV tar file so that it can be
# subsequently un-tar'd
 
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
  THAT=$SMV_TAR
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
  echo Extracting the file embedded in this updater to \$THAT
  tail -n +\$SKIP \$THISSCRIPT > \$THAT
  exit 0
fi

#--- get FDS root directory

echo ""
echo "Where would you like to put the Smokeview update?"
EOF
  if [ "$ostype" == "OSX" ]; then
cat << EOF >> $INSTALLER
    echo "  Press 1 to update in /Applications/$INSTALLDIR"
    echo "  Press 2 to update in \$HOME/$INSTALLDIR"
EOF
  else
cat << EOF >> $INSTALLER
    echo "  Press 1 to update in \$HOME/$INSTALLDIR"
    echo "  Press 2 to update in /opt/$INSTALLDIR"
    echo "  Press 3 to update in /usr/local/bin/$INSTALLDIR"
EOF
  fi
cat << EOF >> $INSTALLER
echo "  Enter directory path to update elsewhere"

if [ "\$OVERRIDE" == "y" ] 
then
  answer="1"
else
  read answer
fi
EOF

if [ "$ostype" == "OSX" ]; then
cat << EOF >> $INSTALLER
  if [[ "\$answer" == "1" || "\$answer" == "" ]]; then
    eval SMV_root=/Applications/$INSTALLDIR
  elif [ "\$answer" == "2" ]; then
    eval SMV_root=\$HOME/$INSTALLDIR
  else
    eval SMV_root=\$answer
  fi
EOF
else
cat << EOF >> $INSTALLER
  if [[ "\$answer" == "1" || "\$answer" == "" ]]; then
    eval SMV_root=\$HOME/$INSTALLDIR
  elif [ "\$answer" == "2" ]; then
    SMV_root=/opt/$INSTALLDIR
  elif [ "\$answer" == "3" ]; then
    SMV_root=/usr/local/bin/$INSTALLDIR
  else
    eval SMV_root=\$answer
  fi
EOF
fi
cat << EOF >> $INSTALLER

if [ ! -d \$SMV_root ]
then
   echo FDS/Smokeview installation directory \$SMV_root does not exist.
   echo Install FDS and Smokeview before proceeding with update.
   echo Update cancelled.
   exit
fi
touch \$SMV_root/test_write >& /dev/null
if [ ! -e \$SMV_root/test_write ]
then
   echo Unable to write to installation directory.
   echo Update cancelled
   exit
fi
rm \$SMV_root/test_write

#--- do we want to proceed

while true; do
    echo ""
    echo "Update directory: \$SMV_root"
    if [ "\$OVERRIDE" == "y" ] 
    then
      yn="y"
    else
      read -p "Do you wish to begin the update? (yes/no) " yn
    fi
    case \$yn in
        [Yy]* ) break;;
        [Nn]* ) echo "Update cancelled";exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
 
echo "Copying updated Smokeview files to"  \$SMV_root
cd \$SMV_root
tail -n +\$SKIP \$THISSCRIPT | tar -xz
echo "Update complete."
exit 0


__TARFILE_FOLLOWS__
EOF
chmod +x $INSTALLER
cat $SMV_TAR >> $INSTALLER
echo "Updater created."
