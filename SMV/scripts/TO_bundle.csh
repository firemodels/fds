#!/bin/csh -f

# script to copy Smokeview files to a "holding" area in 
# preparation for a bundle

set FOR_BUNDLE=~/FDS-SMV/SMV/for_bundle
set FROM_LINUX_32=$FOR_BUNDLE/smv_5.3.10_3188_linux
set FROM_OSX_32=$FOR_BUNDLE/smv_5.3.10_3188_osx

# --------------------------------------
# shouldn't need to edit lines below unless the files 
# in the distribution change

set MANUALS=~/FDS-SMV/Manuals/All_PDF_Files

set BRYAN_SHARE=~/FDS_SV_SHARE/NIST
set TO_LINUX_32=$BRYAN_SHARE/Linux/NIST/Smokeview
set TO_WIN_32=$BRYAN_SHARE/Windows/Smokeview
set TO_OSX_32=$BRYAN_SHARE/OSX/NIST/Smokeview

# -------------  Copy Linux files ---------------

echo copying LINUX_32 Smokeview files 
echo   from dir: $FROM_LINUX_32
echo     to dir: $TO_LINUX_32
cp $FOR_BUNDLE/devices.svo $TO_LINUX_32/.
cp $FOR_BUNDLE/readme.html $TO_LINUX_32/.
cp $FOR_BUNDLE/smokeview.ini $TO_LINUX_32/.
cp $MANUALS/SMV_User_Guide.pdf $TO_LINUX_32/Documentation/.
cp $FROM_LINUX_32/smokezip_linux $TO_LINUX_32/.
cp $FROM_LINUX_32/smv5_linux_32 $TO_LINUX_32/smv5_linux
cp ~/bin/libimf.so $TO_LINUX_32/.
cp $FOR_BUNDLE/readme_linux.first $TO_LINUX_32/.

# -------------  Copy Windows files ---------------

echo
echo copying WIN_32 Smokeview files
echo   from dir: $FOR_BUNDLE
echo     to dir: $TO_WIN_32
cp $FOR_BUNDLE/devices.svo $TO_WIN_32/.
cp $FOR_BUNDLE/glew32.dll $TO_WIN_32/.
cp $FOR_BUNDLE/pthreadVC.dll $TO_WIN_32/.
cp $FOR_BUNDLE/readme.html $TO_WIN_32/.
cp $FOR_BUNDLE/smokeview.ini $TO_WIN_32/.
cp $MANUALS/SMV_User_Guide.pdf $TO_LINUX_32/Documentation/.
cp $FOR_BUNDLE/smokeview_release.exe $TO_WIN_32/smokeview.exe
cp $FOR_BUNDLE/smokezip_release.exe $TO_WIN_32/smokezip.exe

# -------------  Copy OSX files ---------------

echo
echo copying OSX_32 Smokeview files 
echo   from dir: $FROM_OSX_32
echo     to dir: $TO_OSX_32
cp $FOR_BUNDLE/devices.svo $TO_OSX_32/.
cp $FOR_BUNDLE/readme.html $TO_OSX_32/.
cp $FOR_BUNDLE/smokeview.ini $TO_OSX_32/.
cp $MANUALS/SMV_User_Guide.pdf $TO_OSX_32/Documentation/.
cp $FROM_OSX_32/smokezip_osx $TO_OSX_32/.
cp $FROM_OSX_32/smv5_osx_32 $TO_OSX_32/smv5_osx
