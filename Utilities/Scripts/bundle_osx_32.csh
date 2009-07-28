#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
set fds_smvroot=~/$1
set scp_fds_smvroot=$1
set makedir=$scp_fds_smvroot/FDS_Compilation
set googledir=$fds_smvroot/Utilities/to_google
set bundledir=$2
set bundle_setup=$fds_smvroot/Utilities/Scripts/bundle_setup
set mandir=$fds_smvroot/Manuals/All_PDF_Files
set smvbindir=$scp_fds_smvroot/SMV_5/bin
set osxhost=devi1.nist.gov

cd $googledir
rm -rf $bundledir
mkdir $bundledir
mkdir $bundledir/bin
mkdir $bundledir/Documentation
mkdir $bundledir/Examples

echo Copying files
scp $osxhost\:$makedir/intel_osx_32/fds5_linux_osx_32 $bundledir/bin/.
scp $osxhost\:$makedir/mpi_intel_osx_32/fds5_mpi_intel_osx_32 $bundledir/bin/.
scp tiger.cfr.nist.gov\:$smvbindir/smv5_osx_32 $bundledir/bin/.
scp tiger.cfr.nist.gov\:$smvbindir/smokezip_osx $bundledir/bin/.
cp ~/$smvbindir/smokeview.ini $bundledir/bin/.

cp $bundle_setup/readme_docs.html $bundledir/Documentation/.
cp $mandir/FDS_5_User_Guide.pdf $bundledir/Documentation/.
cp $mandir/SMV_5_User_Guide.pdf $bundledir/Documentation/.

cp $bundle_setup/readme_examples.html $bundledir/Examples/.

echo Building archive
rm -rf $googledir/$bundledir.tar
rm -rf $googledir/$bundledir.tar.gz
tar cvf $googledir/$bundledir.tar $bundledir/.
echo Compressing archive
gzip    $googledir/$bundledir.tar
