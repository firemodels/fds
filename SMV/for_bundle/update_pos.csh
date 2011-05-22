#!/bin/csh -f
#
set here=`pwd`
set MERGEPO=../mergepo/INTEL_LINUX_64/mergepo_linux_64
set MAKEPO=../makepo/INTEL_LINUX_64/makepo_linux_64
set SOURCEDIR=../source/smokeview

echo updating smokeview_template.po
cat $SOURCEDIR/*.c $SOURCEDIR/*.cpp | $MAKEPO | sort -u | $MAKEPO -a > smokeview_template.po

foreach pofile (smokeview_??.po)
echo updating $pofile
$MERGEPO $pofile smokeview_template.po > xxx.po
mv xxx.po $pofile
end
