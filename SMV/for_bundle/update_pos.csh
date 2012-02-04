#!/bin/csh -f
#
set here=`pwd`
set MERGEPO=../mergepo/intel_linux_64/mergepo_linux_64

foreach pofile (smokeview_??.po)
echo updating $pofile
$MERGEPO $pofile smokeview_template.po > xxx.po
mv xxx.po $pofile.revised
end
