#!/bin/csh -f
#
set here=`pwd`
set SOURCEDIR=../source/smokeview
set MAKEPO=../makepo/INTEL_LINUX_64/makepo_linux_64

foreach lang (pl)
echo
echo creating $lang translation
set infile=smokeview_template.po
set outfile=smokeview_$lang_initial.po
curl -F pofile=@$infile \
        -F language=$lang \
        -F output=pofile  \
        http://bryanklein.com/trans/index.php \
        --output $outfile
end
