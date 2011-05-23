#!/bin/csh -f
#
set here=`pwd`
set SOURCEDIR=../source/smokeview
set MAKEPO=../makepo/INTEL_LINUX_64/makepo_linux_64

echo updating smokeview_template.po
cat $SOURCEDIR/*.c $SOURCEDIR/*.cpp $SOURCEDIR/*.h | $MAKEPO | sort -u | $MAKEPO -a > smokeview_template.po

foreach lang (de fr it se)
echo
echo creating $lang translation
set infile=smokeview_template.po
set outfile=smokeview_$lang.po
curl -F pofile=@$infile \
        -F language=$lang \
        -F output=pofile  \
        http://bryanklein.com/trans/index.php \
        --output $outfile
end
