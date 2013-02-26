#!/bin/bash
#
MERGEPO=../mergepo/gcc_linux_64/mergepo_linux_64

for pofile in smokeview_??.po
do
  echo updating $pofile
  $MERGEPO -c $pofile smokeview_template.po > xxx.po
  mv xxx.po $pofile.revised
done
