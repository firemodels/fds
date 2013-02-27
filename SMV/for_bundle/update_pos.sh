#!/bin/bash
#
MERGEPO=../mergepo/intel_linux_64_db/mergepo_linux_64_db

for pofile in smokeview_??.po
do
  echo updating $pofile
  $MERGEPO -c $pofile smokeview_template.po > $pofile.revised
done
