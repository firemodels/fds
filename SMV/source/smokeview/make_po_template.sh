#!/bin/bash

MAKEPO=../../makepo/intel_linux_64/makepo_linux_64
#cat *.c *.cpp | $MAKEPO | sort -u | $MAKEPO -a > smokeview_template.po
cat *.c *.cpp | $MAKEPO | sort -u  > smokeview_template.po
