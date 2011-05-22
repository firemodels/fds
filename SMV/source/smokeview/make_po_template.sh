#!/bin/bash

MAKEPO=../../makepo/INTEL_LINUX_64/makepo_linux_64
cat *.c *.cpp | $MAKEPO | sort -u | $MAKEPO -a > smokeview_template.po
