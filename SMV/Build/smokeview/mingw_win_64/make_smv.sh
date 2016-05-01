#!/bin/bash
rm -rf *.o *.mod
eval make -f ../Makefile mingw_win_64
