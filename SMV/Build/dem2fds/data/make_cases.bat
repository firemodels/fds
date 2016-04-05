@echo off
set case=blodget
dem2fds -e %case% < %case%_elevs.csv > %case%.fds

set case=nist
dem2fds -e %case% < %case%_elevs.csv > %case%.fds

set case=sugarloaf
dem2fds -e %case% < %case%_elevs.csv > %case%.fds

set case=trails
dem2fds -e %case% < %case%_elevs.csv > %case%.fds