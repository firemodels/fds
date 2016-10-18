@echo off

call BUILD_prog background
call BUILD_prog dem2fds
call BUILD_prog fds2ascii
call BUILD_prog smokediff
call BUILD_prog smokezip
call BUILD_prog wind2fds