@echo off
set fdsrun=fds
set bg=call background -d 0 -u 100 
call makecase bench1a
%bg% %fdsrun% bench1a.fds