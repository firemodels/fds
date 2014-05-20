@echo off
set fdsrun=fds
set bg=call background -d 0 -u 100 
call makecase bench2a
call makecase bench2b
%bg% %fdsrun% bench2a.fds
%bg% %fdsrun% bench2b.fds