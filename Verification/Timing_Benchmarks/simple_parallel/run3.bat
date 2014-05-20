@echo off
set fdsrun=fds
set bg=call background -d 0 -u 100 
call makecase bench3a
call makecase bench3b
call makecase bench3c
%bg% %fdsrun% bench3a.fds
%bg% %fdsrun% bench3b.fds
%bg% %fdsrun% bench3c.fds