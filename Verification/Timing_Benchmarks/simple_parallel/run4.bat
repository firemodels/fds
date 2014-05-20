@echo off
set fdsrun=fds
set bg=call background -d 0 -u 100 
call makecase bench4a
call makecase bench4b
call makecase bench4c
call makecase bench4d
%bg% %fdsrun% bench4a.fds
%bg% %fdsrun% bench4b.fds
%bg% %fdsrun% bench4c.fds
%bg% %fdsrun% bench4d.fds