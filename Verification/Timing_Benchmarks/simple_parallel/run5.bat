@echo off
set fdsrun=fds
set bg=call background -d 0 -u 100 
call makecase bench5a
call makecase bench5b
call makecase bench5c
call makecase bench5d
call makecase bench5e
%bg% %fdsrun% bench5a.fds
%bg% %fdsrun% bench5b.fds
%bg% %fdsrun% bench5c.fds
%bg% %fdsrun% bench5d.fds
%bg% %fdsrun% bench5e.fds
