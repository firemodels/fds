@echo off
set fdsrun=fds
set bg=call background -d 0 -u 100 
call makecase bench7a
call makecase bench7b
call makecase bench7c
call makecase bench7d
call makecase bench7e
call makecase bench7f
call makecase bench7g
%bg% %fdsrun% bench7a.fds
%bg% %fdsrun% bench7b.fds
%bg% %fdsrun% bench7c.fds
%bg% %fdsrun% bench7d.fds
%bg% %fdsrun% bench7e.fds
%bg% %fdsrun% bench7f.fds
%bg% %fdsrun% bench7g.fds