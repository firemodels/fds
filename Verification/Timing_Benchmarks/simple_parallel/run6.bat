@echo off
call setopts 
call makecase bench6a
call makecase bench6b
call makecase bench6c
call makecase bench6d
call makecase bench6e
call makecase bench6f
%bg% %fdsrun% bench6a.fds
%bg% %fdsrun% bench6b.fds
%bg% %fdsrun% bench6c.fds
%bg% %fdsrun% bench6d.fds
%bg% %fdsrun% bench6e.fds
%bg% %fdsrun% bench6f.fds