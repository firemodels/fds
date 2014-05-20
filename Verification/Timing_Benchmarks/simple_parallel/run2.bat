@echo off
call setopts 
call makecase bench2a
call makecase bench2b
%bg% %fdsrun% bench2a.fds
%bg% %fdsrun% bench2b.fds