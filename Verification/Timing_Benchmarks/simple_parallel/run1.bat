@echo off
call setopts 
call makecase bench1a
%bg% %fdsrun% bench1a.fds