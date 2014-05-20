@echo off
call setopts 
call makecase bench8a
call makecase bench8b
call makecase bench8c
call makecase bench8d
call makecase bench8e
call makecase bench8f
call makecase bench8g
call makecase bench8h
%bg% %fdsrun% bench8a.fds
%bg% %fdsrun% bench8b.fds
%bg% %fdsrun% bench8c.fds
%bg% %fdsrun% bench8d.fds
%bg% %fdsrun% bench8e.fds
%bg% %fdsrun% bench8f.fds
%bg% %fdsrun% bench8g.fds
%bg% %fdsrun% bench8h.fds