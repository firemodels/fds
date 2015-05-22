@echo off

:: expand keywords Revision, RevisionDate and CompileDate in file

if "%1" == "" (
  echo usage: expand_file dir file
  echo        expand all occurrences of the keywords Revision, 
  echo        RevisionDate and CompileDate in file using properties
  echo        of the directory dir
  goto eof
)

set bindir=%~p0
set dir=%1
set file=%2

if NOT exist %file% (
  goto eof
)
call "%bindir%\get_repo_properties" %dir%

call "%bindir%\expand_keyword" Revision %revision% %file%
call "%bindir%\expand_keyword" RevisionDate "%revision_date%  %revision_time%" %file%
call "%bindir%\expand_keyword" CompileDate "%build_date%  %build_time%" %file%

:eof
