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

set fullfile=%dir%\%file%

if NOT exist %fullfile% (
  exit /b 1
)

call "%bindir%\get_repo_properties" %dir%

call "%bindir%\expand_keyword" Revision %revision% %fullfile%
call "%bindir%\expand_keyword" RevisionDate "%revision_date% %revision_time%" %fullfile%
call "%bindir%\expand_keyword" CompileDate "%build_date% %build_time%" %fullfile%

:eof
