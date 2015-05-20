@echo off

:: contract keywords Revision, RevisionDate and CompileDate in file

if "%1" == "" (
  echo usage: contract_file file
  echo        contract all occurrences of the keywords Revision, 
  echo        RevisionDate and CompileDate in file
  goto eof
)

set bindir=%~p0
set file=%1

if NOT exist %file% (
  exit /b 1
)

call "%bindir%\contract_keyword" Revision %file%
call "%bindir%\contract_keyword" RevisionDate %file%
call "%bindir%\contract_keyword" CompileDate %file%

:eof
