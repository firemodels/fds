@echo off
set filein=%1
set fileout=%filein%.md5

call :is_file_installed certutil || exit /b 1
if exist %filein% (
  certutil -hashfile %filein% md5 > %fileout%
)
goto eof

:: -------------------------------------------------------------
:is_file_installed
:: -------------------------------------------------------------

  set program=%1
  %program% --help 1>> file_exist.txt 2>&1
  type file_exist.txt | find /i /c "not recognized" > file_exist_count.txt
  set /p nothave=<file_exist_count.txt
  erase file_exist.txt
  erase file_exist_count.txt
  if %nothave% == 1 (
    exit /b 1
  )
  exit /b 0

:eof