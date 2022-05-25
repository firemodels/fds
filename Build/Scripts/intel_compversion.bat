@echo off
setlocal enabledelayedexpansion
set prog=%1

call :is_file_installed %prog% || exit /b 1

%prog% /version > f_version.txt 2>&1
set /p vers=<f_version.txt

set nargs=0
for %%a in ("%vers: =" "%") do (
   set vector2[!nargs!]=%%~a
   set /A nargs+=1
   set vector[!nargs!]=%%a
)

FOR /L %%i IN (0,1,%nargs%) DO (
  if !vector[%%i]! == "Version"     echo "Intel %prog% !vector2[%%i]!"
)

erase f_version.txt
goto eof

:: -------------------------------------------------------------
:is_file_installed
:: -------------------------------------------------------------

  set program=%1
  %program% 1> output.txt 2>&1
  type output.txt | find /i /c "not recognized" > output_count.txt
  set /p nothave=<output_count.txt
  if %nothave% == 1 (
    echo unknown
    erase output.txt output_count.txt
    exit /b 1
  )
  erase output.txt output_count.txt
  exit /b 0

:eof
