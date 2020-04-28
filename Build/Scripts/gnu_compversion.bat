@echo off
call :is_file_installed head|| exit /b 1
call :is_file_installed sed || exit /b 1
call :is_file_installed gawk|| exit /b 1
call :is_file_installed gfortran || exit /b 1

gfortran --version > f_version.txt 2>&1
head -1 f_version.txt | sed "s/^.*\(Version.*\).*$/\1/" | gawk "{print $4}" > vers.out
set /p vers=<vers.out
echo "Gnu gfortran %vers%"
erase f_version.txt vers.out
goto eof

:: -------------------------------------------------------------
:is_file_installed
:: -------------------------------------------------------------

  set program=%1
  %program% --help 1> output.txt 2>&1
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


