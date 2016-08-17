@echo off

call :is_file_installed magick|| exit /b 1
echo             found jpeg conversion program magick

for %%f in (*.jp2) do (
  echo.
  echo converting %%~nf.jp2 from jpeg 2000 to jpeg
  magick %%~nf.jp2  %%~nf.jpg
)

goto eof

:: -------------------------------------------------------------
:is_file_installed
:: -------------------------------------------------------------

  set program=%1
  %program% --help 1> %temp%\file_exist.txt 2>&1
  type %temp%\file_exist.txt | find /i /c "not recognized" > %temp%\file_exist_count.txt
  set /p nothave=<%temp%\file_exist_count.txt
  if %nothave% == 1 (
    echo ***Fatal error: %program% not present
    echo                 jpeg2 to jpeg image file conversion aborted
    exit /b 1
  )
  exit /b 0

:eof
