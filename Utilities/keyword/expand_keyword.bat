@echo off

:: expand keyword found in file using newvalue
::   (note that the ~ in %~2 removes surrounding quotes (allowing imbedded blanks in newvalue)

if "%1" == "" (
  echo usage: expand_keyword keyword newvalue file
  echo        in file, convert all occurrences of $keyword: .... $ to
  echo        $keyword: newvalue $
  goto eof
)

set keyword=%1
set newvalue=%~2
set file=%3

if NOT exist %file% (
  exit 1 /b
)
sed -e "s/$%keyword%:.*\$/$%keyword%: %newvalue% $/g" %file% | sed "s/$/\r/" > %temp%\temp.txt
copy %temp%\temp.txt %file% 1> Nul 2>&1

:eof