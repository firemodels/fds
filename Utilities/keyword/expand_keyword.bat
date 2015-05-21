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
  goto eof
)

:: looking for sed

set temp1="%temp%\temp.txt"
set temp1c="%temp%\tempc.txt"

sed 1> %temp1% 2>&1
type %temp1% | find /i /c "not recognized" > %temp1c%
set flag=
set /p flag=<%temp1c%
if %flag% == 1 (
  goto eof
)

sed -e "s/$%keyword%:.*\$/$%keyword%: %newvalue% $/g" %file% | sed "s/$/\r/" > %temp%\temp.txt
copy %temp%\temp.txt %file% 1> Nul 2>&1

:eof