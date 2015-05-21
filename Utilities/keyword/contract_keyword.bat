@echo off

:: contract keyword located in file

if "%1" == "" (
  echo usage: contract_keyword keyword file
  echo        change all occurences of $keyword: .... $ to
  echo        $keyword: unknown $ in file
  goto eof
)

set keyword=%1
set file=%2

if NOT exist %file% (
  exit 1 /b
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

sed -e "s/$%keyword%:.*\$/$%keyword%: unknown $/g" %file% | sed "s/$/\r/" > %temp%\temp2.txt
copy %temp%\temp2.txt %file% 1> Nul 2>&1

:eof