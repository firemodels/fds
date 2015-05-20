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

sed -e "s/$%keyword%:.*\$/$%keyword%: unknown $/g" %file% | sed "s/$/\r/" > %temp%\temp2.txt
copy %temp%\temp2.txt %file% 1> Nul 2>&1

:eof