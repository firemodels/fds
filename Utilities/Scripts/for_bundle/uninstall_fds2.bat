goto eof

:is_cfast_installed
cfast 1> %temp%\file_exist.txt 2>&1
type %temp%\file_exist.txt | find /i /c "not recognized" > %temp%\file_exist_count.txt
set /p nothave=<%temp%\file_exist_count.txt
set cfastinstalled=1
if %nothave% == 1 (
  set cfastinstalled=0
)
exit /b 0

:eof
