goto eof

:is_file_installed
set program=%1
%program% 1>> %temp%\file_exist.txt 2>&1
type %temp%\file_exist.txt | find /i /c "not recognized" > %temp%\file_exist_count.txt
set /p nothave=<%temp%\file_exist_count.txt
exit /b 0

:eof
