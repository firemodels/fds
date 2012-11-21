@echo off
echo copying fds_smv_env.bat to %homedrive%%homepath%
echo press any key to proceed or ^<CTRL^> C to abort
pause > NUL
copy /-Y ..\..\SMV\scripts\fds_smv_env.bat "%homedrive%%homepath%"
echo copy complete.  press any key to finish.
pause > NUL
