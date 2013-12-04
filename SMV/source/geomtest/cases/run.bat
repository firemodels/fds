@echo  off
set BASEDIR=%CD%

cd %BASEDIR%\..\..\..\..
set SVNROOT=%CD%

set GEOMTEST=%SVNROOT%\SMV\source\geomtest\intel_win_64\geomtest
set CASEDIR=%SVNROOT%\SMV\source\geomtest\cases
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

cd %CASEDIR%

%SH2BAT% geom_cases.sh geom_cases.bat


call geom_cases.bat

