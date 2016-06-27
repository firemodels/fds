@echo off

echo Creating figures for the FDS+Evac User's and Verification guide

set evac="%CD%\..\Evacuation"
set evacug="%CD%\..\..\Manuals\Evac_5_Guide"

cd %evacug%

erase SCRIPT_FIGURES\*.png
erase SCRIPT_FIGURES\*.help
erase SCRIPT_FIGURES\*.version

smokeview -help > SCRIPT_FIGURES\smokeview.help
smokeview -version > SCRIPT_FIGURES\smokeview.version

cd %evac%
smokeview -runscript evac_smv_testcase1 > evac_smv_testcase1.log 2>&1
smokeview -runscript evac_smv_testcase2 > evac_smv_testcase2.log 2>&1

echo SMV Evacuation pictures done
pause Press any key...
