
:: Firebot variables
set FIREBOT_DIR=%FDS_GITROOT%/Utilities/Structural_Interaction/twowaycode/scripts
set TWOWAY_DIR=%FDS_GITROOT%/Utilities/Structural_Interaction/twowaycode
set OUTPUT_DIR=%TWOWAY_DIR%/scripts/output
set ERROR_LOG=%OUTPUT_DIR%/errors
set WARNING_LOG=%OUTPUT_DIR%/warnings

:: Get Git Hash
for /f "delims=" %%i in ('git describe --long --dirty') do set GIT_HASH=%%i 
echo %GIT_HASH%

:: Clean outputs
cd %FIREBOT_DIR%
echo Y | rmdir output

:: Create output dir
mkdir output

:: Clean previous results
cd %TWOWAY_DIR%/examples/simply_beam
del *.csv 
cd %FIREBOT_DIR%/SCRIPT_FIGURES
del *.pdf

:: Compile fds_win_64
cd %FDS_GITROOT%/FDS_Compilation/intel_win_64
echo Y | make_fds.bat

:: Compile twowaycode_win_64
cd %TWOWAY_DIR%/intel_win_64
echo Y | make_twowaycode.bat

:: Print the FDS revision number on User Guide
cd %TWOWAY_DIR%
sed -i "s:.*Git Hash.*:\\path{%GIT_HASH%}:" twowaycode_user_guide.tex

:: Print the FDS revision number on python scripts
cd %FIREBOT_DIR%
sed -i "s:.*GIT=.*:GIT='%GIT_HASH%':" generate_plots.py

:: Run Verification Cases
cd %TWOWAY_DIR%/examples/simply_beam
simply_beam.bat

:: Run python scripts
cd %FIREBOT_DIR%
python generate_plots.py

:: Compile User Guide
cd %TWOWAY_DIR%
set TEXINPUTS=%TEXINPUTS%:.:../../../Manuals/LaTeX_Style_Files/
pdflatex twowaycode_user_guide.tex
bibtex twowaycode_user_guide
pdflatex twowaycode_user_guide.tex

:: Revert the FDS revision number on User Guide
cd %TWOWAY_DIR%
git checkout -- twowaycode_user_guide.tex
del sed*

:: Revert the FDS revision number on python scripts
cd %FIREBOT_DIR%
git checkout -- generate_plots.py
del sed*