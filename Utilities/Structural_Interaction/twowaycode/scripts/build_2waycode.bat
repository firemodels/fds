
:: Firebot variables
set FIREBOT_DIR=%FDS_SVNROOT%/Utilities/Structural_Interaction/twowaycode/scripts
set TWOWAY_DIR=%FDS_SVNROOT%/Utilities/Structural_Interaction/twowaycode
set OUTPUT_DIR=%TWOWAY_DIR%/scripts/output
set ERROR_LOG=%OUTPUT_DIR%/errors
set WARNING_LOG=%OUTPUT_DIR%/warnings

:: Update Repository
set SVN_REVISION="%1"
echo %SVN_REVISION%
cd %FDS_SVNROOT%/FDS_Source
svn up -r %SVN_REVISION%

:: Clean outputs
cd %FIREBOT_DIR%
echo Y | del output

:: Create output dir
mkdir output

:: Clean previous results
cd %TWOWAY_DIR%/examples/simply_beam
del *.csv 

:: Compile fds_win_64
cd %FDS_SVNROOT%/FDS_Compilation/intel_win_64
echo Y | make_fds.bat

:: Compile twowaycode_win_64
cd %TWOWAY_DIR%/intel_win_64
echo Y | make_twowaycode.bat

:: Print the FDS revision number on User Guide
cd %TWOWAY_DIR%
sed -i "s:.*SVN Repository Revision.*:SVN Repository Revision %SVN_REVISION%:" twowaycode_user_guide.tex

:: Print the FDS revision number on python scripts
cd %FIREBOT_DIR%
sed -i "s:.*SVN=.*:SVN='%SVN_REVISION%':" generate_plots.py

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
svn revert twowaycode_user_guide.tex
del sed*

:: Revert the FDS revision number on python scripts
cd %FIREBOT_DIR%
svn revert generate_plots.py
del sed*