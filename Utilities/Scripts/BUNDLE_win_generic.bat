@echo off

set fdsversion=%fds_edition%
set smvversion=SMV6

set SVNROOT=%svn_root%
set fdsdir=%svn_root%\FDS_Compilation\intel_win_%platform%
set fdsmpidir=%svn_root%\FDS_Compilation\mpi_intel_win_%platform%
set basename=FDS_%fds_version%-SMV_%smv_version%_win%platform%

set in_pdf=%userprofile%\FDS_Guides
set in_intel_dll=%userprofile%\FIRE-LOCAL\LIBS\WINDOWS
set in_fds2ascii=%svn_root%\Utilities\fds2ascii
set in_setpath=%svn_root%\SMV\Build\set_path\intel_win_64
set in_shortcut=%userprofile%\Fire-LOCAL\repo_exes
set in_smokediff=%svn_root%\SMV\Build\smokediff
set in_smokezip=%svn_root%\SMV\Build\smokezip
set in_dem2fds=%svn_root%\SMV\Build\dem2fds
set in_smvscriptdir=%svn_root%\SMV\scripts
set in_wind2fds=%svn_root%\SMV\Build\wind2fds
set in_testmpi=%svn_root%\Utilities\test_mpi\impi_intel_win
set in_background=%svn_root%\SMV\Build\background
set in_smv=%svn_root%\SMV\Build\smokeview\intel_win_%platform%
set in_for_bundle=%svn_root%\SMV\for_bundle
set in_sh2bat=%svn_root%\SMV\Build\sh2bat\intel_win_64
set in_impi=%userprofile%\FIRE-LOCAL\LIBS\RUNTIME\WINDOWS_HYDRA2fix

set uploads=%svn_root%\Utilities\uploads
set basedir=%uploads%\%basename%

set out_bundle=%basedir%\firemodels
set out_bin=%out_bundle%\%fdsversion%\bin
set out_uninstall=%out_bundle%\%fdsversion%\Uninstall
set out_doc=%out_bundle%\%fdsversion%\Documentation
set out_guides="%out_doc%\Guides_and_Release_Notes"
set out_web="%out_doc%\FDS_on_the_Web"
set out_examples=%out_bundle%\%fdsversion%\Examples
set out_examples2=%svn_root%\Verification

set out_smv=%out_bundle%\%smvversion%
set out_textures=%out_smv%\textures

set fds_casessh=%svn_root%\Verification\FDS_Cases.sh
set fds_casesbat=%svn_root%\Verification\FDS_Cases.bat
set smv_casessh=%svn_root%\Verification\scripts\SMV_Cases.sh
set smv_casesbat=%svn_root%\Verification\scripts\SMV_Cases.bat
set wui_casessh=%svn_root%\Verification\scripts\WUI_Cases.sh
set wui_casesbat=%svn_root%\Verification\scripts\WUI_Cases.bat

set copyFDScases=%svn_root%\Utilities\Scripts\copyFDScases.bat
set copyCFASTcases=%svn_root%\Utilities\Scripts\copyCFASTcases.bat

set bundleinfo=%svn_root%\Utilities\Scripts\for_bundle

:: erase the temporary bundle directory if it already exists

if exist %basedir% rmdir /s /q %basedir%

mkdir %basedir%
mkdir %out_bundle%
mkdir %out_bundle%\%fdsversion%
mkdir %out_bundle%\%smvversion%
mkdir %out_bin%
mkdir %out_textures%
mkdir %out_doc%
mkdir %out_guides%
mkdir %out_web%
mkdir %out_examples%
mkdir %out_uninstall%

set release_version=%fdssmv_major_version%_win_%platform%
set release_version=

echo.
echo --- filling distribution directory ---
echo.


copy %in_for_bundle%\*.po                                                        %out_bin%\.>Nul

CALL :COPY  %fdsmpidir%\fds_mpi_win_%platform%.exe                               %out_bin%\fds.exe

CALL :COPY  %in_testmpi%\test_mpi.exe                                            %out_bin%\test_mpi.exe

CALL :COPY  %in_smv%\smokeview_win_%platform%.exe                                %out_smv%\smokeview.exe

CALL :COPY  %in_smokediff%\intel_win_%platform%\smokediff_win_%platform%.exe     %out_smv%\smokediff.exe

CALL :COPY  %in_smokezip%\intel_win_%platform%\smokezip_win_%platform%.exe       %out_smv%\smokezip.exe 

CALL :COPY  %in_dem2fds%\intel_win_%platform%\dem2fds_win_%platform%.exe         %out_smv%\dem2fds.exe 

CALL :COPY  %in_smvscriptdir%\jp2conv.bat                                        %out_smv%\jp2conv.bat 

CALL :COPY  %in_wind2fds%\intel_win_%platform%\wind2fds_win_%platform%.exe       %out_smv%\wind2fds.exe 

CALL :COPY  %in_fds2ascii%\intel_win_%platform%\fds2ascii_win_%platform%.exe     %out_bin%\fds2ascii.exe

CALL :COPY  %in_background%\intel_win_64\background.exe                          %out_bin%\background.exe

CALL :COPY %in_impi%\impi.dll          %out_bin%\impi.dll
CALL :COPY %in_impi%\mpiexec.hydra.exe %out_bin%\mpiexec.exe
CALL :COPY %in_impi%\pmi_proxy.exe     %out_bin%\pmi_proxy.exe
CALL :COPY %in_impi%\hydra_service.exe %out_bin%\hydra_service2.exe

CALL :COPY  %in_sh2bat%\sh2bat.exe %out_bin%\sh2bat.exe

echo.
echo --- copying auxillary files ---
echo.
CALL :COPY  %in_for_bundle%\objects.svo            %out_smv%\.
CALL :COPY  %in_for_bundle%\volrender.ssf          %out_smv%\.

CALL :COPY %in_intel_dll%\LIB64\libiomp5md.dll     %out_bin%\.
CALL :COPY  %in_for_bundle%\pthreadVC2_x64.dll     %out_smv%\.
CALL :COPY  %in_for_bundle%\glew32_x64.dll         %out_smv%\.

CALL :COPY  %in_for_bundle%\smokeview.ini          %out_smv%\.

echo copying textures
copy %in_for_bundle%\textures\*.jpg          %out_textures%\.>Nul
copy %in_for_bundle%\textures\*.png          %out_textures%\.>Nul

echo.
echo --- copying uninstaller ---
echo.
CALL :COPY  "%bundleinfo%\uninstall_fds.bat" "%out_uninstall%\uninstall_base.bat"
CALL :COPY  "%bundleinfo%\uninstall_fds2.bat" "%out_uninstall%\uninstall_base2.bat"
CALL :COPY  "%bundleinfo%\uninstall.bat"     "%out_uninstall%\uninstall.bat"
echo @echo off > "%out_uninstall%\uninstall.vbs"

CALL :COPY  "%in_setpath%\set_path64.exe"      "%out_uninstall%\set_path.exe"

echo.
echo --- copying FDS documentation ---
echo.

CALL :COPY  "%bundleinfo%\FDS_Release_Notes.htm" "%out_guides%\FDS_Release_Notes.htm"

CALL :COPY  %in_pdf%\FDS_Config_Management_Plan.pdf %out_guides%\.

CALL :COPY  %in_pdf%\FDS_User_Guide.pdf %out_guides%\.

CALL :COPY  %in_pdf%\FDS_Technical_Reference_Guide.pdf %out_guides%\.

CALL :COPY  %in_pdf%\FDS_Validation_Guide.pdf %out_guides%\.

CALL :COPY  %in_pdf%\FDS_Verification_Guide.pdf %out_guides%\.

echo.
echo --- copying Smokeview documentation ---
echo.

CALL :COPY %in_pdf%\SMV_User_Guide.pdf %out_guides%\.

CALL :COPY %in_pdf%\SMV_Technical_Reference_Guide.pdf %out_guides%\.

CALL :COPY %in_pdf%\SMV_Verification_Guide.pdf %out_guides%\.

echo.
echo --- copying starup shortcuts ---
echo.
 
CALL :COPY "%userprofile%\FDS-SMVwebpages\smv_readme.html" "%out_guides%\Smokeview_release_notes.html"

CALL :COPY "%bundleinfo%\Overview.html"                "%out_doc%\Overview.html"

CALL :COPY "%bundleinfo%\FDS_Web_Site.url"             "%out_web%\Official_Web_Site.url"

CALL :COPY "%bundleinfo%\Updates.url"                  "%out_web%\Software_Updates.url"

CALL :COPY "%bundleinfo%\Docs.url"                     "%out_web%\Documentation_Updates.url"

CALL :COPY "%bundleinfo%\FDS_Development_Web_Site.url" "%out_web%\Developer_Web_Site.url"

CALL :COPY "%bundleinfo%\discussion_group.url"         "%out_web%\Discussion_Group.url"

CALL :COPY "%bundleinfo%\issue_tracker.url"            "%out_web%\Issue_Tracker.url"

 
CALL :COPY %bundleinfo%\readme_examples.html           "%out_examples%\Examples notes.html"

echo.
echo --- copying example files ---
echo.

set outdir=%out_examples%
set QFDS=call %copyFDScases%
set RUNTFDS=call %copyFDScases%
set RUNCFAST=call %copyCFASTcases%

cd %out_examples2%
%in_sh2bat%\sh2bat %fds_casessh% %fds_casesbat%
call %fds_casesbat%>Nul

echo.
cd %out_examples2%
%in_sh2bat%\sh2bat %smv_casessh% %smv_casesbat%
call %smv_casesbat%>Nul
%in_sh2bat%\sh2bat %wui_casessh% %wui_casesbat%
call %wui_casesbat%>Nul

echo.
echo --- copying installer wrapup scripts ---
echo.

CALL :COPY  "%bundleinfo%\wrapup_fds_install.bat" "%out_bundle%\%fdsversion%\wrapup_fds_install.bat"

CALL :COPY  "%bundleinfo%\setup_fds_firewall.bat" "%out_bundle%\%fdsversion%\setup_fds_firewall.bat"

CALL :COPY  "%in_shortcut%\shortcut.exe"           "%out_bundle%\%fdsversion%\shortcut.exe"

CALL :COPY  "%in_setpath%\set_path64.exe"           "%out_bundle%\%fdsversion%\set_path.exe"

echo.
echo --- compressing distribution ---
echo.

cd %uploads%
if exist %basename%.zip erase %basename%.zip
cd %out_bundle%\%fdsversion%
wzzip -a -r -xExamples\*.csv -P ..\..\..\%basename%.zip * ..\%smvversion% > Nul

:: create an installation file from the zipped bundle directory

echo.
echo --- creating installer ---
echo.

cd %uploads%
echo Setup is about to install FDS %fds_version% and Smokeview %smv_version% > %bundleinfo%\message.txt
echo Press Setup to begin installation. > %bundleinfo%\main.txt
if exist %basename%.exe erase %basename%.exe
wzipse32 %basename%.zip -runasadmin -a %bundleinfo%\about.txt -st"FDS %fds_version% Smokeview %smv_version% Setup" -d "c:\Program Files\firemodels\FDS6" -c wrapup_fds_install.bat

echo.
echo --- installer built ---
echo.

cd %CURDIR%>Nul

GOTO :EOF

:COPY
set label=%~n1%~x1
set infile=%1
set infiletime=%~t1
set outfile=%2
IF EXIST %infile% (
   echo copying %label% %infiletime%
   copy %infile% %outfile% >Nul
) ELSE (
   echo.
   echo *** warning: %infile% does not exist
   echo.
   pause
)
exit /b
