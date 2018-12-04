echo off
set FDSMAJORVERSION=6
set FDSEDITION=FDS6
set SMVEDITION=SMV6

set fdsversion=%FDSEDITION%
set smvversion=%SMVEDITION%

if "%env_defined%" == "1" goto endif_env_defined
set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist2
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use smv/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist2

call %envfile%
:endif_env_defined

set      in_impi=%userprofile%\fire-notes\INSTALL\LIBS\RUNTIME\MPI_%INTELVERSION%
set in_intel_dll=%userprofile%\fire-notes\INSTALL\LIBS\WINDOWS\%INTELVERSION%

set SVNROOT=%svn_root%
set fdsdir=%svn_root%\fds\Build\intel_win_64
set fdsmpidir=%svn_root%\fds\Build\impi_intel_win_64
set fdsmpidirdb=%svn_root%\fds\Build\impi_intel_win_64_db
set basename=%fds_version%-%smv_version%_win64
set hashfile=%svn_root%\smv\Build\hashfile\intel_win_64\hashfile_win_64.exe
if NOT exist %hashfile% (
  echo ***warning: %hashfile% does not exist
  echo Bundle will not contain hashes of application files
  pause
)

set in_pdf=%userprofile%\.bundle\pubs
set in_shortcut=%userprofile%\fire-notes\INSTALL\repoexes
set in_for_bundle=%svn_root%\smv\Build\Bundle\for_bundle


set uploads=%svn_root%\fds\Build\Bundle\uploads
set basedir=%uploads%\%basename%

set out_bundle=%basedir%\firemodels
set out_bin=%out_bundle%\%fdsversion%\bin
set out_fdshash=%out_bundle%\%fdsversion%\bin\hash
set out_uninstall=%out_bundle%\%fdsversion%\Uninstall
set out_doc=%out_bundle%\%fdsversion%\Documentation
set out_guides="%out_doc%\Guides_and_Release_Notes"
set out_web="%out_doc%\FDS_on_the_Web"
set out_examples=%out_bundle%\%fdsversion%\Examples
set fds_examples=%svn_root%\fds\Verification
set smv_examples=%svn_root%\smv\Verification

set out_smv=%out_bundle%\%smvversion%
set out_textures=%out_smv%\textures
set out_smvhash=%out_smv%\hash

set fds_casessh=%svn_root%\fds\Verification\FDS_Cases.sh
set fds_casesbat=%svn_root%\fds\Verification\FDS_Cases.bat
set smv_casessh=%svn_root%\smv\Verification\scripts\SMV_Cases.sh
set smv_casesbat=%svn_root%\smv\Verification\scripts\SMV_Cases.bat
set wui_casessh=%svn_root%\smv\Verification\scripts\WUI_Cases.sh
set wui_casesbat=%svn_root%\smv\Verification\scripts\WUI_Cases.bat

set copyFDScases=%svn_root%\fds\Utilities\Scripts\copyFDScases.bat
set copyCFASTcases=%svn_root%\fds\Utilities\Scripts\copyCFASTcases.bat

set bundleinfo=%svn_root%\fds\Build\Bundle\for_bundle

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

mkdir %out_fdshash%
mkdir %out_smvhash%

set release_version=%FDSMAJORVERSION%_win_64
set release_version=

echo.
echo --- filling distribution directory ---
echo.


copy %in_for_bundle%\*.po                                                     %out_bin%\.>Nul

if "%fds_debug%" == "1" (
  CALL :COPY  %fdsmpidirdb%\fds_impi_win_64_db.exe                            %out_bin%\fds_db.exe
)
CALL :COPY  %fdsmpidir%\fds_impi_win_64.exe                                   %out_bin%\fds.exe
CALL :COPY  %svn_root%\fds\Utilities\fds2ascii\intel_win_64\fds2ascii_win_64.exe %out_bin%\fds2ascii.exe
CALL :COPY  %svn_root%\smv\Build\background\intel_win_64\background.exe       %out_bin%\background.exe
CALL :COPY  %svn_root%\fds\Utilities\test_mpi\impi_intel_win\test_mpi.exe     %out_bin%\test_mpi.exe

CALL :COPY  %svn_root%\smv\Build\smokeview\intel_win_64\smokeview_win_64.exe  %out_smv%\smokeview.exe
CALL :COPY  %svn_root%\smv\Build\smokediff\intel_win_64\smokediff_win_64.exe  %out_smv%\smokediff.exe
CALL :COPY  %svn_root%\smv\Build\smokezip\intel_win_64\smokezip_win_64.exe    %out_smv%\smokezip.exe 
CALL :COPY  %svn_root%\smv\Build\dem2fds\intel_win_64\dem2fds_win_64.exe      %out_smv%\dem2fds.exe 
CALL :COPY  %svn_root%\smv\Build\hashfile\intel_win_64\hashfile_win_64.exe    %out_smv%\hashfile.exe 
CALL :COPY  %svn_root%\smv\Build\wind2fds\intel_win_64\wind2fds_win_64.exe    %out_smv%\wind2fds.exe 
CALL :COPY  %svn_root%\smv\Build\hashfile\intel_win_64\hashfile_win_64.exe    %out_smv%\hashfile.exe 
CALL :COPY  %svn_root%\smv\scripts\jp2conv.bat                                %out_smv%\jp2conv.bat

set curdir=%CD%
cd %out_bin%
%hashfile% fds.exe        >  hash\fds_%fds_version%.exe.sha1
%hashfile% fds2ascii.exe  >  hash\fds2ascii_%fds_version%.exe.sha1
%hashfile% background.exe >  hash\background_%fds_version%.exe.sha1
%hashfile% test_mpi.exe   >  hash\test_mpi_%fds_version%.exe.sha1
cd hash
cat *.sha1              >  %uploads%\%basename%.sha1

cd %out_smv%
%hashfile% hashfile.exe   >  hash\hashfile_%smv_version%.exe.sha1
%hashfile% smokeview.exe  >  hash\smokeview_%smv_version%.exe.sha1
%hashfile% smokediff.exe  >  hash\smokediff_%smv_version%.exe.sha1
%hashfile% smokezip.exe   >  hash\smokezip_%smv_version%.exe.sha1
%hashfile% dem2fds.exe    >  hash\dem2fds_%smv_version%.exe.sha1
%hashfile% wind2fds.exe   >  hash\wind2fds_%smv_version%.exe.sha1
cd hash
cat *.sha1              >>  %uploads%\%basename%.sha1

cd %curdir%
CALL :COPY %in_intel_dll%\libiomp5md.dll     %out_bin%\libiomp5md.dll
CALL :COPY %in_intel_dll%\libfabric.dll      %out_bin%\libfabric.dll

CALL :COPY %in_impi%\impi.dll                %out_bin%\impi.dll
CALL :COPY %in_impi%\mpiexec.hydra.exe       %out_bin%\mpiexec.exe
CALL :COPY %in_impi%\pmi_proxy.exe           %out_bin%\pmi_proxy.exe 
CALL :COPY %in_impi%\hydra_service.exe       %out_bin%\hydra_service2.exe

CALL :COPY  %svn_root%\smv\Build\sh2bat\intel_win_64\sh2bat.exe %out_bin%\sh2bat.exe

echo.
echo --- copying auxillary files ---
echo.
CALL :COPY  %in_for_bundle%\objects.svo            %out_smv%\.
CALL :COPY  %in_for_bundle%\volrender.ssf          %out_smv%\.

CALL :COPY  %in_for_bundle%\smokeview.ini          %out_smv%\.

echo copying textures
copy %in_for_bundle%\textures\*.jpg          %out_textures%\.>Nul
copy %in_for_bundle%\textures\*.png          %out_textures%\.>Nul

echo.
echo --- copying uninstaller ---
echo.
CALL :COPY  "%bundleinfo%\uninstall_fds.bat"  "%out_uninstall%\uninstall_base.bat"
CALL :COPY  "%bundleinfo%\uninstall_fds2.bat" "%out_uninstall%\uninstall_base2.bat"
CALL :COPY  "%bundleinfo%\uninstall.bat"      "%out_uninstall%\uninstall.bat"
echo @echo off > "%out_uninstall%\uninstall.vbs"

CALL :COPY  "%svn_root%\smv\Build\set_path\intel_win_64\set_path64.exe"      "%out_uninstall%\set_path.exe"

echo.
echo --- copying FDS documentation ---
echo.

CALL :COPY  "%svn_root%\webpages\FDS_Release_Notes.htm"  %out_guides%\FDS_Release_Notes.htm
CALL :COPY  %in_pdf%\FDS_Config_Management_Plan.pdf      %out_guides%\.
CALL :COPY  %in_pdf%\FDS_User_Guide.pdf                  %out_guides%\.
CALL :COPY  %in_pdf%\FDS_Technical_Reference_Guide.pdf   %out_guides%\.
CALL :COPY  %in_pdf%\FDS_Validation_Guide.pdf            %out_guides%\.
CALL :COPY  %in_pdf%\FDS_Verification_Guide.pdf          %out_guides%\.

echo.
echo --- copying Smokeview documentation ---
echo.

CALL :COPY %in_pdf%\SMV_User_Guide.pdf                %out_guides%\.
CALL :COPY %in_pdf%\SMV_Technical_Reference_Guide.pdf %out_guides%\.
CALL :COPY %in_pdf%\SMV_Verification_Guide.pdf        %out_guides%\.

echo.
echo --- copying startup shortcuts ---
echo.
 
CALL :COPY "%svn_root%\webpages\smv_readme.html"       "%out_guides%\Smokeview_release_notes.html"
CALL :COPY "%bundleinfo%\Overview.html"                "%out_doc%\Overview.html"
CALL :COPY "%bundleinfo%\FDS_Web_Site.url"             "%out_web%\Official_Web_Site.url"

echo.
echo --- copying example files ---

set outdir=%out_examples%
set QFDS=call %copyFDScases%
set RUNTFDS=call %copyFDScases%
set RUNCFAST=call %copyCFASTcases%

cd %fds_examples%
%svn_root%\smv\Build\sh2bat\intel_win_64\sh2bat %fds_casessh% %fds_casesbat%
call %fds_casesbat%>Nul

cd %smv_examples%
%svn_root%\smv\Build\sh2bat\intel_win_64\sh2bat %smv_casessh% %smv_casesbat%
call %smv_casesbat%>Nul
%svn_root%\smv\Build\sh2bat\intel_win_64\sh2bat %wui_casessh% %wui_casesbat%
call %wui_casesbat%>Nul

echo.
echo --- copying installer wrapup scripts ---
echo.

CALL :COPY  "%bundleinfo%\wrapup_fds_install.bat" "%out_bundle%\%fdsversion%\wrapup_fds_install.bat"

CALL :COPY  "%bundleinfo%\setup_fds_firewall.bat" "%out_bundle%\%fdsversion%\setup_fds_firewall.bat"

CALL :COPY  "%in_shortcut%\shortcut.exe"           "%out_bundle%\%fdsversion%\shortcut.exe"

CALL :COPY  "%svn_root%\smv\Build\set_path\intel_win_64\set_path64.exe"           "%out_bundle%\%fdsversion%\set_path.exe"

echo.
echo --- compressing distribution ---

cd %uploads%
if exist %basename%.zip erase %basename%.zip
cd %out_bundle%\%fdsversion%
wzzip -a -r -xExamples\*.csv -P ..\..\..\%basename%.zip * ..\%smvversion% > Nul

:: create an installation file from the zipped bundle directory

echo.
echo --- creating installer ---

cd %uploads%
echo Setup is about to install FDS %fds_version% and Smokeview %smv_version% > %bundleinfo%\message.txt
echo Press Setup to begin installation. > %bundleinfo%\main.txt
if exist %basename%.exe erase %basename%.exe
wzipse32 %basename%.zip -runasadmin -a %bundleinfo%\about.txt -st"FDS %fds_version% Smokeview %smv_version% Setup" -d "c:\Program Files\firemodels\FDS6" -c wrapup_fds_install.bat
%hashfile% %basename%.exe   >>  %uploads%\%basename%.sha1

echo.
echo --- installer built ---

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
