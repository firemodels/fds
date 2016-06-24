:: http://stackoverflow.com/a/15885133/1683264

@echo off
setlocal enabledelayedexpansion
set prog=%1
set progdesc=%2

set progdir=C:\Program Files\firemodels\%prog%
if NOT exist "%progdir%" mkdir "%progdir%"

cd %progdir%
set chooser=%temp%\folder_chooser.exe
if exist !chooser! del !chooser!
>"%temp%\c.cs" echo using System;using System.Windows.Forms;
>>"%temp%\c.cs" echo class dummy{[STAThread]
>>"%temp%\c.cs" echo public static void Main^(^){
>>"%temp%\c.cs" echo FolderBrowserDialog f=new FolderBrowserDialog^(^);
>>"%temp%\c.cs" echo f.RootFolder=Environment.SpecialFolder.ProgramFiles;
>>"%temp%\c.cs" echo f.SelectedPath=System.Environment.CurrentDirectory;
>>"%temp%\c.cs" echo f.Description="Choose location to install %progdesc%.";
>>"%temp%\c.cs" echo f.ShowNewFolderButton=true;
>>"%temp%\c.cs" echo if^(f.ShowDialog^(^)==DialogResult.OK^){Console.Write^(f.SelectedPath^);}}}
for /f "delims=" %%I in ('dir /b /s "%windir%\microsoft.net\*csc.exe"') do (
    if not exist "!chooser!" "%%I" /nologo /out:"!chooser!" "%temp%\c.cs" 2>NUL
)
call %chooser%

del "%temp%\c.cs"
del "%temp%\folder_chooser.exe"
