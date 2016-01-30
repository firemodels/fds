@echo off
set testimage=xplume5c_test.png
set curdir=%CD%
if exist %testimage% erase %testimage%
cd ..\..
set reporoot=%CD%
set visdir=%reporoot%\Verification\Visualization
set SMV=smokeview
cd %visdir%
%SMV% -script plume5c2.ssf plume5c
cd %curdir%
if exist %testimage% echo test image creation succeeded
if not exist %%testimage% echo test image creation failed

