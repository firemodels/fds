@echo off
set curdir=%CD%
set fromdir=%user%@blaze.nist.gov:/var/www/html/firebot/manuals
set todir="%userprofile%\Google Drive\FDS-SMV Newest Manuals"

cd %todir%
pscp %fromdir%/FDS_Configuration_Management_Plan.pdf .
pscp %fromdir%/FDS_Technical_Reference_Guide.pdf .
pscp %fromdir%/FDS_User_Guide.pdf .
pscp %fromdir%/FDS_Validation_Guide.pdf .
pscp %fromdir%/FDS_Verification_Guide.pdf .
pscp %fromdir%/SMV_Technical_Reference_Guide.pdf .
pscp %fromdir%/SMV_User_Guide.pdf .
pscp %fromdir%/SMV_Verification_Guide.pdf .

cd %curdir%