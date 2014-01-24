@echo off
set to=%1
set subject=%2
set message=%3

:: SMTP_USER_NAME, SMTP_USER_NAME_BASE, SMTP_USER_PASS are predefined environment variables

mailsend -to %to% -from %SMTP_USER_NAME% -ssl -smtp smtp.gmail.com -port 465 -sub %subject% -M %message% -q -auth-plain -user %SMTP_USER_NAME_BASE%

