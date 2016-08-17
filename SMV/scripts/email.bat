@echo off
set to=%1
set subject=%2
set file=%3

:: SMTP_xxx variables are predefined environment variables

set SSL=
if %SMTP_PORT% == 465 (
  set SSL=-ssl
)

mailsend -to %to% -from %SMTP_USER_NAME% -smtp %SMTP_SERVER% %SSL% -port %SMTP_PORT% -sub %subject% -attach %file%,text/plain,i -q -auth-plain -user %SMTP_USER_NAME_BASE%
