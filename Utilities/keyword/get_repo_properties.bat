@echo off

:: obtain various properties of the directory repo_dir
:: These properties are the revision and revision date 
:: of the directory and the current date.

set repo_dir=%1

set revision=unknown

set revision_date=unknown
set revision_time=unknown

set build_date=unknown
set build_time=unknown

set havesvn=1
set havegit=1

set validsvn=0
set validgit=0

set temp1="%temp%\temp.txt"
set temp1c="%temp%\tempc.txt"

set CURDIR=%CD%

:: check to see if %repo_dir% exists

if NOT exist %repo_dir% (
  echo *** warning: The directory %repo_dir% does not exist.
  exit /b 1
)

:: ----------------- make sure various required software tools are available --------------------------

:: looking for svn

svn 1> %temp1% 2>&1
type %temp1% | find /i /c "not recognized" > %temp1c%
set flag=
set /p flag=<%temp1c%
if %flag% == 1 (
  set havesvn=0
)

cd %repo_dir%

:: is this a valid svn repository

if %havesvn% == 0 goto skiphavesvn
  set validsvn=1
  svn info 1> %temp1% 2>&1
  type %temp1% | find /i /c "not a working copy" > %temp1c%
  set glag=
  set /p glag=<%temp1c%
  if %glag% == 1 (
    set svn_revision=invalid
    cd %CURDIR%
    set validsvn=0
    if %havegit% == 0 (
      echo "*** warning: %repo_dir% is not a valid svn repository"
      exit /b 1
    )
  )
:skiphavesvn

:: looking for git

git 1> %temp1% 2>&1
type %temp1% | find /i /c "not recognized" > %temp1c%
set flag=
set /p flag=<%temp1c%
if %flag% == 1 (
  set havegit=0
  if %havesvn% == 0 (
    echo *** warning: both svn and git were not found
    exit /b 1
  )
)

:: is this a valid git repository

if %validsvn% == 0 (
  if %havegit% == 1 (
    set validgit=1
    git log . 1> %temp1% 2>&1
    type %temp1% | find /i /c "Not a git repository" > %temp1c%
    set flag=
    set /p flag=<%temp1c%
    if %flag% == 1 (
      set svn_revision=invalid
      cd %CURDIR%
      set validgit=0
      echo *** warning: %repo_dir% is not a valid git repository
      exit /b 1
    ) 
  )
)

:: looking for head (used only if this is notonly used with git)

if %validgit% == 1 (
  head -h 1> %temp1% 2>&1
  type %temp1% | find /i /c "not recognized" > %temp1c%
  set flag=
  set /p flag=<%temp1c%
  if %flag% == 1 (
    echo *** warning: head was not found.
    exit /b 1
  ) 
)

:: looking for tail (only used with git)

if %validgit% == 1 (
  tail -h 1> %temp1% 2>&1
  type %temp1% | find /i /c "not recognized" > %temp1c%
  set flag=
  set /p flag=<%temp1c%
  if %flag% == 1 (
    echo *** warning: tail was not found.
    exit /b 1
  )
)

:: looking for gawk

gawk 1> %temp1% 2>&1
type %temp1% | find /i /c "not recognized" > %temp1c%
set flag=
set /p flag=<%temp1c%
if %flag% == 1 (
  echo *** warning: gawk was not found.
  exit /b 1
)

:: ----------------- get properties --------------------------

:: get revision number

if %validsvn% ==1 (
  svn info 2>&1 | find /i "Last Changed Rev:" | gawk -F" " "{print $4}" > %temp1%
  set /p revision=<%temp1%
)
if %validgit% ==1 (
  git log --abbrev-commit . 2>&1 | head -1 | gawk -F" " "{print $2}" > %temp1%
  set /p revision=<%temp1%
)

:: get date and time of latest repository commit

if %validsvn% ==1 (
  svn info 2>&1 | find /i "Last Changed Date:" | gawk -F" " "{print $4}" > %temp1%
  set /p revision_date=<%temp1%
  svn info 2>&1 | find /i "Last Changed Date:" | gawk -F" " "{print $5}" |gawk -F":" "{print $1\":\"$2}"  > %temp1%
  set /p revision_time=<%temp1%
)
if %validgit% ==1 (
  git log --date=short . 2>&1 | head -3 | tail -1 | gawk -F" " "{print $2}" > %temp1%
  set /p revision_date=<%temp1%
  git log . 2>&1 | head -3 | tail -1 | gawk -F" " "{print $5}" |gawk -F":" "{print $1\":\"$2}"  > %temp1%
  set /p revision_time=<%temp1%
)

:: get current date time

echo %date% 2>&1 | gawk -F" " "{print $2}" | gawk -F"/" "{print $3\"-\"$1\"-\"$2}" > %temp1% 
set /p build_date=<%temp1%
echo %time% 2>&1 | gawk -F":" "{print $1\":\"$2}" > %temp1%
set /p build_time=<%temp1%

cd %CURDIR%
