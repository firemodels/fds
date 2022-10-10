@echo off
echo @echo off
for %%f in (*.smv) do (
  echo.
  for %%g in (%%~nf*.ini) do (
    echo smokeview -update_ini %%~ng.ini %%~nf.smv
  )
)
