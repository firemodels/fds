@echo off
for %%f in (*.ini) do (
  smokeview -update_ini %%f update
)