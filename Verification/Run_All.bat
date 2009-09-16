cd Atmospheric_Effects
for %%i in (*.fds) do %1 %%i
cd ..
cd Controls
for %%i in (*.fds) do %1 %%i
cd ..
cd Decaying_Isotropic_Turbulence
for %%i in (*.fds) do %1 %%i
cd ..
cd Detectors
for %%i in (*.fds) do %1 %%i
cd ..
cd Evacuation
for %%i in (*.fds) do %1 %%i
cd ..
cd Fires
for %%i in (*.fds) do %1 %%i
cd ..
cd Flowfields
for %%i in (*.fds) do %1 %%i
cd ..
cd Heat_Transfer
for %%i in (*.fds) do %1 %%i
cd ..
cd Miscellaneous
for %%i in (*.fds) do %1 %%i
cd ..
cd NS_Analytical_Solution
for %%i in (*.fds) do %1 %%i
cd ..
cd Pressure_Effects
for %%i in (*.fds) do %1 %%i
cd ..
cd Pyrolysis
for %%i in (*.fds) do %1 %%i
cd ..
cd Radiation
for %%i in (*.fds) do %1 %%i
cd ..
cd SMV_Script_Example
for %%i in (*.fds) do %1 %%i
cd ..
cd Species
for %%i in (*.fds) do %1 %%i
cd ..
cd Sprinklers_and_Sprays
for %%i in (*.fds) do %1 %%i
cd ..
cd Timing_Benchmarks
for %%i in (*.fds) do %1 %%i
cd ..
cd Visualization
%1 colorconv.fds
%1 plume5a.fds
%1 plume5b.fds
%1 plume5c.fds
%1 sillytexture.fds
%1 script_test.fds
%1 smoke_sensor.fds
%1 smoke_test.fds
%1 smoke_test2.fds
%1 thouse5.fds
cd ..
cd Wui
%1 fire_line.fds
