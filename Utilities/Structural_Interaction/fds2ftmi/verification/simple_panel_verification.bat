..\..\..\..\Build\impi_intel_win_64\fds_impi_win_64 simple_panel.fds 
cd ..
cd intel_win_64
echo Y | make_fds2ftmi.bat
cd ..
cd verification
..\intel_win_64\fds2ftmi_win_64 simple_panel 0.05 2 0 600 2 0 0 1 simple_panel_verification 
python plots_verification.py

