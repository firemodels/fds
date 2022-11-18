..\..\..\..\Build\impi_intel_win\fds_impi_intel_win simple_panel.fds 
cd ..
cd intel_win
echo Y | make_fds2ftmi.bat
cd ..
cd verification
..\intel_win\fds2ftmi_win simple_panel 0.05 2 0 600 2 0 0 1 simple_panel_verification 
python plots_verification.py

