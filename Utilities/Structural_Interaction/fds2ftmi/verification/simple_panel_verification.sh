../../../../Build/impi_intel_linux_64/fds_impi_intel_linux_64 simple_panel.fds 
cd ..
cd intel_linux_64
echo Y | make_fds2ftmi.bat
cd ..
cd verification
../intel_linux_64/fds2ftmi_linux_64 simple_panel 0.05 2 0 598 2 0 0 1 simple_panel_verification 
python plots_verification.py