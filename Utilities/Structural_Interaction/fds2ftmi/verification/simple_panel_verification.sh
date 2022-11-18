../../../../Build/impi_intel_linux/fds_impi_intel_linux simple_panel.fds 
cd ..
cd intel_linux
echo Y | make_fds2ftmi.bat
cd ..
cd verification
../intel_linux/fds2ftmi_linux simple_panel 0.05 2 0 598 2 0 0 1 simple_panel_verification 
python plots_verification.py