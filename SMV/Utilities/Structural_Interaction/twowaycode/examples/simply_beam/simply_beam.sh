ansys150 -b nolist -j simply_beam_ansys -i read_geometry_simply_beam.ans
python simply_beam.py
ansys150 -b -j simply_beam_step1 -i simply_beam_step1.ans
ansys150 -b -j simply_beam_step2 -i simply_beam_step2.ans
../../intel_linux_64/twowaycode_linux_64 input_simply_beam.2way
ansys150 -b -j simply_beam_ansys_mechanical -i export_data_simply_beam.ans
ansys150 -b -j simply_beam_ansys_mechanical -i plot_simply_beam.ans
export DISPLAY=localhost:0.0 #
display150 -d POSTSCRIPT -j simply_beam_ansys #
epstopdf --outfile=simply_beam_temp.pdf pscr000.eps  # 
mv simply_beam_temp.pdf ../../scripts/SCRIPT_FIGURES/simply_beam_temp.pdf #
rm *.lock #
rm *.dat #
rm *_ansys* #
rm *.tmp #
rm *.eps #
rm *.log #
rm *.err #
rm *.db #
#rm *.dbb #