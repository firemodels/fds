## Generating the input files for USFS/Catchpole

There are 354 USFS/Catchpole experiments listed in the test matrix in the Appendix of

W.R. Catchpole et al., "Rate of Spread of Free-Burning Fires in Woody Fuels in a Wind Tunnel," Combustion Science and Technology, Vol. 131, pp. 1-37, 1998.

These parameters are listed in the file called Test_Matrix.csv. A short Python program called generate_files.py will take the generic base input file in template.txt and generate 354 individual input files using replacements for the key parameters.

