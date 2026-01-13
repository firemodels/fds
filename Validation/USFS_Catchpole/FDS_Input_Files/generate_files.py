# Python script to generate input files for Catchpole cases
# Uses template.txt input file and replaces keywords based on parameters in Test_Matrix.csv

import pandas as pd

# Read the template text file
with open('template.txt', 'r') as template_file:
   template = template_file.read()

# Iterate through table of experimental cases and replace keywords
metadata = pd.read_csv('Test_Matrix.csv')
for i, row in metadata.iterrows():
   textout = template

   textout = textout.replace('JOBID',row['Test'])
   t_end = round(min(999.,12/row['R']),0)
   textout = textout.replace('endtime',str(t_end))
   textout = textout.replace('windspeed',str(row['U']))
   textout = textout.replace('vegheight',str(row['delta']))
   textout = textout.replace('moistfrac',str(row['M']))
   textout = textout.replace('svratio',str(row['s']))

   if (row['Test'][:2]=='MF'):
      fuel_type = 'Pine Sticks'
      density = 442.
      heat_of_reaction = 659.
      ignition_end = 40.
      ignition_p1 = ignition_end + 1
   elif (row['Test'][:4]=='EXSC'):
      fuel_type = 'Coarse Excelsior'
      density = 398.
      heat_of_reaction = 711.
      ignition_end = 10.
      ignition_p1 = ignition_end + 1
   elif (row['Test'][:4]=='PPMC'):
      fuel_type = 'Pine Needles'
      density = 510.
      heat_of_reaction = 609.
      ignition_end = 20.
      ignition_p1 = ignition_end + 1
   elif (row['Test'][:2]=='EX'):
      fuel_type = 'Regular Excelsior'
      density = 398.
      heat_of_reaction = 711.
      ignition_end = 10.
      ignition_p1 = ignition_end + 1
   
   mpv=round(row['beta']*density,2)
   textout = textout.replace('masspervolume',str(mpv))
   textout = textout.replace('vegdensity',str(density))
   textout = textout.replace('heatofreaction',str(heat_of_reaction))
   textout = textout.replace('ignitionend',str(ignition_end))
   textout = textout.replace('ignitionp1',str(ignition_p1))
   textout = textout.replace('fueltype',fuel_type)

   with open(row['Test']+'.fds', 'w') as output_file:
      output_file.write(textout)

