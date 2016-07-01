from __future__ import division

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from pylab import *
from matplotlib import rcParams
from scipy import integrate

filenames = ['simply_input.fds', '2fds_input.geo']
with open('simply_beam.fds', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
	
with open('simply_beam.fds', 'a') as outfile:
	outfile.write("MATL_ID='BEAM', SURF_ID='BEAM' / \n")
	outfile.write("&TAIL / \n")