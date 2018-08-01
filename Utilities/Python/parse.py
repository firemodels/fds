import numpy as np
# McDermott
# 6-08-2009
# parse.m

# [S] = parse(character_string)

# Uses textscan to parse a character string that is delimited with a "|".


def parse(character_string=None):

    S=np.asarray([], dtype='object')
# temp/parse.m:11
    cell_array=np.asarray([character_string.split('|')],dtype='object').T
# temp/parse.m:12
    S=cell_array[:].T
# temp/parse.m:13
    return S