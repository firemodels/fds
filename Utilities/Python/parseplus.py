import numpy as np
# McDermott
# 12-30-2015
# parseplus.m

# [S] = parseplus(character_string)

# Uses textscan to parse a character string that is delimited with a "+".


def parseplus(character_string=None):

    S=np.asarray([], dtype='object')
# temp/parseplus.m:11
    cell_array=np.asarray([str(character_string).split('+')],dtype='object').T
# temp/parseplus.m:12
    S=cell_array[:].T
# temp/parseplus.m:13
    return S