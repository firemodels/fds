import numpy as np
import smop_util
# McDermott
# 6-11-2009
# stripcell.m

# Stip cell array of empty cells


def stripcell(CELL_ARRAY=None):

    C=np.asarray([], dtype='object')
# temp/stripcell.m:9
    i=1
# temp/stripcell.m:10
    for j in range(1,CELL_ARRAY.size+1):
        if isinstance(CELL_ARRAY[j-1], str):
            C = smop_util.safe_set(C,(i-1,),CELL_ARRAY[j-1])
# temp/stripcell.m:13
            i=i + 1
# temp/stripcell.m:14
    
    return C