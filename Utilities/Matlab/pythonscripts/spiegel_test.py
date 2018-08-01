import math as m
import numpy as np
def spiegel_test(x=None):

    # compute pvalue under null of x normally distributed;
# x should be a vector;
    xm=np.mean(x)
# temp/spiegel_test.m:4
    xs=np.std(x, ddof = 1)
# temp/spiegel_test.m:5
    xz=(x - xm) / xs
# temp/spiegel_test.m:6
    xz2=xz ** 2
# temp/spiegel_test.m:7
    N=np.sum(np.multiply(xz2,np.log(xz2)))
# temp/spiegel_test.m:8
    n=x.size
# temp/spiegel_test.m:9
    ts=(N - np.dot(0.73,n)) / (np.dot(0.8969,m.sqrt(n)))
# temp/spiegel_test.m:10
    
    pval=1 - abs(m.erf(ts / m.sqrt(2)))
# temp/spiegel_test.m:11
    
    return pval