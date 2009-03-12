"""
pyroplot_calcs.py

"""

import sys
import os
from math import sqrt

def mu_2sigma(x_data_set,y_data_set,diagnostic_level):
    epsilion_vals = []
    for i in range(len(x_data_set)):
        for j in range(len(x_data_set[i])):
            epsilion_vals = epsilion_vals+[(y_data_set[i][j]-x_data_set[i][j])/x_data_set[i][j]]
            
    mu_val = sum(epsilion_vals)/len(epsilion_vals)
    if diagnostic_level >= 3:
        print 'Mu:', mu_val
    
    sigma_val = sqrt(sum([(e_val-mu_val)**2 for e_val in epsilion_vals])/len(epsilion_vals))
    sigma_2_val = sigma_val*2
    if diagnostic_level >= 3:
        print '2 Sigma:', sigma_2_val
    
    return[mu_val,sigma_2_val]

def calc_min(d1_data,d2_data,d1_initial_value,d2_initial_value,diagnostic_level):
    if diagnostic_level >= 3:
        print "*** Compute Drop ***"
    
    temp_d1_data_values = [x for x in d1_data if x != -9999.0]
    d1_drop_value = float(d1_initial_value) - min(temp_d1_data_values)
    temp_d2_data_values = [x for x in d2_data if x != -9999.0]
    d2_drop_value = float(d2_initial_value) - min(temp_d2_data_values)
    
    if diagnostic_level >= 3:
        print "Data Set 1, Initial Value is:", d1_initial_value
        print "Data Set 1, Drop Value is:", d1_drop_value
        print "Data Set 2, Initial Value is:", d2_initial_value
        print "Data Set 2, Drop Value is:", d2_drop_value
        
    if diagnostic_level >= 3:
        print "\n*** Computing Drop Relative Difference ***"
    try:
        relative_difference = ((d2_drop_value-d1_drop_value)/d1_drop_value)
        if diagnostic_level >= 3:
            print "Min Relative Difference is:", relative_difference
    except:
        print "!!! Computation of Min relative_difference failed. !!!\nCheck source data for columns listed above."
        exit()
        
    return [d1_drop_value,d2_drop_value,relative_difference]

def calc_max(d1_data,d2_data,d1_initial_value,d2_initial_value,diagnostic_level):
    if diagnostic_level >= 3:
        print "*** Compute Rise ***"
        
    temp_d1_data_values = [x for x in d1_data if x != -9999.0]
    d1_rise_value = max(temp_d1_data_values) - float(d1_initial_value)
    temp_d2_data_values = [x for x in d2_data if x != -9999.0]
    d2_rise_value = max(temp_d2_data_values) - float(d2_initial_value)
    
    if diagnostic_level >= 3:
        print "Data Set 1, Initial Value is:", d1_initial_value
        print "Data Set 1, Rise Value is:", d1_rise_value
        print "Data Set 2, Initial Value is:", d2_initial_value
        print "Data Set 2, Rise Value is:", d2_rise_value
        
    if diagnostic_level >= 3:
        print "\n*** Computing Rise Relative Difference ***"
    try:
        relative_difference = ((d2_rise_value-d1_rise_value)/d1_rise_value)
        if diagnostic_level >= 3:
            print "Rise Relative Difference is:", relative_difference
    except:
        print "!!! Computation of Rise relative_difference failed. !!!\nCheck source data for columns listed above."
        exit()
    
    return [d1_rise_value,d2_rise_value,relative_difference]

def test_mu_2sigma():
    x_data_set = [[40.00,60.00,125.00,165.00,245.00],[300.00, 430.00, 500.00, 555.00, 600.00],[620.00, 615.00, 600.00, 653.00, 660.00, 646.00]]
    y_data_set = [[50.80,69.23,119.71,167.59,256.20],[356.02, 425.82, 504.18, 553.45, 602.47],[624.12, 624.60, 656.12, 670.82, 668.81, 660.29]]
    
    mu_2sigma(x_data_set,y_data_set,diagnostic_level)

def main():
    test_mu_2sigma()

if __name__ == '__main__':
    main()
