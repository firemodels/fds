"""
pyrograph_calcs.py

"""

import sys
import os
from math import log, sqrt, exp
import csv

tolerance = 1.0e-8

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
    

def delta_sigma(ind_data_set,dep_data_set,sigma_e,diagnostic_level):
    E_values = []
    M_values = []
    
    for i in range(len(ind_data_set)):
        for E_value in ind_data_set[i]:
            E_values.append(E_value)
        for M_value in dep_data_set[i]:
            M_values.append(M_value)
    
    tempfile = open('tempfile.csv', 'a')
    tempWriter = csv.writer(tempfile,delimiter=',')
    tempWriter.writerow(E_values)
    tempWriter.writerow(M_values)
    tempfile.close()
    
    #print 'E_Values:', E_values
    #print 'M_Values:', M_values
    
    # Find Length of data set.
    n = len(E_values)
    #print 'Length of set:', n
    
    if n == 1:
        print 'Error: Not enough scatter data to compute statistics.'
        exit()
    
    # Compute natural log of each value in data sets.
    ln_E_set = [log(abs(float(E_values[i]))) for i in range(len(E_values))]
    ln_M_set = [log(abs(float(M_values[i]))) for i in range(len(M_values))]
    #print 'ln_E_set:', ln_E_set
    #print 'ln_M_set:', ln_M_set
    
    # Compute 'M_hat' for each value in ln_E_set.
    M_hat = [sum(ln_M_set)/len(ln_M_set)-sum(ln_E_set)/len(ln_E_set)+ln_E_set[i] for i in range(len(ln_E_set))]
    #print 'M_hat:', M_hat
    
    # Compute 'u' values.
    u_temp = [((((ln_M_set[i])-(M_hat[i]))**2.0)/(n-1.0)) for i in range(len(ln_M_set))]
    u = sqrt(sum(u_temp))
    #print 'u:', u
    
    omega = sqrt(abs((u**2)-((sigma_e/2)**2)))
    #print 'omega:', omega
    
    # Compute 'delta'
    delta=exp(sum(ln_M_set)/len(ln_M_set)-sum(ln_E_set)/len(ln_E_set)+(omega**2/2))
    #print 'delta:', delta
    
    #Compute sigma
    sigma  = omega*delta
    #print 'sigma', sigma
    
    return[delta,sigma]

def calc_min(d1_data,d2_data,d1_initial_value,d2_initial_value,diagnostic_level):
    if diagnostic_level >= 3:
        print "*** Compute Drop ***"
    
    temp_d1_data_values = [x for x in d1_data if x != -9999.0]
    d1_drop_value = float(d1_initial_value) - min(temp_d1_data_values)
    temp_d2_data_values = [x for x in d2_data if x != -9999.0]
    d2_drop_value = float(d2_initial_value) - min(temp_d2_data_values)
    
    if d1_drop_value == 0.0:
        if diagnostic_level >= 3:
            print 'Data set 1 has a zero drop value, setting value to lower tolerance.'
        d1_drop_value = tolerance
    if d2_drop_value == 0.0:
        if diagnostic_level >= 3:
            print 'Data set 2 has a zero drop value, setting value to lower tolerance.'
        d2_drop_value = tolerance
    
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
        print "!!! Computation of Min relative_difference failed. !!!\n \
            Check source data for columns listed above."
        exit()
        
    return [d1_drop_value,d2_drop_value,relative_difference]

def calc_max(d1_data,d2_data,d1_initial_value,d2_initial_value,diagnostic_level):
    if diagnostic_level >= 3:
        print "*** Compute Rise ***"
        
    temp_d1_data_values = [x for x in d1_data if x != -9999.0]
    temp_d2_data_values = [x for x in d2_data if x != -9999.0]
    
    if len(temp_d1_data_values) == 0:
        print "Error in d1 scatter data set, length of array is 0."
        exit()
    elif len(temp_d2_data_values) == 0:
        print "Error in d2 scatter data set, length of array is 0."
        exit()
    else:
        d1_rise_value = max(temp_d1_data_values) - float(d1_initial_value)
        d2_rise_value = max(temp_d2_data_values) - float(d2_initial_value)
    
    if d1_rise_value == 0.0:
        if diagnostic_level >= 3:
            print 'Data set 1 has a zero rise value, setting to lower tolerance.'
        d1_rise_value = tolerance
    if d2_rise_value == 0.0:
        if diagnostic_level >= 3:
            print 'Data set 2 has a zero rise value, setting to lower tolerance.'
        d2_rise_value = tolerance
        
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
        print "!!! Computation of Rise relative_difference failed. !!!\n \
            Check source data for columns listed above."
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
