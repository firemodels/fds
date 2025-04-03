# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 13:21:19 2024

@author: Julio Cesar Silva - Rio on Fire
"""
import os
import time
import pandas as pd
import numpy as np
import sys
import subprocess

# Function to check if a file exists
def check_file(file_path):
    return os.path.exists(file_path)

# Function for linear interpolation
def interpolate_temperature(time_input,curve):
    temperature_interpolated = np.interp(time_input, curve['Time'], curve['Temperature (°C)'])
    return temperature_interpolated


def read_files(filepath,head):
    values = pd.read_csv(filepath,header=head)
    return values

def check_file_readability(file_path,check):
    try:
        with open(file_path, 'r') as file:
            check=1
        print(f"File '{file_path}' can be read successfully.")
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except PermissionError:
        print(f"Permission denied to read file '{file_path}'.")
    except Exception as e:
        print(f"An error occurred while trying to read file '{file_path}': {e}")
    return check

restart='false'
if len(sys.argv) <= 2:
    print('few arguments!')
if len(sys.argv) == 3:
    chid = sys.argv[1]
    inp = sys.argv[2]
    wdir = os.getcwd()
    print('Arguments read: '+ chid + ',' + inp + ','+ wdir)
elif len(sys.argv) == 4:
    chid = sys.argv[1]
    inp = sys.argv[2]
    if sys.argv[3] == 'restart':
        restart='true'
        wdir=os.getcwd()
        print('Arguments read: '+ chid + ',' + inp + ', restart '+ wdir)
    else:
        wdir=sys.argv[3]
        print('Arguments read: '+ chid + ',' + inp + ','+ wdir)
elif len(sys.argv) >= 5:
    chid = sys.argv[1]
    inp = sys.argv[2]
    if sys.argv[3] != 'restart':
        print('wrong arguments! Restarting or not?')
        print("Exiting the program...")
        sys.exit(0)
    else:
        restart='true'
        wdir=sys.argv[4]
        print('Arguments read: '+ chid + ',' + inp + ', restart '+ wdir)
elif len(sys.argv) >= 6:
    print('too many arguments!')
    print("Exiting the program...")
    sys.exit(0)

if check_file(wdir+'/fds_wake.txt') == True: os.remove(wdir+'/fds_wake.txt')
file_path1=wdir+'/'+inp+'.csv'
input_values= read_files(file_path1, None)
file_path=wdir+"/"+chid+"_devc.csv"
   
# Run FDS in the background
subprocess.Popen(['mpiexec', '-n', '1', '../../Build/impi_intel_linux_openmp/fds_impi_intel_linux_openmp', chid+'.fds'])

start=0
dt_devc=2
t_end=600
temp_ini=20
nexttime=start+dt_devc

devc=pd.DataFrame({"Firesize": ['FTC_1','FTC_2','FTC_3','FTC_4'],
                    })

curves=read_files('ext_heartbeat_std_curve.csv', 0)

if restart == 'true':
    stored_values=read_files(wdir+"/stored_values_"+chid+".csv", 0)
    last_time=stored_values['Time'][len(stored_values['Time'])-1]
    nexttime=last_time+dt_devc
    time_interpolated = np.arange(0, t_end+dt_devc, dt_devc)  # Generate time values at dt_devc-second intervals
    curves_interpolated = pd.DataFrame({'Time': time_interpolated})
    i=0
    while i < len(devc.columns):
        # Interpolate temperature values at 30-second intervals
        temperature_interpolated = np.interp(time_interpolated, curves['time'], curves[str(devc.columns[i])])
        # Append to DataFrame with interpolated values
        curves_interpolated = pd.concat([curves_interpolated, pd.DataFrame(temperature_interpolated)], axis=1)
        curves_interpolated=curves_interpolated.rename(columns={0:str(devc.columns[i])})
        i=i+1
else:
    stored_values=pd.DataFrame(columns=['Time'])
    stored_values.at[0,'Time']=0
    stored_values['std_Curve']=temp_ini

    time_interpolated = np.arange(0, t_end+dt_devc, dt_devc)  # Generate time values at dt_devc-second intervals
    #Create a new DataFrame with interpolated values
    curves_interpolated = pd.DataFrame({'Time': time_interpolated})
    i=0
    while i < len(devc.columns):
        stored_values[str(devc.columns[i])]=curves[str(devc.columns[i])][0]
        stored_values[str(devc.columns[i])+'_FTC']=temp_ini
        stored_values[str(devc.columns[i])]=input_values[2][i]
        # Interpolate temperature values at 30-second intervals
        temperature_interpolated = np.interp(time_interpolated, curves['time'], curves[str(devc.columns[i])])
        # Append to DataFrame with interpolated values
        curves_interpolated = pd.concat([curves_interpolated, pd.DataFrame(temperature_interpolated)], axis=1)
        curves_interpolated=curves_interpolated.rename(columns={0:str(devc.columns[i])})
        i=i+1
    

# Main loop
if restart=='true':
    i=int(nexttime/dt_devc)
else:
    i=1
while nexttime <= t_end+dt_devc:
    if check_file(file_path):
        # Perform math operations here if the file exists
        print("File found! Reading...")
        check=0
        check=check_file_readability(file_path,check)
        if check == 1:
            result = read_files(file_path, 1)
            if result['Time'][len(result)-1] < nexttime:
                print("Next time step still not achieved. Waiting for 2 seconds...")
                time.sleep(2)
            elif result.iloc[[i]]['Time'][i] >= nexttime:
                j=0 
                while j < len(devc.columns)  :
                    k=0
                    n=0
                    quantity=0

                    while k < len(devc):
                        if devc.iloc[k,j] != 0:
                            quantity=quantity+result.iloc[i][devc.iloc[k,j]]
                            n=n+1
                        k=k+1
                    quantity=quantity/n
                    quantity_target= np.interp(nexttime+dt_devc,curves_interpolated['Time'],curves_interpolated[str(devc.columns[j])])                    
                    # Perform linear interpolation to get the temperature at the input time
                    interpolated_temperature = np.interp(nexttime,curves_interpolated['Time'],curves_interpolated[str(devc.columns[j])])
                    if quantity_target < quantity:
                        input_values.at[j,2]=input_values.at[j,2] * 0.95
                    elif quantity_target > quantity :
                        input_values.at[j,2]=input_values.at[j,2] * 1.05
                    elif quantity_target == quantity :
                        input_values.at[j,2]=input_values.at[j,2]
                    stored_values.at[i,str(devc.columns[j])+'_FTC']=quantity
                    stored_values.at[i,str(devc.columns[j])]=input_values.at[j,2]
                    j=j+1
                j=0
                while j < len(devc.columns)  :
                    print("Updated "+input_values.at[j,1]+" values to "+str(input_values.at[j,2])+" !")
                    j=j+1
                input_values.to_csv(file_path1, index=False,header=None)
                f = open(wdir+"/fds_wake.txt", "w")
                f.close()
                stored_values.at[i,'Time']=nexttime
                stored_values.at[i,'std_Curve']=interpolated_temperature
                if check_file(wdir+"/stored_values_"+chid+".csv"):os.remove(wdir+"/stored_values_"+chid+".csv")
                stored_values.to_csv(wdir+"/stored_values_"+chid+".csv",index=None)
                i=i+1
                nexttime=nexttime+dt_devc
                if nexttime >= t_end:
                    print("Exiting the program...")
                    sys.exit(0)
                                
        else:
            time.sleep(3)
    else:
        print("File not found. Waiting for 2 seconds...")
        time.sleep(2)  # Pause for 2 second
        