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
import re

import sys
sys.stdout = sys.stderr  # send all print() to stderr

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
        #print(f"File '{file_path}' can be read successfully.")
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except PermissionError:
        print(f"Permission denied to read file '{file_path}'.")
    except Exception as e:
        print(f"An error occurred while trying to read file '{file_path}': {e}")
    return check


# Function to check if DT_DEVC and DT_FLUSH exist and their value
def update_dump_params(fds_file_path, dt_devc):
    modified = False
    with open(fds_file_path, 'r') as file:
        lines = file.readlines()

    # Regex pattern to capture DUMP parameters
    dump_pattern = re.compile(
        r'^(\s*&DUMP\b)(.*?)(/.*)$', 
        re.IGNORECASE
    )
    param_pattern = re.compile(
        r'\b(DT_FLUSH|DT_DEVC)\s*=\s*([\d.]+)\b', 
        re.IGNORECASE
    )

    for i in range(len(lines)):
        line = lines[i]
        if line.strip().upper().startswith('&DUMP'):
            match = dump_pattern.match(line)
            if not match:
                continue

            prefix = match.group(1)
            params = match.group(2).strip()
            suffix = match.group(3)
            params_dict = {}
            needs_update = False

            # Parse existing parameters
            for name, value in param_pattern.findall(params):
                params_dict[name.upper()] = value

            # Check/update DT_FLUSH
            if 'DT_FLUSH' in params_dict:
                if float(params_dict['DT_FLUSH']) != dt_devc:
                    needs_update = True
            else:
                needs_update = True

            # Check/update DT_DEVC
            if 'DT_DEVC' in params_dict:
                if float(params_dict['DT_DEVC']) != dt_devc:
                    needs_update = True
            else:
                needs_update = True

            if needs_update:
                # Remove existing parameters
                clean_params = re.sub(
                    r'\b(DT_FLUSH|DT_DEVC)\s*=\s*[\d.]+\b,?\s*', 
                    '', 
                    params, 
                    flags=re.IGNORECASE
                ).strip(', ')

                # Add new parameters
                new_params = f" DT_DEVC={dt_devc}, DT_FLUSH={dt_devc}"
                if clean_params:
                    new_params = f"{clean_params}, {new_params}"

                # Rebuild the line with original formatting
                lines[i] = f"{prefix}{new_params} {suffix}\n"
                modified = True
                print(f"Line {i+1}: Updated DT_DEVC and DT_FLUSH to {dt_devc}")

    if modified:
        with open(fds_file_path, 'w') as file:
            file.writelines(lines)
        print(f"\nFile '{fds_file_path}' was updated to consider DT_FLUSH <= DT_DEVC in DUMP line")
    else:
        print("\nNo changes needed - DT_DEVC and DT_FLUSH already match prescribed values")


restart='false'
chid = 'ext_heartbeat_std_curve' # sys.argv[1]
inp = 'ext_heartbeat_mass_flux' #sys.argv[2]
wdir = '.' #os.getcwd()
#print('Arguments read: '+ chid + ',' + inp + ','+ wdir)

if check_file(wdir+'/fds_wake.txt') == True: os.remove(wdir+'/fds_wake.txt')
file_path1=wdir+'/'+inp+'.csv'
file_path=wdir+"/"+chid+"_devc.csv"
file_pathsmv=wdir+"/"+chid+".smv"


# if devc and smv files already exists from other run, they will be deleted to avoid read results from that
if check_file(file_path) == True: os.remove(file_path)
if check_file(file_pathsmv) == True: os.remove(file_pathsmv)
# if input files already exists from other run, it will be deleted to avoid read results from that. Also, a file will be created with the correct initial value
if check_file(file_path1) == True: os.remove(file_path1)
pd.DataFrame(columns=('RAMP','FIRESIZE', 0.05)).to_csv(file_path1, index=False)

input_values= read_files(file_path1, None)

start=0
dt_devc=0.2
t_end=60
temp_ini=20
nexttime=start+dt_devc

update_dump_params(wdir+"/"+chid+".fds", dt_devc)

print('Alright, now you can run FDS with '+chid+'.fds')
   
# Run FDS in the background
#subprocess.Popen(['mpiexec', '-n', '1', '../../Build/impi_intel_linux_openmp/fds_impi_intel_linux_openmp', chid+'.fds'])
#subprocess.Popen(['$QFDS', '-d', 'Miscellaneous', chid+'.fds'])


devc=pd.DataFrame({"Firesize": ['FTC_1','FTC_2','FTC_3','FTC_4'],})
curves=read_files('ext_heartbeat_std_curve.csv', 0)

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
    
# a loop to check with FDS started to run
i=0
while check_file(file_pathsmv) == False:
    if i==0: print("FDS not started yet, smv file not detected")
    time.sleep(5)
    i=i+1
    #wait for FDS start for how many seconds? i*5 
    if i > 1440: # 120 minutes
        print("FDS did not started in "+int(i*5/60)+" minutes")
        print("Exiting the program...")
        sys.exit(0)

print("FDS started detected - Main loop started")

# Main loop
i=1
counter=0
while nexttime <= t_end+dt_devc:
    check=0
    check=check_file_readability(file_path,check)
    if check == 1:
        result = read_files(file_path, 1)
        if result['Time'][len(result)-1] < nexttime:
            if counter == 0: print("Next time step still not achieved. Waiting...")
            time.sleep(1)
            counter=counter+1
            if counter == 7200: # 7200 seconds or 2 hours
                print("Waited "+str(counter/60)+" minutes for the Next time step and nothing happened, maybe FDS crashed")
                print("Exiting the program...")
                sys.exit(0)
        elif result.iloc[[i]]['Time'][i] >= nexttime:
            counter=0
            print(str(result.iloc[[i]]['Time'][i])+" time step read in devc file")
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
                print("Updated "+input_values.at[j,1]+" values to "+str(input_values.at[j,2]))
                j=j+1
            input_values.to_csv(file_path1, index=False,header=None)
            f = open(wdir+"/fds_wake.txt", "w")
            f.close()
            stored_values.at[i,'Time']=nexttime
            stored_values.at[i,'std_Curve']=interpolated_temperature
            if check_file(wdir+"/"+chid+"_stored_values.csv"):os.remove(wdir+"/"+chid+"_stored_values.csv")
            stored_values.to_csv(wdir+"/"+chid+"_stored_values.csv",index=None)
            i=i+1
            nexttime=nexttime+dt_devc
            if nexttime >= t_end:
                print("Exiting the program...")
                sys.exit(0)
                            
    else:
        time.sleep(2)
    
