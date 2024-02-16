''' 
J. Hodges
05-13-2023
generate_scaling_pyrolysis_cases.py

This script will automatically generate FDS input files using the
scaling pyrolysis model, import the predictions, and compare them
with experimental results.

Detailed instructions on the input file can be found at:
<enter wiki url>

Output:
    FDS predictions in: out/Scaling_Pyrolysis/
    Validation guide table .tex in: <enter here>
    Validation dataplot input in: <enter here>
    Figures in: fds/Validation/Scaling_Pyrolysis/figures/

Input:
    specificationFile - string specifying the location of the input file
    figoutdir - string documenting the local path to save figures
    FDSDIR - location of FDS for running the models. Note this must be
        set as an environmental variable. If not set, the script will
        try to find FDS within the system PATH.

Special switch_id tags:
    Note copied from dataplot.m, currently only d and s are implemented.

    'd' -- Process this data line as usual (exception: see 'o' below)
    's' -- Skip this line
    'o' -- Add 'o' in the switch_id column (first column) of FDS_validation_dataplot_inputs.csv to process "only" these lines.
    'f' -- Follow the previous line and "hold on" the figure window, adding this line to the current plot.
    'g' -- Generate plot, but ignore in scatplot.  Good for cases under development.

'''
import numpy as np
import os, sys, argparse, shutil, glob
import subprocess
import pandas as pd
from collections import defaultdict
from matplotlib.lines import Line2D
import scipy

def getColors():
    ''' This function defines colors for plotting.
    '''
    colors = np.array([[0, 65, 101],
                       [229, 114, 0],
                       [136, 139, 141],
                       [170, 39, 44],
                       [119, 197, 213],
                       [161, 216, 132],
                       [255, 200, 69],
                       [101, 0, 65],
                       [0, 229, 114],
                       [141, 136, 139],
                       [44, 170, 39],
                       [213, 119, 197],
                       [132, 161, 216],
                       [69, 255, 200],
                       [65, 101, 0],
                       [114, 0, 229]], dtype=np.float32)
    colors = colors/255
    return colors

def runModel(outdir, outfile, mpiProcesses, fdsdir, fdscmd, printLiveOutput=False):
    ''' This function will run fds with an input file
    '''
    my_env = os.environ.copy()
    my_env['I_MPI_ROOT'] = fdsdir+"\\mpi"
    my_env['PATH'] = fdsdir + ';' + my_env['I_MPI_ROOT'] + ';' + my_env["PATH"]
    my_env['OMP_NUM_THREADS'] = '1'
    
    if printLiveOutput is True:
        process = subprocess.Popen([fdscmd, outfile, ">", "log.err"], cwd=r'%s'%(outdir), env=my_env, shell=False)
    else:
        process = subprocess.Popen([fdscmd, outfile, ">", "log.err"], cwd=r'%s'%(outdir), env=my_env, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
    
    out, err = process.communicate()
    errcode = process.returncode   
    return out, err, errcode

def plotResults_exp(data, exp_data, fluxes, validationTimeColumns, validationHrrpuaColumns, lw=3, fs=16):
    ''' This function will generate a plot comparing experimental
    results with the scaling predictions.
    '''
    fig = plt.figure(figsize=(12,8))
    time = data['Time']
    colors = getColors()
    tMax = 5
    for i, flux in enumerate(fluxes):
        scaling = data['"HRRPUA-%02d"'%(flux)].values
        if type(exp_data) is dict:
            plt.plot(exp_data[validationTimeColumns[i]]/60, exp_data[validationHrrpuaColumns[i]], '-', linewidth=lw, label='exp', color=colors[i])
            ind = np.where(exp_data[validationHrrpuaColumns[i]] > 0)[-1]
            tMax = max([tMax, np.nanmax(exp_data[validationTimeColumns[i]][ind])/60])
        else:
            plt.plot(exp_data[validationTimeColumns[i]].values/60, exp_data[validationHrrpuaColumns[i]].values, '-', linewidth=lw, label='exp', color=colors[i])
            ind = np.where(exp_data[validationHrrpuaColumns[i]].values > 0)[-1]
            tMax = max([tMax, np.nanmax(exp_data[validationTimeColumns[i]].values[ind])/60])
        plt.plot(time.values/60, scaling, '--', linewidth=lw, label='Scaling-%02d'%(flux), color=colors[i])
    tMax = np.ceil(tMax/5)*5
    plt.xlabel('Time (min)', fontsize=fs)
    plt.ylabel('HRRPUA (kW/m2)', fontsize=fs)
    plt.xlim(0, tMax)
    plt.tick_params(labelsize=fs)
    plt.legend(fontsize=fs)
    plt.grid()
    plt.tight_layout()
    return fig

def extractAnalysisData(filterWidth, data, exp_data, fluxes, validationTimeColumns, validationHrrpuaColumns):
    time = data['Time'].values
    model_dt = np.median(time[1:]-time[:-1])
    
    filterWidth_mod = int(filterWidth/model_dt)
    fil_mod = np.ones(filterWidth_mod)/filterWidth_mod
    
    uncertainty_dict = dict()
    
    for i, flux in enumerate(fluxes):
        scaling = data['"HRRPUA-%02d"'%(flux)].values
        scaling_filtered = np.convolve(scaling, fil_mod, mode='same')
        if type(exp_data) is dict:
            expTime = exp_data[validationTimeColumns[i]]
            expHRR = exp_data[validationHrrpuaColumns[i]]
        else:
            expTime = exp_data[validationTimeColumns[i]].values
            expHRR = exp_data[validationHrrpuaColumns[i]].values
        exp_dt = np.min(expTime[1:]-expTime[:-1])
        exp_times2 = np.linspace(0, expTime.max(), int(np.round(expTime.max()/exp_dt)+1))
        exp_hrr2 = np.interp(exp_times2, expTime, expHRR)
        
        filterWidth_exp = int(filterWidth/exp_dt)
        fil_exp = np.ones(filterWidth_exp)/filterWidth_exp
        exp_filtered = np.convolve(exp_hrr2, fil_exp, mode='same')
        
        uncertainty_dict['%03d-EXP'%(flux)] = np.nanmax(exp_filtered)
        uncertainty_dict['%03d-MOD'%(flux)] = np.nanmax(scaling_filtered)
    return uncertainty_dict
        
def findHeaderLength(lines):
    ''' This is a helper function to dynamically find the
    length of a header in csv data
    '''
    counter = 0
    headerCheck = True
    while headerCheck and counter < 100:
        line = (lines[counter].decode('utf-8')).replace('\r\n','')
        while line[-1] == ',': line = line[:-1]
        try:
            [float(y) for y in line.split(',')]
            counter = counter - 1
            headerCheck = False
        except:
            counter = counter + 1
    if counter < 100:
        return counter
    else:
        print("Unable to find header length, returning 0")
        return 0

def cleanDataLines(lines2, headerLines):
    ''' This is a helper function to clean data rows
    '''
    lines = lines2[headerLines+1:]
    for i in range(0, len(lines)):
        line = (lines[i].decode('utf-8')).replace('\r\n','')
        while line[-1] == ',': line = line[:-1]
        lines[i] = [float(y) for y in line.split(',')]
    return lines

def load_csv(modeldir, chid, suffix='_devc', labelRow=-1):
    ''' This function imports a csv output by FDS
    '''
    file = "%s%s%s%s.csv"%(modeldir, os.sep, chid, suffix)
    f = open(file, 'rb')
    lines = f.readlines()
    f.close()
    headerLines = findHeaderLength(lines)
    if labelRow == -1:
        header = (lines[headerLines].decode('utf-8')).replace('\r\n','').replace('\n','').split(',')
    else:
        header = (lines[labelRow].decode('utf-8')).replace('\r\n','').replace('\n','').split(',')
    dataLines = cleanDataLines(lines, headerLines)
    data = pd.DataFrame(dataLines, columns=header,)
    return data


def findLimits(times, HRRs, energyCutoff1=0.001, energyCutoff2=1.01):
    ''' This function extracts the burning duration data from
    a cone calorimeter dataset. This is based on two cutoffs for the
    total energy released. The energy cutoff for the time to ignition
    is the same as used in findIgnitionTime and is arbitrary. The
    energy cutoff used on the trailing end is dynamically calculated
    based on the data curve. A seed value can be set as a start cutoff.
    By default, finds the last time where HRPPUA > 0.
    '''
    v = np.cumsum(HRRs)
    ind1 = 0 
    counter = 0
    while ind1 == 0:
        try:
            ind1 = np.where(v < np.nanmax(v)*energyCutoff1)[0][-1]
        except:
            energyCutoff1 = energyCutoff1*2
        counter += 1
        if counter > 20:
            ind1 = 0
            break
    ind2 = v.shape[0]
    counter = 0
    while ind2 == v.shape[0]:
        try:
            ind2 = np.where(v > np.nanmax(v)*energyCutoff2)[0][0]
        except:
            energyCutoff2 = energyCutoff2*0.99
        counter += 1
        if counter > 20:
            ind2 = v.shape[0]
            break
    times_trimmed = times[ind1:ind2]
    hrrs_trimmed = HRRs[ind1:ind2]
    tign = times[ind1]
    return tign, times_trimmed, hrrs_trimmed

def findIgnitionTime(times, HRRs, energyCutoff1=0.001):
    ''' This function finds the time to ignition in a cone calorimeter
    dataset. Ignition is assumed to occur when a threshold of the
    total energy has been released. The default threshold of 0.001
    is arbitrary but provides reasonable agreement with the tests
    analyzed here.
    '''
    v = np.cumsum(HRRs)
    ind1 = 0
    counter = 0
    while ind1 == 0 and counter < 100:
        try:
            ind1 = np.where(v < np.nanmax(v)*energyCutoff1)[0][-1]
        except:
            energyCutoff1 = energyCutoff1*2
        counter += 1
    tign = times[ind1]
    return tign

def interpolateExperimentalData(times, HRRs, targetDt=False, filterWidth=False,
                                hrrThresh=1., numPoints=10000):
    ''' This function interpolates the experimental data to a
    fixed time interval.
    '''
    dt = np.nanmedian(times[1:]-times[:-1])
    if filterWidth is not False:
        filterWidth = int(filterWidth/dt)
        fil = np.ones(filterWidth)/filterWidth
        HRRs = np.convolve(HRRs, fil, mode='same')
    
    if targetDt is not False:
        dt = targetDt
    else:
        dt = np.nanmedian(times[1:]-times[:-1])
    tmax = np.round(times.max()/dt)*dt
    tmin = np.round(times.min()/dt)*dt
    targetTimes = np.linspace(tmin, tmax, int((tmax-tmin)/dt + 1))
    targetHRRs = np.interp(targetTimes, times, HRRs)
    
    t_peak = times[np.argmax(HRRs)]
    q_peak = np.nanmax(HRRs)
    
    
    targetTimes = np.append(targetTimes, t_peak)
    targetHRRs = np.append(targetHRRs, q_peak)
    
    targetHRRs = targetHRRs[np.argsort(targetTimes)]
    targetTimes = np.sort(targetTimes)
    
    '''
    totalEnergy = np.trapz(HRRs, times)
    totalEnergy2 = np.trapz(targetHRRs, targetTimes)
    
    
    while abs(1-(totalEnergy/totalEnergy2)) > 0.01:
        targetDt = targetDt*0.9
        targetTimes = np.linspace(tmin, tmax, int((tmax-tmin)/targetDt + 1))
        targetHRRs = np.interp(targetTimes, times, HRRs)
        targetTimes = np.append(targetTimes, t_peak)
        targetHRRs = np.append(targetHRRs, q_peak)
        
        targetHRRs = targetHRRs[np.argsort(targetTimes)]
        targetTimes = np.sort(targetTimes)
        totalEnergy2 = np.trapz(targetHRRs, targetTimes)
    '''
    
    reconstructedHRR = np.interp(times, targetTimes, targetHRRs)
    while (np.max(abs(reconstructedHRR - HRRs)) > hrrThresh) and (targetTimes.shape[0] < numPoints):
        ind = np.argmax(abs(reconstructedHRR - HRRs))
        t_diff = times[ind]
        q_diff = HRRs[ind]
        
        targetTimes = np.append(targetTimes, t_diff)
        targetHRRs = np.append(targetHRRs, q_diff)
        
        targetHRRs = targetHRRs[np.argsort(targetTimes)]
        targetTimes = np.sort(targetTimes)
        reconstructedHRR = np.interp(times, targetTimes, targetHRRs)
    #print('%d, %0.2f'%(targetTimes.shape[0], np.max(abs(reconstructedHRR - HRRs))))
    return targetTimes, targetHRRs

def getRepresentativeHrrpua(HRRPUA, time, factor=0.5):
    ''' This function calculates a representative HRRPUA for use
    in estimating the flame heat flux. HRR contains a time series
    of HRRPUA data from a cone calorimeter experiment. All data
    points above a threshold percentile (default 50%) are used
    to calculate the average HRRPUA. The threshold is intended
    to remove leading and trailing parts of the curve which
    would artificially lower the HRRPUA. The threshold of 50%
    is arbitrary but provides reasonable agreement on the cases
    evaluated here.
    '''
    dts = time[1:]-time[:-1]
    dt = np.min(dts[np.nonzero(dts)])
    tmin = time.min()
    tmax = time.max()
    time_i = np.linspace(tmin, tmax, int((tmax-tmin)/dt + 1))
    hrrpua_i = np.interp(time_i, time, HRRPUA)
    representativeHRRPUA = hrrpua_i[hrrpua_i > HRRPUA.max()*factor].mean()
    return representativeHRRPUA

def estimateExposureFlux(coneExposure, representativeHRRPUA):
    ''' Estimates the exposure flux at a cone exposure by using
    an empirical estimate of the flame heat flux. The empirical
    estimate of the flame heat flux was developed by running
    steady-state FDS simulations at fixed HRRPUAs in a cone
    calorimeter configuration. The simulations used a fixed
    surface heat transfer coefficient of 15 W/m2-K and a fixed
    gas phase radiative fraction of 0.35.
    '''
    exposureFlux = 55*(representativeHRRPUA**0.065) - 50 + coneExposure
    return exposureFlux

def buildFdsFile(chid, coneExposure, e, k, rho, cp, Tign, d, time, hrrpua, tend,
                 HFs, front_h, hfs_tign=False,
                 ignitionMode='Temperature', outputTemperature=False,
                 calculateDevcDt=True, devc_dt=1.):
    ''' Generate a solid phase only FDS input file representing cone
    experimens at different exposures given a reference curve and
    material properties. The configuration can support up to 9
    thermal exposures, configured in a 3x3 grid.
    
    Notes:
        1. Zero heat transfer coefficient at the surface. This is
           intentional because the calculated flame heat flux which
           is included in the IGNITION-RAMPs includes convection.
        2. The estimated flame heat flux currently assumes a surface
           heat transfer coefficient of 15 W/m2-K and a gas phase
           radiative fraction of 0.35.
        3. If the ignition temperature is unknown, it can be calculated
           by comparing with experimental times to ignition. Changing
           the input of Tign to 'Calculated' will tell FDS to save
           out the WALL TEMPERATURE data needed to extract this
           information.
        4. All samples are assumed to have 0.5in / 12.7 mm of ceramic
           fiber insulation behind them.
    '''
    hrrpua_ref = getRepresentativeHrrpua(hrrpua, time)
    qref = estimateExposureFlux(coneExposure, hrrpua_ref)
    
    tempOutput = '.TRUE.' if outputTemperature else '.FALSE.'
    DT_DEVC = devc_dt
    if calculateDevcDt:
        NFRAMES = 1200/1.
        DT_DEVC = tend/NFRAMES
    if ignitionMode == 'Time': Tign = 20
    txt = "&HEAD CHID='%s', /\n"%(chid)
    txt = txt+"&TIME DT=1., T_END=%0.1f /\n"%(tend)
    txt = txt+"&DUMP DT_CTRL=%0.1f, DT_DEVC=%0.1f, DT_HRR=%0.1f, SIG_FIGS=4, SIG_FIGS_EXP=2, /\n"%(DT_DEVC, DT_DEVC, DT_DEVC)
    txt = txt+"&MISC SOLID_PHASE_ONLY=.TRUE., TMPA=27., /\n"
    txt = txt+"&MESH ID='MESH', IJK=3,3,3, XB=0.,0.3,0.,0.3,0.,0.3, /\n"
    txt = txt+"&REAC ID='PROPANE', FUEL='PROPANE', /\n"
    txt = txt+"&MATL ID='BACKING', CONDUCTIVITY=0.10, DENSITY=65., EMISSIVITY=0.9, SPECIFIC_HEAT=1.14, /\n"
    #txt = txt+"&MATL ID='BACKING', CONDUCTIVITY=0.2, DENSITY=585., EMISSIVITY=1., SPECIFIC_HEAT=0.8, /\n"
    txt = txt+"&MATL ID='SAMPLE', CONDUCTIVITY=%0.4f, DENSITY=%0.1f, EMISSIVITY=%0.4f, SPECIFIC_HEAT=%0.4f, /\n"%(k, rho, e, cp)
    
    prevTime=-1e6
    for i in range(0, len(time)):
        if (time[i]-prevTime) < 0.0001:
            #txt = txt+"&RAMP ID='CONE-RAMP', T=%0.4f, F=%0.1f, /\n"%(time[i]-time[0]+0.0001, hrrpua[i])
            pass
        else:
            txt = txt+"&RAMP ID='CONE-RAMP', T=%0.4f, F=%0.1f, /\n"%(time[i]-time[0], hrrpua[i])
        prevTime = time[i]
    y = -0.05
    for i, hf in enumerate(HFs):
        hf_ign = estimateHrrpua(coneExposure, hrrpua_ref, hf)
        if i%3 == 0: y = y + 0.1
        XYZ = [((i % 3))*0.1+0.05, y, 0.0]
        XB = [XYZ[0]-0.05, XYZ[0]+0.05, XYZ[1]-0.05, XYZ[1]+0.05, 0.0,0.0]
        
        txt = txt+"&SURF ID='SAMPLE-%02d', EXTERNAL_FLUX=1., "%(hf)
        txt = txt+"HEAT_TRANSFER_COEFFICIENT=%0.4f, HEAT_TRANSFER_COEFFICIENT_BACK=10., "%(front_h)
        txt = txt+"HRRPUA=1., IGNITION_TEMPERATURE=%0.1f, MATL_ID(1:2,1)='SAMPLE','BACKING', "%(Tign)
        txt = txt+"RAMP_EF='IGNITION_RAMP-%02d', RAMP_Q='CONE-RAMP', "%(hf)
        txt = txt+"REFERENCE_HEAT_FLUX=%0.4f, REFERENCE_HEAT_FLUX_TIME_INTERVAL=1.,"%(qref)
        txt = txt+'THICKNESS(1:2)=%0.8f,%0.8f, /\n'%(d, 0.0254/2)
        
        if ignitionMode == 'Temperature':
            txt = txt+"&RAMP ID='IGNITION_RAMP-%02d', T=%0.1f, F=%0.4f, DEVC_ID='IGNITION_DEVC-%02d', /\n"%(hf, 0.0, hf, hf)
            txt = txt+"&RAMP ID='IGNITION_RAMP-%02d', T=%0.1f, F=%0.4f, /\n"%(hf, 1.0, hf_ign)
        else:
            txt = txt+"&RAMP ID='IGNITION_RAMP-%02d', T=%0.1f, F=%0.4f, /\n"%(hf, 0.0, hf_ign)
            txt = txt+"&RAMP ID='IGNITION_RAMP-%02d', T=%0.1f, F=%0.4f, /\n"%(hf, 1.0, hf_ign)
        
        txt = txt+"&OBST ID='SAMPLE-%02d', SURF_ID='SAMPLE-%02d', XB="%(hf, hf)
        for x in XB:
            txt = txt+"%0.4f,"%(x)
        if ignitionMode == 'Time':
            txt = txt+"DEVC_ID='TIGN-%02d'"%(hf)
        txt = txt+', /\n'
        
        txt = txt+"&DEVC ID='WALL TEMPERATURE-%02d', INITIAL_STATE=.FALSE., IOR=3, OUTPUT=%s, "%(hf, tempOutput)
        txt = txt+"QUANTITY='WALL TEMPERATURE', SETPOINT=%0.1f, XYZ=%0.4f,%0.4f,%0.4f, /\n"%(Tign, XYZ[0], XYZ[1], XYZ[2])
        
        txt = txt+"&CTRL ID='IGNITION-CTRL-%02d', FUNCTION_TYPE='ANY', INPUT_ID='WALL TEMPERATURE-%02d', /\n"%(hf, hf)
        if ignitionMode == 'Time':
            txt = txt+"&DEVC ID='TIGN-%02d', XYZ=0,0,0, SETPOINT=%0.4f, QUANTITY='TIME', INITIAL_STATE=.FALSE., /\n"%(hf, hfs_tign[hf])
            
        txt = txt+"&DEVC ID='IGNITION_DEVC-%02d', CTRL_ID='IGNITION-CTRL-%02d', IOR=3, OUTPUT=.FALSE., QUANTITY='CONTROL', "%(hf,hf)
        txt = txt+"XYZ=%0.4f,%0.4f,%0.4f, /\n"%(XYZ[0], XYZ[1], XYZ[2])
        
        txt = txt+"&DEVC ID='HRRPUA-%02d', IOR=3, QUANTITY='HRRPUA', SPEC_ID='PROPANE', "%(hf)
        txt = txt+"XYZ=%0.4f,%0.4f,%0.4f, /\n"%(XYZ[0], XYZ[1], XYZ[2])
        
        txt = txt+"&DEVC ID='IGNITION-TIME-%02d', NO_UPDATE_DEVC_ID='IGNITION_DEVC-%02d', OUTPUT=.FALSE., "%(hf,hf)
        txt = txt+"QUANTITY='TIME', XYZ=%0.4f,%0.4f,%0.4f, /\n"%(XYZ[0], XYZ[1], XYZ[2])
        
                        
    return txt



def estimateHrrpua(cone_ref, hrrpua_ref, cone_flux):
    ''' Calculate the scaled heat flux based on the reference hrrpua,
    reference cone exposure, and scaled cone exposure. An iterative
    process is used since the flame heat flux depends on the scaled
    hrrpua.
    '''
    
    exposureFlux = estimateExposureFlux(cone_ref, hrrpua_ref)
    scaledFlux = exposureFlux - cone_ref + cone_flux
    
    prevFlux = scaledFlux
    scaled_hrrpua = (scaledFlux/exposureFlux)*hrrpua_ref
    scaledFlux = estimateExposureFlux(cone_flux, scaled_hrrpua)
    
    while abs(prevFlux - scaledFlux) > 0.01:
        prevFlux = scaledFlux
        scaled_hrrpua = (scaledFlux/exposureFlux)*hrrpua_ref
        scaledFlux = estimateExposureFlux(cone_flux, scaled_hrrpua)
    return scaledFlux


def getSeriesLatexReference(series):
    if series == 'Aalto_Woods':
        ref = '\label{Aalto_Woods_Description}'
    elif series == 'JH_FRA':
        ref ='\label{Aalto_Woods_Description}'
    elif series == 'JH_NIJ':
        ref = '\label{Aalto_Woods_Description}'
    elif series == "FAA_Polymers":
        ref ='\label{Aalto_Woods_Description}'
    elif series == "FPL_Materials":
        ref ='\label{Aalto_Woods_Description}'
    elif series == "FSIR_Materials":
        ref ='\label{Aalto_Woods_Description}'
    elif series == 'RISE_Materials':
        ref ='\label{Aalto_Woods_Description}'

def findFds():
    ''' First check if FDSDIR environmental variable is defined. If not, print
    warning then use which to look for a checklist. If not found anywhere
    error out.
    '''
    fdsDir = os.getenv('FDSDIR')
    if fdsDir is not None: return fdsDir, 'fds'
    print("Warning FDSDIR environmental variable not set. Trying to find FDS in path.")
    checklist = ['fds', 'fds_ompi_gnu_linux']
    for check in checklist:
        fdsPath = shutil.which(check)
        if fdsPath is not None:
            fdsDir = os.sep.join(fdsPath.split(os.sep)[:-1]) + os.sep
            print("FDS found in %s"%(fdsDir))
            return fdsDir, check
    print("Warning FDS not found")
  
def calculateUncertainty(x, y):
    sigma_e = 0.075
    mask = np.logical_and(~np.isnan(np.array(x, dtype=float)),
                          ~np.isnan(np.array(y, dtype=float)))
    if np.var(np.log(y[mask] / x[mask])) - (sigma_e**2) < 0:
        sigma_m = sigma_e
    else:
        sigma_m = (np.var(np.log(y[mask] / x[mask])) - (sigma_e**2))**0.5
    #sigma_m2 = np.var(np.log(y / x)) / 2
    sigma_m = np.nanmax([sigma_m, sigma_e])
    delta = np.exp(np.mean(np.log(y[mask] / x[mask])) + (sigma_m**2)/2 - (sigma_e**2)/2)
    return delta, sigma_m, sigma_e, np.log(y/x)

def calculateUncertaintyBounds(flatx, flaty, flatFlux, split=False):
    d = pd.DataFrame(np.array([flatx, flaty, flatFlux]).T, columns=['exp','mod','flux'])
    d[d == 0] = np.nan
    #d2[d2 < 0] = np.nan
    mask = np.logical_and(~np.isnan(np.array(d.values[:,0], dtype=float)),
                          ~np.isnan(np.array(d.values[:,1], dtype=float)))
    d2 = d.loc[mask]
    if split:
        uniqueFluxes = np.unique(flatFlux)
        delta = dict()
        sigma_m = dict()
        num_points = dict()
        points = dict()
        for flux in uniqueFluxes:
            x = np.array(d.loc[d['flux'] == flux, 'exp'].values, dtype=float)
            y = np.array(d.loc[d['flux'] == flux, 'mod'].values, dtype=float)
            delta[flux], sigma_m[flux], sigma_e, points[flux] = calculateUncertainty(x, y)
            num_points[flux] = x.shape[0]
    else:
        (x, y) = (np.array(d['exp'].values, dtype=float), np.array(d['mod'].values, dtype=float))
        delta, sigma_m, sigma_e, points = calculateUncertainty(x, y)
        num_points = d2.shape[0]
    return delta, sigma_m, sigma_e, num_points, points

def getNewColors():
    colors2 = np.array([[216, 27, 96],
                        [30, 136, 229],
                        [0, 77, 64],
                        [255, 193, 7],
                        [216, 27, 216],
                        [27, 216, 27],
                        ]) / 255
    return colors2

def plotMaterialExtraction(x, y, f, label, diff=None, axmin=None, axmax=None, loglog=False, labelName=None):
    
    axmin2 = min([np.min(x), np.min(y)])
    axmax2 = min([np.max(x), np.max(y)])
    
    delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x, y, f, split=False)
    
    #print(label, delta, sigma_m, sigma_e)
    
    if axmin is not None:
        axmin2 = axmin
    if axmax is not None:
        axmax2 = axmax
        
    fig, ax = plt.subplots(figsize=(12, 10))
    if loglog:
        ax.set_yscale('log')
        ax.set_xscale('log')
    fs=24
    
    xcoords = np.array([axmin2, axmax2])
    ycoords = np.array([axmin2, axmax2])
    dashes=(10, 10)
    ax.plot(xcoords, ycoords, 'k', linewidth=2)
    ax.plot(xcoords, ycoords*(1+2*sigma_e), '--k', linewidth=2, dashes=dashes)
    ax.plot(xcoords, ycoords/(1+2*sigma_e), '--k', linewidth=2, dashes=dashes)
    
    ax.plot(xcoords, ycoords*delta, 'r', linewidth=2)
    ax.plot(xcoords, ycoords*delta*(1+2*sigma_m), '--r', linewidth=2, dashes=dashes)
    ax.plot(xcoords, ycoords*delta/(1+2*sigma_m), '--r', linewidth=2, dashes=dashes)
    
    markers = ['o', 's', 'd', '>', '<', '^']
    cmap = plt.cm.viridis
    
    colors2 = getNewColors()
    
    cinds = {25:0, 50:1, 75:2}
    
    
    mew = 3
    if diff is not None:
        cases = np.array(list(set(diff)))
        cases.sort()
        for j in range(0, len(f)):
            caseInd = np.where(cases == diff[j])[0][0]
            #c = 0 if diff[j] > 0 else 1
            ax.scatter(x[j], y[j], marker=markers[caseInd], s=100, facecolors='none', edgecolors=colors2[caseInd], linewidths=mew)
        customMarkers = []
        for caseInd, case in enumerate(cases):
            if labelName is None:
                customMarkers.append(Line2D([0],[0],marker=markers[caseInd], color='w', markeredgecolor=colors2[caseInd], markerfacecolor='w', label=case, markersize=15, markeredgewidth=mew))
            else:
                customMarkers.append(Line2D([0],[0],marker=markers[caseInd], color='w', markeredgecolor=colors2[caseInd], markerfacecolor='w', label=labelName[case], markersize=15, markeredgewidth=mew))
        
        #minFlux = -50
        #maxFlux = 50
        #for j in range(0, len(f)):
        #    c = (diff[j]-minFlux)/(maxFlux-minFlux)
        #    ax.scatter(x[j], y[j], s=100, color=cmap(c))
        #customMarkers = []
        #for i, f in enumerate([-50,-25,0,25,50]):
        #    v = (f-minFlux)/(maxFlux-minFlux)
        #    customMarkers.append(Line2D([0],[0],marker=markers[0], color='w', markerfacecolor=cmap(v), label='%0.0f $\mathrm{kW/m^{2}}$'%(f), markersize=15))
        
        ax.legend(handles=customMarkers, fontsize=fs)
    else:
        ax.scatter(x, y, s=100)
    '''
    for j in range(0, len(f)):
        try:
            ind = np.where(f[j] == uniqueFluxes)[0][0]
            c = (plotFlux[j]-minFlux)/(maxFlux-minFlux)
            print(i, j, c, plotFlux[j], x[j], y[j])
            #ax.scatter(x[j], y[j], marker=markers[ind], label=material, s=100, color=cmap(c))
            ax.scatter(x[j], y[j], marker=markers[ind], label=material, s=100, color=colors2[cinds[plotFlux[j]]])
        except:
            pass
    
    customMarkers = []
    for i, f in enumerate(uniqueFluxes):
        v = (f-minFlux)/(maxFlux-minFlux)
        customMarkers.append(Line2D([0],[0],marker=markers[i], color='w', markerfacecolor=colors2[cinds[f]], label='%0.0f $\mathrm{kW/m^{2}}$'%(f), markersize=15))
    
    ax.legend(handles=customMarkers, fontsize=fs)
    '''
    plt.xlabel('Measured %s'%(label), fontsize=fs)
    plt.ylabel('Predicted %s'%(label), fontsize=fs)
    
    plt.xlim([axmin2, axmax2])
    plt.ylim([axmin2, axmax2])
    plt.tick_params(labelsize=fs)
    plt.tick_params(which='major', length=16, width=1, direction='in', top=True,right=True)
    plt.tick_params(which='minor', length=8, width=1, direction='in', top=True,right=True)
    
    #annotation = '%s\n'%(label)
    annotation = ''
    annotation = '%s%s %0.2f\n'%(annotation, 'Exp. Rel. Std. Dev.:', sigma_e)
    annotation = '%s%s %0.2f\n'%(annotation, 'Model Rel. Std. Dev.:', sigma_m)
    annotation = '%s%s %0.2f\n'%(annotation, 'Model Bias Factor:', delta)
    plt.annotate(annotation, (0.5, 0.1), size=fs, xycoords='figure fraction', textcoords='figure fraction', xytext=(0.56,0.1))
    plt.tight_layout()
    return fig, sigma_m, delta


def getNormalStats(delta, sigma_m):
    mu = 1/delta
    sig = ((sigma_m**2)*(mu**2))**0.5
    var = sig ** 2
    
    return mu, var

def preprocessConeData(time, hrrpua, truncateTime=False, truncateBelow=10, 
                       energyThreshold=0.01, filterWidth=101):
    if truncateTime is not False:  hrrpua[time > truncateTime, :] = 0
    if truncateBelow is not False: hrrpua[hrrpua < truncateBelow] = 0
    NT = int(1+np.ceil(time.max())/0.1)
    time_interp = np.linspace(0, np.ceil(time.max()), NT)
    hrrpuas_interp = np.interp(time_interp, time, hrrpua)
    fil = np.ones(filterWidth)/filterWidth
    hrrpuas_interp_filtered = np.convolve(hrrpuas_interp, fil, mode='same')
    
    inds = np.where(~np.isnan(hrrpuas_interp_filtered))
    totalEnergy = np.trapz(hrrpuas_interp_filtered[inds], time_interp[inds])
    
    totalEnergy2 = 0
    ind2 = 1
    while totalEnergy2 < energyThreshold * totalEnergy:
        ind2 = ind2 + 1
        totalEnergy2 = np.trapz(hrrpuas_interp_filtered[:ind2], time_interp[:ind2])
    tign = time_interp[ind2-1]
    hrrpuas_interp_filtered[time_interp < tign] = 0
    return time_interp, hrrpuas_interp_filtered

def getMaterialsDatabase(systemPath):
    
    from git import Repo  # pip install gitpython
    repo = Repo.clone_from('https://github.com/johodges/materials.git', os.path.join(systemPath, "materials"), branch="main")
    
    #cloned_repo = repo.clone('https://github.com/johodges/materials.git')
    repo.submodule_update(recursive=True)

def adjust_tmax_qmax_by_material(material, tMax_by_thickness1, qMax_by_thickness1):
    tMax_by_thickness = dict(tMax_by_thickness1)
    qMax_by_thickness = dict(qMax_by_thickness1)
    initial_length = len(tMax_by_thickness.keys())
    if material == 'Aalto_Pine_Flaming':
        thickness = 20.00000000
        tMax_by_thickness[thickness] = 30
        qMax_by_thickness[thickness] = 250
    elif material == 'Aalto_Spruce_Flaming':
        thickness = 20.00000000
        tMax_by_thickness[thickness] = 30
        qMax_by_thickness[thickness] = 250
    
    elif material == 'FAA_HDPE':
        thickness = 3.2
        tMax_by_thickness[thickness] = 12.0
        qMax_by_thickness[thickness] = 2350.0
        thickness = 8.1
        tMax_by_thickness[thickness] = 28.0
        qMax_by_thickness[thickness] = 2550.0
        thickness = 27.0
        tMax_by_thickness[thickness] = 68.0
        qMax_by_thickness[thickness] = 2800.0
    elif material == 'FAA_HIPS':
        thickness = 3.2
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 1370.0
        thickness = 8.1
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 1425.0
        thickness = 27.0
        tMax_by_thickness[thickness] = 54.0
        qMax_by_thickness[thickness] = 1370.0
    elif material == 'FAA_PBT':
        thickness = 4.00000000
        tMax_by_thickness[thickness] = 7
        qMax_by_thickness[thickness] = 1975.0
    elif material == 'FAA_PBTGF':
        thickness = 4.00000000
        tMax_by_thickness[thickness] = 7.0
        qMax_by_thickness[thickness] = 985.0
    elif material == 'FAA_PC':
        thickness = 3.0
        tMax_by_thickness[thickness] = 5.0
        qMax_by_thickness[thickness] = 950.0
        thickness = 5.5
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 765.0
        thickness = 9.0
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 655.0
    elif material == 'FAA_PEEK':
        thickness = 3.90000000
        tMax_by_thickness[thickness] = 8.0
        qMax_by_thickness[thickness] = 600.0
    elif material == 'FAA_PMMA':
        thickness = 3.2
        tMax_by_thickness[thickness] = 8.0
        qMax_by_thickness[thickness] = 1400.0
        thickness = 8.1
        tMax_by_thickness[thickness] = 17.0
        qMax_by_thickness[thickness] = 1650.0
        thickness = 27.0
        tMax_by_thickness[thickness] = 50.0
        qMax_by_thickness[thickness] = 1600.0
    elif material == 'FAA_PVC':
        thickness = 3.0
        tMax_by_thickness[thickness] = 6.0
        qMax_by_thickness[thickness] = 350.0
        thickness = 6.0
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 400.0
        thickness = 9.0
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 325.0


    elif material == 'FPL_hardboard_6mm':
        thickness = 7.25000000
        tMax_by_thickness[thickness] = 15
        qMax_by_thickness[thickness] = 655.0
    elif material == 'FPL_lumber_redoak_20mm':
        thickness = 19.75000000
        tMax_by_thickness[thickness] = 66.0
        qMax_by_thickness[thickness] = 380.0
    elif material == 'FPL_osb_12mm':
        thickness = 11.50000000
        tMax_by_thickness[thickness] = 25
        qMax_by_thickness[thickness] = 490.0
    elif material == 'FPL_plywood_douglas_fir_12mm':
        thickness = 11.75000000
        tMax_by_thickness[thickness] = 30
        qMax_by_thickness[thickness] = 270.0
    elif material == 'FPL_plywood_douglas_fir_frt_12mm':
        thickness = 12.50000000
        tMax_by_thickness[thickness] = 20
        qMax_by_thickness[thickness] = 270.0
    elif material == 'FPL_plywood_oak_13mm':
        thickness = 12.70000000
        tMax_by_thickness[thickness] = 35
        qMax_by_thickness[thickness] = 325.0
    elif material == 'FPL_plywood_southern_pine_frt_11mm':
        thickness = 11.25000000
        tMax_by_thickness[thickness] = 30
        qMax_by_thickness[thickness] = 270.0
    elif material == 'FPL_waferboard_13mm':
        thickness = 13.00000000
        tMax_by_thickness[thickness] = 25
        qMax_by_thickness[thickness] = 380.0

    elif material == 'FSRI_Asphalt_Shingle':
        thickness = 3.0033
        tMax_by_thickness[thickness] = 5
        qMax_by_thickness[thickness] = 1040.0
    elif material == 'FSRI_Cellulose_Insulation':
        thickness = 37.7206
        tMax_by_thickness[thickness] = 5
        qMax_by_thickness[thickness] = 350.0
    elif material == 'FSRI_Cotton_Rug':
        thickness = 6.1034
        tMax_by_thickness[thickness] = 4
        qMax_by_thickness[thickness] = 930.0
    elif material == 'FSRI_Cotton_Sheet':
        thickness = 0.2201
        tMax_by_thickness[thickness] = 1.0
        qMax_by_thickness[thickness] = 270.0
    elif material == 'FSRI_EPDM_Membrane':
        thickness = 7.7797
        tMax_by_thickness[thickness] = 15
        qMax_by_thickness[thickness] = 930.0
    elif material == 'FSRI_Excelsior':
        thickness = 2.0204
        tMax_by_thickness[thickness] = 1
        qMax_by_thickness[thickness] = 270.0
    elif material == 'FSRI_FDNY_LDF':
        thickness = 12.5935
        tMax_by_thickness[thickness] = 10
        qMax_by_thickness[thickness] = 400.0
    elif material == 'FSRI_FRP_Panel':
        thickness = 6.3791
        tMax_by_thickness[thickness] = 4
        qMax_by_thickness[thickness] = 1260.0
    

    elif material == 'FSRI_Face_Shield':
        thickness = 19.2521
        tMax_by_thickness[thickness] = 40.0
        qMax_by_thickness[thickness] = 435.0
    elif material == 'FSRI_Feather_Pillow_Feathers':
        thickness = 0.8007
        tMax_by_thickness[thickness] = 1.5
        qMax_by_thickness[thickness] = 490.0
    elif material == 'FSRI_Fiberglass_Insulation_R30_Paper_Faced':
        thickness = 2.0204
        tMax_by_thickness[thickness] = 1
        qMax_by_thickness[thickness] = 270.0
    elif material == 'FSRI_Gypsum_Wallboard':
        thickness = 13.0028
        tMax_by_thickness[thickness] = 7
        qMax_by_thickness[thickness] = 325.0
    elif material == 'FSRI_Hemp_Sheet':
        thickness = 0.3804
        tMax_by_thickness[thickness] = 1
        qMax_by_thickness[thickness] = 435.0
    elif material == 'FSRI_House_Wrap':
        thickness = 0.1110
        tMax_by_thickness[thickness] = 5
        qMax_by_thickness[thickness] = 80.0
    elif material == 'FSRI_Latex_Pillow_Foam':
        thickness = 3.2577
        tMax_by_thickness[thickness] = 5
        qMax_by_thickness[thickness] = 1550.0
    elif material == 'FSRI_Lightweight_Gypsum_Wallboard':
        thickness = 12.9049
        tMax_by_thickness[thickness] = 6
        qMax_by_thickness[thickness] = 325.0
        
    elif material == 'FSRI_Overstuffed_Chair_Assembly':
        thickness = 4.2084
        tMax_by_thickness[thickness] = 6
        qMax_by_thickness[thickness] = 1250.0
    elif material == 'FSRI_PE_Foam_Pipe_Insulation':
        thickness = 14.2924
        tMax_by_thickness[thickness] = 4
        qMax_by_thickness[thickness] = 950.0
    elif material == 'FSRI_Polyisocyanurate_Foam_Board':
        thickness = 13.7129
        tMax_by_thickness[thickness] = 4
        qMax_by_thickness[thickness] = 380.0
    elif material == 'FSRI_Pressure_Treated_Deck':
        thickness = 19.2521
        tMax_by_thickness[thickness] = 40.0
        qMax_by_thickness[thickness] = 435.0
    elif material == 'FSRI_Roof_Felt':
        thickness = 1.0731
        tMax_by_thickness[thickness] = 3.0
        qMax_by_thickness[thickness] = 1205.0
    elif material == 'FSRI_Rubber_Band':
        thickness = 4.6658
        tMax_by_thickness[thickness] = 6.0
        qMax_by_thickness[thickness] = 1040.0
    elif material == 'FSRI_Rubber_Foam_Pipe_Insulation':
        thickness = 12.7742
        tMax_by_thickness[thickness] = 3.0
        qMax_by_thickness[thickness] = 380.0
    elif material == 'FSRI_Rug_Pad':
        thickness = 4.0605
        tMax_by_thickness[thickness] = 3.0
        qMax_by_thickness[thickness] = 600.0
        
    elif material == 'FSRI_Wool_Rug':
        thickness = 16.6327
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 550.0
    elif material == 'FSRI_XPS_Foam_Board':
        thickness = 27.0723
        tMax_by_thickness[thickness] = 5.0
        qMax_by_thickness[thickness] = 1450.0
        
    elif material == 'FSRI_ABS':
        thickness = 2.9545
        tMax_by_thickness[thickness] = 6
        qMax_by_thickness[thickness] = 2150.0
    elif material == 'FSRI_Black_PMMA':
        thickness = 8.6925
        tMax_by_thickness[thickness] = 13.0
        qMax_by_thickness[thickness] = 1700.0
    elif material == 'FSRI_Cotton_Raw':
        thickness = 2.4410
        tMax_by_thickness[thickness] = 3
        qMax_by_thickness[thickness] = 710.0
    elif material == 'FSRI_HDPE':
        thickness = 3.2414
        tMax_by_thickness[thickness] = 10
        qMax_by_thickness[thickness] = 2965.0
    elif material == 'FSRI_HIPS':
        thickness = 2.9531
        tMax_by_thickness[thickness] = 7.0
        qMax_by_thickness[thickness] = 1920.0
    elif material == 'FSRI_High_Temperature_SCBA_Facepiece':
        thickness = 37.1994
        tMax_by_thickness[thickness] = 152.0
        qMax_by_thickness[thickness] = 350.0
    elif material == 'FSRI_LDPE':
        thickness = 3.2459
        tMax_by_thickness[thickness] = 7.0
        qMax_by_thickness[thickness] = 3075.0
    elif material == 'FSRI_Memory_Foam_Carpet_Pad':
        thickness = 12.4020
        tMax_by_thickness[thickness] = 3.0
        qMax_by_thickness[thickness] = 1850.0
    
    elif material == 'FSRI_Nylon':
        thickness = 3.4008
        tMax_by_thickness[thickness] = 16
        qMax_by_thickness[thickness] = 1950.0
    elif material == 'FSRI_Nylon_Carpet_High_Pile':
        thickness = 13.7257
        tMax_by_thickness[thickness] = 18
        qMax_by_thickness[thickness] = 900.0
    elif material == 'FSRI_Overstuffed_Chair_Polyester_Batting':
        thickness = 1.4901
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 710.0
    elif material == 'FSRI_Overstuffed_Chair_Polyester_Fabric':
        thickness = 0.5245
        tMax_by_thickness[thickness] = 11.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'FSRI_Overstuffed_Chair_Polyurethane_Foam':
        thickness = 0.8163
        tMax_by_thickness[thickness] = 2
        qMax_by_thickness[thickness] = 1535.0
    elif material == 'FSRI_PC':
        thickness = 5.3279
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 820.0
    elif material == 'FSRI_PET':
        thickness = 6.5311
        tMax_by_thickness[thickness] = 17.0
        qMax_by_thickness[thickness] = 1315.0
    elif material == 'FSRI_PETG':
        thickness = 2.6256
        tMax_by_thickness[thickness] = 8.0
        qMax_by_thickness[thickness] = 2745.0
        
    elif material == 'FSRI_PMMA':
        thickness = 2.8173
        tMax_by_thickness[thickness] = 6
        qMax_by_thickness[thickness] = 1920.0
    elif material == 'FSRI_PP':
        thickness = 3.2256
        tMax_by_thickness[thickness] = 7.0
        qMax_by_thickness[thickness] = 3130.0
    elif material == 'FSRI_PVC':
        thickness = 3.1772
        tMax_by_thickness[thickness] = 28.0
        qMax_by_thickness[thickness] = 400.0
    elif material == 'FSRI_Pallet_Wood':
        thickness = 2.9531
        tMax_by_thickness[thickness] = 7
        qMax_by_thickness[thickness] = 1920.0
    elif material == 'FSRI_PlasticC':
        thickness = 2.9782
        tMax_by_thickness[thickness] = 7
        qMax_by_thickness[thickness] = 1975.0
    elif material == 'FSRI_Plastic_Laminate_Countertop':
        thickness = 37.1994
        tMax_by_thickness[thickness] = 152.0
        qMax_by_thickness[thickness] = 350.0
    elif material == 'FSRI_Plywood':
        thickness = 8.6925
        tMax_by_thickness[thickness] = 13
        qMax_by_thickness[thickness] = 1700.0
    elif material == 'FSRI_Polyester_Bed_Skirt':
        thickness = 1.1950
        tMax_by_thickness[thickness] = 8
        qMax_by_thickness[thickness] = 655.0
        
    elif material == 'FSRI_Polyester_Microfiber_Sheet':
        thickness = 1.0386
        tMax_by_thickness[thickness] = 10
        qMax_by_thickness[thickness] = 655.0
    elif material == 'FSRI_Polyolefin_Carpet_Low_Pile':
        thickness = 7.1574
        tMax_by_thickness[thickness] = 5
        qMax_by_thickness[thickness] = 1150.0
    elif material == 'FSRI_Rebond_Foam_Carpet_Pad':
        thickness = 9.0117
        tMax_by_thickness[thickness] = 2
        qMax_by_thickness[thickness] = 1650.0
    elif material == 'FSRI_Vinyl_Plank_Flooring':
        thickness = 2.4410
        tMax_by_thickness[thickness] = 3
        qMax_by_thickness[thickness] = 710.0
    elif material == 'FSRI_Vinyl_Siding':
        thickness = 1.1668
        tMax_by_thickness[thickness] = 8.0
        qMax_by_thickness[thickness] = 435.0
    elif material == 'FSRI_Vinyl_Tile':
        thickness = 7.8814
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 380.0
        
        
    elif material == 'FSRI_Basswood_Panel':
        thickness = 19.8317
        tMax_by_thickness[thickness] = 30.0
        qMax_by_thickness[thickness] = 490.0
    elif material == 'FSRI_Composite_Deck_Board':
        thickness = 13.2801
        tMax_by_thickness[thickness] = 58.0
        qMax_by_thickness[thickness] = 1350.0
    elif material == 'FSRI_Engineered_Flooring':
        thickness = 9.0160
        tMax_by_thickness[thickness] = 20
        qMax_by_thickness[thickness] = 600.0
    elif material == 'FSRI_Engineered_Wood_Furniture':
        thickness = 12.0809
        tMax_by_thickness[thickness] = 30
        qMax_by_thickness[thickness] = 710.0
    elif material == 'FSRI_Engineered_Wood_Table':
        thickness = 89.839
        tMax_by_thickness[thickness] = 60
        qMax_by_thickness[thickness] = 500.0
    elif material == 'FSRI_Eucalyptus_Flooring':
        thickness = 15.6166
        tMax_by_thickness[thickness] = 50
        qMax_by_thickness[thickness] = 710.0
    elif material == 'FSRI_Homasote':
        thickness = 13.2714
        tMax_by_thickness[thickness] = 20
        qMax_by_thickness[thickness] = 400.0
    elif material == 'FSRI_Luan_Panel':
        thickness = 5.8665
        tMax_by_thickness[thickness] = 6
        qMax_by_thickness[thickness] = 490.0

    elif material == 'FSRI_MDF':
        thickness = 19.4120
        tMax_by_thickness[thickness] = 40.0
        qMax_by_thickness[thickness] = 600.0
    elif material == 'FSRI_Masonite_Board':
        thickness = 3.0701
        tMax_by_thickness[thickness] = 8
        qMax_by_thickness[thickness] = 1040.0
    elif material == 'FSRI_OSB':
        thickness = 16.9322
        tMax_by_thickness[thickness] = 32.0
        qMax_by_thickness[thickness] = 490.0
    elif material == 'FSRI_Oak_Flooring':
        thickness = 19.8678
        tMax_by_thickness[thickness] = 40
        qMax_by_thickness[thickness] = 545.0
    elif material == 'FSRI_Particleboard':
        thickness = 20.3250
        tMax_by_thickness[thickness] = 52.0
        qMax_by_thickness[thickness] = 550.0
    elif material == 'FSRI_Pine_Siding':
        thickness = 18.8606
        tMax_by_thickness[thickness] = 34.0
        qMax_by_thickness[thickness] = 490.0
    elif material == 'FSRI_Wood_Stud':
        thickness = 46.1926
        tMax_by_thickness[thickness] = 88.0
        qMax_by_thickness[thickness] = 490.0

    elif material == 'JH_Acrylic':
        thickness = 4.50000000
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 1400.0
    elif material == 'JH_Black PMMA':
        thickness = 9.20000000
        tMax_by_thickness[thickness] = 13.0
        qMax_by_thickness[thickness] = 1535.0
    elif material == 'JH_Cardboard':
        thickness = 4.10000000
        tMax_by_thickness[thickness] = 2.0
        qMax_by_thickness[thickness] = 450.0
    elif material == 'JH_CPS Balsa Facesheet':
        thickness = 15.90000000
        tMax_by_thickness[thickness] = 34.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'JH_CPS Plywood Facesheet':
        thickness = 12.70000000
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 400.0
    elif material == 'JH_FRP':
        thickness = 12.70000000
        tMax_by_thickness[thickness] = 50.0
        qMax_by_thickness[thickness] = 435.0
    elif material == 'JH_MDF':
        thickness = 19.20000000
        tMax_by_thickness[thickness] = 38.0
        qMax_by_thickness[thickness] = 545.0
    elif material == 'JH_OSB':
        thickness = 16.10000000
        tMax_by_thickness[thickness] = 26.0
        qMax_by_thickness[thickness] = 435.0
    elif material == 'JH_PC Blend':
        thickness = 4.50000000
        tMax_by_thickness[thickness] = 16.0
        qMax_by_thickness[thickness] = 350.0
    elif material == 'JH_Phenolic Resin Fiberglass Composite':
        thickness = 3.27000000
        tMax_by_thickness[thickness] = 30.0
        qMax_by_thickness[thickness] = 215.0
    elif material == 'JH_Plywood':
        thickness = 6.35000000
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 600.0
    elif material == 'JH_PVC Blend':
        thickness = 3.27000000
        tMax_by_thickness[thickness] = 22.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'JH_Vinyl Ester Resin FRP':
        thickness = 4.50000000
        tMax_by_thickness[thickness] = 36.0
        qMax_by_thickness[thickness] = 270.0
    elif material == 'JH_White Pine':
        thickness = 19.10000000
        tMax_by_thickness[thickness] = 30.0
        qMax_by_thickness[thickness] = 400.0
    elif material == 'JH_White Spruce':
        thickness = 37.22000000
        tMax_by_thickness[thickness] = 74.0
        qMax_by_thickness[thickness] = 325.0

    elif material == 'RISE_80_wool__20__Nylon_Glue_Plywood-22m':
        thickness = 22.00000000
        tMax_by_thickness[thickness] = 32.0
        qMax_by_thickness[thickness] = 550.0
    elif material == 'RISE_Fabric_Foam-28mm':
        thickness = 28.00000000
        tMax_by_thickness[thickness] = 15
        qMax_by_thickness[thickness] = 550.0
    elif material == 'RISE_Fabric_Protection_layer_Foam-32mm':
        thickness = 32.00000000
        tMax_by_thickness[thickness] = 20
        qMax_by_thickness[thickness] = 545.0
    elif material == 'RISE_Fabric_vandaliz_protected_Foam-42mm':
        thickness = 42.00000000
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 550.0
    elif material == 'RISE_HPL_Melamine_polyester_film_-13mm':
        thickness = 13.40000000
        tMax_by_thickness[thickness] = 5
        qMax_by_thickness[thickness] = 270.0
    elif material == 'RISE_Melami_face_Calcium_silicat_board-1':
        thickness = 12.50000000
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 250.0
    elif material == 'RISE_Needl_punch_carpe_Glue_Recor_sealin':
        thickness = 10.00000000
        tMax_by_thickness[thickness] = 5.0
        qMax_by_thickness[thickness] = 550.0
    elif material == 'RISE_PE_XLPE-40mm':
        thickness = 40.00000000
        tMax_by_thickness[thickness] = 168.0
        qMax_by_thickness[thickness] = 380.0
        
    elif material == 'RISE_POlyolefin_XLPE-45mm':
        thickness = 45.00000000
        tMax_by_thickness[thickness] = 60.0
        qMax_by_thickness[thickness] = 270.0
    elif material == 'RISE_PUR_rigid_Plasti_faced_steel_sheet-':
        thickness = 79.00000000
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 160.0
    elif material == 'RISE_PVC_EPR-32mm':
        thickness = 32.50000000
        tMax_by_thickness[thickness] = 80.0
        qMax_by_thickness[thickness] = 250.0
    elif material == 'RISE_PVC_PE-00mm':
        thickness = 0.01270000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 650.0
    elif material == 'RISE_PVC_PE-05mm':
        thickness = 5.10000000
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 600.0
    elif material == 'RISE_PVC_PE-08mm':
        thickness = 8.20000000
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 550.0
    elif material == 'RISE_PVC_PE-10mm':
        thickness = 10.20000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 545.0
    elif material == 'RISE_PVC_PE-14mm':
        thickness = 14.80000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 900.0
        
    elif material == 'RISE_PVC_PE-20mm':
        thickness = 20.00000000
        tMax_by_thickness[thickness] = 30.0
        qMax_by_thickness[thickness] = 600.0
    elif material == 'RISE_PVC_XLPE-17mm':
        thickness = 17.70000000
        tMax_by_thickness[thickness] = 25.0
        qMax_by_thickness[thickness] = 710.0
    elif material == 'RISE_PVC_XLPE-22mm':
        thickness = 22.00000000
        tMax_by_thickness[thickness] = 35.0
        qMax_by_thickness[thickness] = 450.0
    elif material == 'RISE_PVC_XLPE-35mm':
        thickness = 35.00000000
        tMax_by_thickness[thickness] = 58.0
        qMax_by_thickness[thickness] = 1300.0
    elif material == 'RISE_PVC_XLPE-38mm':
        thickness = 38.60000000
        tMax_by_thickness[thickness] = 80.0
        qMax_by_thickness[thickness] = 600.0
    elif material == 'RISE_PVC_XLPE-46mm':
        thickness = 46.00000000
        tMax_by_thickness[thickness] = 122.0
        qMax_by_thickness[thickness] = 380.0
    elif material == 'RISE_PVC_wall_carpet_paper_plasterboard-':
        thickness = 0.01270000
        tMax_by_thickness[thickness] = 6.0
        qMax_by_thickness[thickness] = 215.0
    elif material == 'RISE_Painted_paper_plasterboa_plasterboa':
        thickness = 0.01270000
        tMax_by_thickness[thickness] = 3
        qMax_by_thickness[thickness] = 350.0

    elif material == 'RISE_Polyolefin_EPR-18mm':
        thickness = 18.10000000
        tMax_by_thickness[thickness] = 86.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'RISE_Polyolefin_EPR-32mm':
        thickness = 32.20000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'RISE_Polyolefin_PA-02mm':
        thickness = 2.50000000
        tMax_by_thickness[thickness] = 8.0
        qMax_by_thickness[thickness] = 435.0
    elif material == 'RISE_Polyolefin_PA-06mm':
        thickness = 6.00000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 400.0
    elif material == 'RISE_Polyolefin_PP-08mm':
        thickness = 8.70000000
        tMax_by_thickness[thickness] = 54.0
        qMax_by_thickness[thickness] = 350.0
    elif material == 'RISE_Polyolefin_XLPE-18mm':
        thickness = 18.10000000
        tMax_by_thickness[thickness] = 76.0
        qMax_by_thickness[thickness] = 435.0
    elif material == 'RISE_Polyolefin_XLPE-25mm':
        thickness = 25.00000000
        tMax_by_thickness[thickness] = 122.0
        qMax_by_thickness[thickness] = 450.0
    elif material == 'RISE_Polyolefin_XLPE-38mm':
        thickness = 38.10000000
        tMax_by_thickness[thickness] = 138.0
        qMax_by_thickness[thickness] = 435.0

    elif material == 'RISE_RPPVC_PEF-04mm':
        thickness = 4.50000000
        tMax_by_thickness[thickness] = 8.0
        qMax_by_thickness[thickness] = 215.0
    elif material == 'RISE_RPPVC_PVC-14mm':
        thickness = 14.00000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'RISE_RPPVC_XLPE-17mm':
        thickness = 17.70000000
        tMax_by_thickness[thickness] = 35.0
        qMax_by_thickness[thickness] = 600.0
    elif material == 'RISE_RPPVC_XLPE-22mm':
        thickness = 22.50000000
        tMax_by_thickness[thickness] = 50.0
        qMax_by_thickness[thickness] = 550.0
    elif material == 'RISE_RPPVC_XLPE-39mm':
        thickness = 39.30000000
        tMax_by_thickness[thickness] = 80.0
        qMax_by_thickness[thickness] = 545.0
    elif material == 'RISE_RPPVC_XLPE-45mm':
        thickness = 45.00000000
        tMax_by_thickness[thickness] = 154.0
        qMax_by_thickness[thickness] = 270.0
    elif material == 'RISE_Synthetic_rubber_Glue_Plywood-15mm':
        thickness = 15.00000000
        tMax_by_thickness[thickness] = 34.0
        qMax_by_thickness[thickness] = 600.0
    elif material == 'RISE_Textile_wall_coverin_paper_plasterb':
        thickness = 0.01270000
        tMax_by_thickness[thickness] = 4.0
        qMax_by_thickness[thickness] = 435.0

    elif material == 'RISE_ZHPolyolefin_PP-08mm':
        thickness = 8.40000000
        tMax_by_thickness[thickness] = 25.0
        qMax_by_thickness[thickness] = 490.0
    elif material == 'RISE_ZHPolyolefin_XLPE-13mm':
        thickness = 13.00000000
        tMax_by_thickness[thickness] = 50.0
        qMax_by_thickness[thickness] = 700.0
    elif material == 'RISE_ZHPolyolefin_XLPE-27mm':
        thickness = 27.00000000
        tMax_by_thickness[thickness] = 70.0
        qMax_by_thickness[thickness] = 450.0
    elif material == 'RISE_fabric__vandaliz_protecte_foam-42mm':
        thickness = 42.00000000
        tMax_by_thickness[thickness] = 6.0
        qMax_by_thickness[thickness] = 300.0
    elif material == 'RISE_synthetic_rubber_glue_plywood-14mm':
        thickness = 14.80000000
        tMax_by_thickness[thickness] = 30.0
        qMax_by_thickness[thickness] = 270.0

    elif material == 'RISE_Alumi_Honey_comb_coated_with_HPL-22':
        thickness = 22.70000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 400.0
    elif material == 'RISE_HPL_compact_-04mm':
        thickness = 4.00000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 300.0
    elif material == 'RISE_Spruce-10mm':
        thickness = 10.00000000
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'RISE_Wool_fabric_Mixed_fabric-00mm':
        thickness = 0.50000000
        tMax_by_thickness[thickness] = 5.0
        qMax_by_thickness[thickness] = 200.0
    elif material == 'RISE_Woolfabric__mixed_fabric-00mm':
        thickness = 0.50000000
        tMax_by_thickness[thickness] = 5.0
        qMax_by_thickness[thickness] = 435.0


    elif material == 'RISE_FR_polycarbonate-16mm':
        thickness = 16.00000000
        tMax_by_thickness[thickness] = 10.0
        qMax_by_thickness[thickness] = 710.0
    elif material == 'RISE_PVC-02mm':
        thickness = 2.90000000
        tMax_by_thickness[thickness] = 5.0
        qMax_by_thickness[thickness] = 380.0
    elif material == 'RISE_PVC_PVC-08mm':
        thickness = 8.20000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'RISE_PVC_PVC-09mm':
        thickness = 9.40000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 380.0
    elif material == 'RISE_PVC_PVC-10mm':
        thickness = 10.00000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 380.0
    elif material == 'RISE_PVC_PVC-14mm':
        thickness = 14.50000000
        tMax_by_thickness[thickness] = 20.0
        qMax_by_thickness[thickness] = 380.0
    elif material == 'RISE_PVC_PVC-18mm':
        thickness = 18.00000000
        tMax_by_thickness[thickness] = 30.0
        qMax_by_thickness[thickness] = 380.0
    elif material == 'RISE_PVC_PVC-21mm':
        thickness = 21.40000000
        tMax_by_thickness[thickness] = 40.0
        qMax_by_thickness[thickness] = 380.0

    elif material == 'RISE_PVDF-01mm':
        thickness = 1.95000000
        tMax_by_thickness[thickness] = 4.0
        qMax_by_thickness[thickness] = 700.0
    elif material == 'RISE_Paint_GFK_polyes_with_Gelcoa_handla':
        thickness = 4.80000000
        tMax_by_thickness[thickness] = 25.0
        qMax_by_thickness[thickness] = 270.0
    elif material == 'RISE_Paint_GRP_polyes_with_gelcoa_handla':
        thickness = 4.80000000
        tMax_by_thickness[thickness] = 12.0
        qMax_by_thickness[thickness] = 160.0
    elif material == 'RISE_Polyester-02mm':
        thickness = 2.10000000
        tMax_by_thickness[thickness] = 5.0
        qMax_by_thickness[thickness] = 490.0

    elif material == 'RISE_Polyolefin-02mm':
        thickness = 2.90000000
        tMax_by_thickness[thickness] = 8.0
        qMax_by_thickness[thickness] = 545.0
    elif material == 'RISE_RPPVC-02mm':
        thickness = 2.90000000
        tMax_by_thickness[thickness] = 6.0
        qMax_by_thickness[thickness] = 325.0
    elif material == 'RISE_Transparent_polycarbonate-02mm':
        thickness = 2.30000000
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 545.0

    elif material == 'RISE_FR_Particle_board-79mm':
        thickness = 79.00000000
        tMax_by_thickness[thickness] = 5.0
        qMax_by_thickness[thickness] = 435.0
    elif material == 'RISE_FR_particle_board-16mm':
        thickness = 16.00000000
        tMax_by_thickness[thickness] = 3.0
        qMax_by_thickness[thickness] = 270.0
    elif material == 'RISE_MDF_board-12mm':
        thickness = 12.00000000
        tMax_by_thickness[thickness] = 15.0
        qMax_by_thickness[thickness] = 490.0
    
    final_length = len(tMax_by_thickness.keys())
    
    if initial_length != final_length:
        print("warning for material %s"%(material))
        print(tMax_by_thickness1.keys())
        print(tMax_by_thickness.keys())
    
    return tMax_by_thickness, qMax_by_thickness 

if __name__ == "__main__":
    
    args = sys.argv
    systemPath = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('call')
    parser.add_argument('--clean', action='store_true', help='Deletes processed data and outputs prior to run')
    parser.add_argument('--donotinitialize', action='store_true', help='Ignores initialization')
    parser.add_argument('--donotrun', action='store_true', help='Do not rerun fds')
    parser.add_argument('--donotcopy', action='store_true', help='Do not copy results to out/exp')
    
    cmdargs = parser.parse_args(args)
    if cmdargs.clean:
        if os.path.exists(os.path.join(systemPath,'materials')):
            print("Cleaning materials repo from %s"%(os.path.join(systemPath, 'materials')))
            shutil.rmtree(os.path.join(systemPath,'materials'))
        for f in glob.glob(os.path.join(systemPath,'..','*.fds')):
            os.remove(f)
    
    if os.path.exists(os.path.join(systemPath,'materials')) is not True:
        print("Cloning materials repo to %s"%(os.path.join(systemPath, 'materials')))
        getMaterialsDatabase(systemPath)
    
    assert os.path.exists(os.path.join(systemPath,'materials')), "Failed to clone materials repo"
    
    my_env = os.environ.copy()
    if cmdargs.donotinitialize is not True:
        process = subprocess.Popen([sys.executable, os.path.join(systemPath, 'materials','scripts','initialize.py')], env=my_env, shell=False)
        out, err = process.communicate()
        errcode = process.returncode
    else:
        print("Ignoring initialization")
    
    runcommand = [sys.executable, os.path.join(systemPath, 'materials','scripts','evaluate_database.py'),'--inputfiles',os.path.join(systemPath,'..')]
    if cmdargs.donotrun is True:
        runcommand.append('--donotrun')
    process = subprocess.Popen(runcommand, env=my_env, shell=False)
    out, err = process.communicate()
    errcode = process.returncode
    
    if cmdargs.donotcopy is False:
        #files = glob.glob(os.path.join(systemPath,'materials','input_files','*','*_devc.csv'))
        #for file in files:
        #    shutil.copy(file, os.path.join(systemPath,'..','..','..','out','Scaling_Pyrolysis',os.path.basename(file)))
        #files = glob.glob(os.path.join(systemPath,'materials','input_files','*','*_git.txt'))
        #for file in files:
        #    shutil.copy(file, os.path.join(systemPath,'..','..','..','out','Scaling_Pyrolysis',os.path.basename(file)))
        files = glob.glob(os.path.join(systemPath,'materials','data','aalto_materials','*.csv'))
        for file in files:
            ignores = ['-29mm-15.csv','-29mm-20.csv','-29mm-25.csv','-29mm-35.csv','-29mm-50.csv','-29mm-75.csv']
            if os.path.basename not in ignores:
                shutil.copy(file, os.path.join(systemPath,'..','..','..','..','..','exp','Scaling_Pyrolysis',os.path.basename(file)))

        files = glob.glob(os.path.join(systemPath,'materials','data','fpl_materials_processed','*.csv'))
        for file in files:
            ignores = ['-29mm-15.csv','-29mm-20.csv','-29mm-25.csv','-29mm-35.csv','-29mm-50.csv','-29mm-75.csv']
            if os.path.basename not in ignores:
                shutil.copy(file, os.path.join(systemPath,'..','..','..','..','..','exp','Scaling_Pyrolysis',os.path.basename(file)))
        
        files = glob.glob(os.path.join(systemPath,'materials','data','fsri_materials_processed','scaling_pyrolysis','*.csv'))
        for file in files:
            ignores = ['-29mm-15.csv','-29mm-20.csv','-29mm-25.csv','-29mm-35.csv','-29mm-50.csv','-29mm-75.csv']
            if os.path.basename not in ignores:
                shutil.copy(file, os.path.join(systemPath,'..','..','..','..','..','exp','Scaling_Pyrolysis',os.path.basename(file)))
        
        files = glob.glob(os.path.join(systemPath,'materials','data','jh_materials','*.csv'))
        for file in files:
            ignores = ['-29mm-15.csv','-29mm-20.csv','-29mm-25.csv','-29mm-35.csv','-29mm-50.csv','-29mm-75.csv']
            if os.path.basename not in ignores:
                shutil.copy(file, os.path.join(systemPath,'..','..','..','..','..','exp','Scaling_Pyrolysis',os.path.basename(file)))
        
        files = glob.glob(os.path.join(systemPath,'materials','data','rise_materials_processed','*.csv'))
        for file in files:
            ignores = ['-29mm-15.csv','-29mm-20.csv','-29mm-25.csv','-29mm-35.csv','-29mm-50.csv','-29mm-75.csv']
            if os.path.basename(file) not in ignores:
                shutil.copy(file, os.path.join(systemPath,'..','..','..','..','..','exp','Scaling_Pyrolysis',os.path.basename(file)))
        
    from materials.scripts.algorithms import load_material_simulation, get_filtered_cases, getMaterials, load_csv, get_ignition_temperature
    #from materials.scripts.evaluate_database import load_material_simulation, get_filtered_cases
    spec_file_dict = getMaterials(dataDirectory=os.path.abspath(os.path.join(systemPath,'materials','data')))
    materials = list(spec_file_dict.keys())
    baseDir = os.path.join(systemPath,'materials','input_files')
    energyThreshold=0.0
    output_statistics = dict()
    material_output_data = defaultdict(bool)
    outTxt = ''
    #materials = ['FAA_PMMA']
    completeMaterials = []
    outTxt3 = ''
    for material in materials:
        cases = get_filtered_cases(spec_file_dict, material, energyThreshold=energyThreshold)
        if cases is False: continue
        if len(list(cases.keys())) <= 1: continue
        completeMaterials.append(material)
        output_statistics[material] = load_material_simulation(material, baseDir, cases)
        Tign = get_ignition_temperature(material, baseDir, cases)
        spec_file_dict[material]['Tign'] = Tign
        
        materialClass = spec_file_dict[material]['materialClass']
        series = spec_file_dict[material]['series']
        properties = spec_file_dict[material]
        density = properties['density']
        conductivity = properties['conductivity']
        emissivity = properties['emissivity']
        specific_heat = properties['specific_heat']
        soot_yield = properties['soot_yield']
        heat_of_combustion = properties['heat_of_combustion']
        fluxes = [cases[c]['cone'] for c in cases.keys()]
        thicknesses = [cases[c]['delta']*1e3 for c in cases.keys()]
        material_output_data[(materialClass,series,material)] = defaultdict(bool)
        material_output_data[(materialClass,series,material)]['Conductivity\n($\mathrm{W/(m\cdot K)}$)'] = conductivity
        material_output_data[(materialClass,series,material)]['Density\n($\mathrm{kg/m^{3}}$)'] = density
        material_output_data[(materialClass,series,material)]['Emissivity\n(-)'] = emissivity
        material_output_data[(materialClass,series,material)]['Specific Heat\n($\mathrm{kJ/(kg\cdot C}$)'] = specific_heat
        material_output_data[(materialClass,series,material)]['Thickness\n($\mathrm{mm}$)'] = '|'.join(['%0.1f'%(cases[c]['delta']*1000) for c in cases.keys()])
        #material_output_data[(materialClass,series,material)]['Ignition Temperaure\n($\mathrm{^{\circ}C}$)'] = Tign
        #material_output_data[(materialClass,series,material)]['Reference Heat Flux\n($\mathrm{kW/m^{2}}$)'] = qref
        #material_output_data[(materialClass,series,material)]['Validation Heat Fluxes\n($\mathrm{kW/m^{2}}$)'] = '|'.join(['"HRRPUA-%02d"'%(flux) for flux in validationFluxes])
        material_output_data[(materialClass,series,material)]['Validation Heat Fluxes\n($\mathrm{kW/m^{2}}$)'] = '|'.join(['%0.0f'%(f) for f in fluxes])
        
        lineColors = ['k','r','g','m','c','y','k','r','g','m','c','y']
        materialClassDict = dict()
        materialClassDict['Wood-Based'] = dict()
        materialClassDict['Wood-Based']['marker'] = '>'
        materialClassDict['Wood-Based']['color'] = 'b'
        materialClassDict['Polymers'] = dict()
        materialClassDict['Polymers']['marker'] = '^'
        materialClassDict['Polymers']['color'] = 'r'
        materialClassDict['Mixtures'] = dict()
        materialClassDict['Mixtures']['marker'] = 'v'
        materialClassDict['Mixtures']['color'] = 'g'
        materialClassDict['Others'] = dict()
        materialClassDict['Others']['marker'] = 'o'
        materialClassDict['Others']['color'] = 'm'
        
        workingDir = os.path.join(baseDir, material.replace(' ','_')) + os.sep
        chid = material.replace(' ','_')
        data = load_csv(workingDir, chid)
        qMax = np.ceil(np.nanmax(data[[c for c in data.columns if 'HRRPUA' in c]].values)/10*1.1)*10
        tInd = np.argwhere(np.sum(data[[c for c in data.columns if 'HRRPUA' in c]].values,axis=1)>0)[-1][0]
        
        tMax_by_thickness = dict()
        qMax_by_thickness = dict()
        for thickness in sorted(list(set(thicknesses))):
            columns = [c for c in data.columns if 'HRRPUA' in c and ('%0.2f'%(thickness)).replace('.','p') in c]
            qMax = (np.ceil(np.nanmax(data[columns])/50)*1.1+1)*50
            thickness_rounded = np.round(thickness, decimals=4)
            qMax_by_thickness[thickness_rounded] = max([qMax, 0])
            tind = np.where(np.sum(abs(data[columns]),axis=1) > 0)[0][-1]
            tMax = (np.ceil(data['Time'].values[tind]/60/2)+1)*2
            tMax_by_thickness[thickness_rounded] = max([tMax, 2])
        
        '''
        tMax = np.ceil(data['Time'].iloc[tInd]/60/5)*5
        qMax = max([qMax, 0])
        tMax = max([tMax, 5])
        '''
        expFiles = spec_file_dict[material]['expFiles']
        headerRows = spec_file_dict[material]['headerRows']
        if type(headerRows) is int:
            headerRows = [headerRows for x in expFiles]
        else:
            headerRows = [int(x) for x in headerRows.split('|')]
        for i, expFile in enumerate(expFiles):
            expf = os.path.join(systemPath,'..','..','..','..','..','exp','Scaling_Pyrolysis',expFile)
            headerRow = headerRows[i]
            with open(expf, 'r') as f:
                d = f.readlines()
            d = np.array([dd.replace('\n','').split(',') for dd in d])
            
            for ii in range(headerRow, len(d)):
                for j in range(0, len(d[ii])):
                    try:
                        d[ii,j] = float(d[ii,j])
                    except:
                        d[ii,j] = np.nan
            columns = [str(c).replace('/','_') for c in d[0]]
            data_exp = pd.DataFrame(np.array(d[headerRow:, :], dtype=float), columns=columns)
            hrrColumn = spec_file_dict[material]['hrrColumns'][i]
            if '.csv' in hrrColumn: hrrColumn = hrrColumn.split('.csv-')[1]
            qMax = np.ceil(np.nanmax(data_exp[hrrColumn].values)*1.1/50+1)*50
            timeColumn = spec_file_dict[material]['timeColumns'][i]
            if '.csv' in timeColumn: timeColumn = timeColumn.split('.csv-')[1]
            tMax = np.ceil(np.nanmax(data_exp[timeColumn].values)/60/2+1)*2
            thickness = thicknesses[i]
            thickness_rounded = np.round(thickness, decimals=4)
            tMax_by_thickness[thickness_rounded] = max([tMax, tMax_by_thickness[thickness_rounded]])
            qMax_by_thickness[thickness_rounded] = max([qMax, qMax_by_thickness[thickness_rounded]])
        tMax_by_thickness, qMax_by_thickness = adjust_tmax_qmax_by_material(material, tMax_by_thickness, qMax_by_thickness)
        for thickness in sorted(list(set(thicknesses))):
            counter = 0
            thickness_rounded = np.round(thickness, decimals=4)
            for iii, flux in enumerate(fluxes):
                if thickness != thicknesses[iii]: continue
                lc = lineColors[counter]
                if counter == 0:
                    switchId = 'd'
                else:
                    switchId = 'f'
                counter += 1
                materialMarker = materialClassDict[materialClass]['marker']
                materialColor = materialClassDict[materialClass]['color']
                matlabelname = material
                while '__' in matlabelname: matlabelname = matlabelname.replace('__','_')
                matlabelname = '_'.join(matlabelname.split('_')[1:])
                if 'FPL' in series:
                    matlabelname = '_'.join(matlabelname.split('_')[:-1])
                matlabelname = matlabelname.replace('_','\_')
                outTxt = outTxt + switchId + ',' + 'Scaling_Pyrolysis_'+materialClass + ',Scaling_Pyrolysis/' + expFiles[iii] + ',' #'%s-%02d.csv"'%(material,flux) + ","
                timeColumn,hrrColumn=spec_file_dict[material]['timeColumns'][iii],spec_file_dict[material]['hrrColumns'][iii]
                if '.csv' in timeColumn: timeColumn = timeColumn.split('.csv-')[1]
                if '.csv' in hrrColumn: hrrColumn = hrrColumn.split('.csv-')[1]
                outTxt = outTxt + "1,2,%s,%s,"%(timeColumn, hrrColumn)
                outTxt = outTxt + "Exp (%02d kW/m²),%s-,0,100000,,0,100000,-10000,10000,0,"%(flux, lc)
                hrrColumn = ('HRRPUA-CONE_%03.2f_%03d'%(thicknesses[iii], flux)).replace('.','p') 
                outTxt = outTxt + 'Scaling_Pyrolysis/' + chid + '_devc.csv,2,3,Time,' + hrrColumn + ',' + 'FDS (%02d kW/m²),%s--,0,100000,,0,100000,-10000,10000,0,'%(flux, lc)
                outTxt = outTxt + '%s %0.1f mm,Time (min),Heat Release Rate (kW/m²),'%(matlabelname, thickness)
                outTxt = outTxt + '0,%0.0f,60,0,%0.0f,1,no,0.05 0.90,NorthEast,,1.0,Scaling_Pyrolysis/'%(tMax_by_thickness[thickness_rounded], qMax_by_thickness[thickness_rounded]) + chid +'_git.txt,linear,FDS_Validation_Guide/SCRIPT_FIGURES/Scaling_Pyrolysis/%s_cone_%s,'%(chid, ('%0.1f'%(thickness)).replace('.','p'))
                outTxt = outTxt + 'Scaling Heat Release Rate Per Unit Area,max,0,' + materialClass + "," + materialMarker+materialColor + "," + materialColor +",TeX\n"
                
                outTxt3 = outTxt3 + "elif material == '%s':\n"%(material)
                outTxt3 = outTxt3 + "    thickness = %0.8f\n"%(thicknesses[iii])
                outTxt3 = outTxt3 + "    tMax_by_thickness[thickness_rounded] = %0.1f\n"%(tMax_by_thickness[thickness_rounded])
                outTxt3 = outTxt3 + "    qMax_by_thickness[thickness_rounded] = %0.1f\n"%(qMax_by_thickness[thickness_rounded])

        
    with open('scaling_dataplot_out.csv', 'w') as f:
        f.write(outTxt)
    
    
    with open(os.path.join(systemPath,'..','..','Run_All.sh'),'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.\n')
        f.write('# See the file Validation/Common_Run_All.sh for more information.\n')
        f.write('export SVNROOT=`pwd`/../..\n')
        f.write('source $SVNROOT/Validation/Common_Run_All.sh\n\n')
        
        files = sorted(glob.glob(os.path.join(systemPath,'..','*.fds')))
        for file in files:
            f.write('$QFDS $DEBUG $QUEUE -d $INDIR %s\n'%(os.path.basename(file)))
        f.write('\necho FDS cases submitted')
    
    
    materialClasses = [spec_file_dict[material]['materialClass'] for material in completeMaterials]
    series = [spec_file_dict[material]['series'] for material in completeMaterials]
    
    outTxt = '\\begin{table}[!h]\n'
    outTxt = outTxt + '\\caption[Summary of Materials]{Summary of Materials}\n'
    outTxt = outTxt + '\\centering\n'
    outTxt = outTxt + '\\begin{tabular}{|l|p{1.6cm}|p{1.6cm}|p{1.6cm}|p{1.6cm}|p{1.6cm}|}\n'
    outTxt = outTxt + '\\hline\n'
    outTxt = outTxt + 'Dataset'.ljust(20) + '& '
    for mc in sorted(list(set(materialClasses))):
        outTxt = outTxt + mc.ljust(12) + '& '
    outTxt = outTxt[:-2] + '& Total'.ljust(14)+'\\\\ \\hline\n'
    material_count = dict()
    for s in sorted(list(set(series))):
        material_count[s] = dict()
        outTxt = outTxt + s.replace('_',' ').ljust(20) + '& '
        counter = 0
        for mc in sorted(list(set(materialClasses))):
            material_count[s][mc] = len([x for x, y in zip(materialClasses, series) if ((mc == x) and (s == y))])
            outTxt = outTxt + ('%d'%(material_count[s][mc])).ljust(12) + '& '
            counter += material_count[s][mc]
        outTxt = outTxt + ('%d'%(counter)).ljust(12) + '\\\\ \\hline\n'
    outTxt = outTxt + 'Total'.ljust(20) + '& '
    counter = 0
    for mc in sorted(list(set(materialClasses))):
        count = len([x for x in materialClasses if (mc == x)])
        outTxt = outTxt + ('%d'%(count)).ljust(12) + '& '
        counter = counter + count
    outTxt = outTxt[:-2] +  ('& %d'%(counter)).ljust(14) +'\\\\ \\hline\n'
    outTxt = outTxt + '\\end{tabular}\n'
    outTxt = outTxt + '\\label{Scaling_Pyrolysis_Materials}\n'
    outTxt = outTxt + '\\end{table}\n'
    
    with open('material_summary.tex', 'w') as f:
        f.write(outTxt)
    
    
    
    
    citation_dict = {'Aalto Woods': 'Rinta-Paavola:2023',
                     'FAA Polymers': 'Stoliarov:CF2009,Stoliarov:FM2012',
                     'FPL Materials': 'FPL:Fire_Database',
                     'FSRI Materials': 'McKinnon:FSRI2023_Data',
                     'JH Materials': 'Luo:FRA2019,Lattimer:NIJ19',
                     'RISE Materials': 'RISE:Fire_Database'}
    class_dict = {'Mixtures': 'mixture', 'Others': 'other', 'Polymers': 'polymer', 'Wood-Based': 'Wood-Based'}
    
    outTxt = ''
    series = [spec_file_dict[material]['series'] for material in completeMaterials]
    for s in sorted(list(set(series))):
        s_name = s.replace('_',' ')
        s_cite = citation_dict[s_name]
        for mc in sorted(list(set(materialClasses))):
            if material_count[s][mc] == 0: continue
            mc_name = class_dict[mc]
            outTxt = outTxt + '\\begin{table}[!h]\n'
            outTxt = outTxt + '\\caption[Properties of %s, %s materials]{Properties of %s, %s materials ~\\cite{%s}.}\n'%(s_name, mc_name, s_name, mc_name, s_cite)
            outTxt = outTxt + '\\centering\n'
            outTxt = outTxt + '\\begin{tabular}{|p{5.5cm}|p{1.0cm}|p{1.0cm}|p{0.8cm}|p{1.4cm}|p{1.0cm}|p{1.0cm}|c|}\n'
            outTxt = outTxt + '\\hline\n'
            outTxt = outTxt + ' '.ljust(12) + '& \\centering$k$'.ljust(12) + '& \\centering$\\rho$'.ljust(12)+'& \\centering$\\varepsilon$'.ljust(12)+'& \\centering$c_{p}$'.ljust(12)+'& \\centering$T_{\\mathrm{ign}}$'.ljust(12) + '&\\centering$\\Delta H_{c}$' + '& $Y_{s}$  \\\\\n'
            outTxt = outTxt + 'Material'.ljust(12) + '& $\\mathrm{\\left(\\frac{W}{m\\cdot K}\\right)}$ & $\\mathrm{\\left(\\frac{kg}{m^{3}}\\right)}$ & $\\mathrm{( - )}$ & $\\mathrm{\\left(\\frac{kJ}{(kg\\cdot ^{\\circ}C)}\\right)}$ &  ($\\mathrm{^{\\circ}C}$)   & $\\left(\\mathrm{\\frac{kJ}{kg}}\\right)$ & $\\mathrm{\\left(-\\right)}$ \\\\ \\hline\n'
            outTxt = outTxt + '\\hline\n'
            for material in sorted(list(completeMaterials)):
                if (spec_file_dict[material]['materialClass'] != mc) or (spec_file_dict[material]['series'] != s): continue
                k = spec_file_dict[material]['conductivity']
                rho = spec_file_dict[material]['density']
                eps = spec_file_dict[material]['emissivity']
                cp = spec_file_dict[material]['specific_heat']
                DHc = spec_file_dict[material]['heat_of_combustion']*1e3
                Ys = spec_file_dict[material]['soot_yield']
                Tign = spec_file_dict[material]['Tign']
                
                outTxt = outTxt + material.ljust(50).replace('_',' ') + '& %0.2f & %0.0f & %0.2f & %0.2f & %0.1f & %0.0f & %0.4f \\\\\\hline\n'%(k, rho, eps, cp, Tign, DHc, Ys)
            outTxt = outTxt + '\\end{tabular}\n'
            outTxt = outTxt + '\\label{Properties_%s_%s}\n'%(s,mc)
            outTxt = outTxt + '\\end{table}\n\n\n'
    
    with open('material_properties.tex', 'w') as f:
        f.write(outTxt)
    
    sending_dict = {'Aalto_Woods': '',
                     'FAA_Polymers': '',
                     'FPL_Materials': '',
                     'FSRI_Materials': 'mc',
                     'JH_Materials': '',
                     'RISE_Materials': 'mc'}
    outTxt = ''
    series = [spec_file_dict[material]['series'] for material in completeMaterials]
    for s in sorted(list(set(series))):
        s_name = s.replace('_',' ')
        s_cite = citation_dict[s_name]
        for mc in sorted(list(set(materialClasses))):
            if material_count[s][mc] == 0: continue
            mc_name = class_dict[mc]
            sending = sending_dict[s]
            if sending == 'mc': sending = ', %s materials'%(mc)
            outTxt = outTxt + '\\begin{figure}[p]\n'
            outTxt = outTxt + '\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}r}\n'
            counter = 0
            for material in sorted(list(completeMaterials)):
                if (spec_file_dict[material]['materialClass'] != mc) or (spec_file_dict[material]['series'] != s): continue
                chid = material.replace(' ','_')
                cases = get_filtered_cases(spec_file_dict, material, energyThreshold=energyThreshold)
                thicknesses = [cases[c]['delta']*1e3 for c in cases.keys()]
                for thickness in sorted(list(set(thicknesses))):
                    fname = ('%s_cone_%s'''%(chid, ('%0.1f'%(thickness)))).replace('.','p')
                    if counter < 8:
                        ending = '&' if (counter%2 == 0) else '\\\\'
                        outTxt = outTxt + '\\includegraphics[height=2.10in]{SCRIPT_FIGURES/Scaling_Pyrolysis/%s.pdf} %s\n'%(fname, ending)
                        counter += 1
                    else:
                        outTxt = outTxt + '\\end{tabular*}\n'
                        outTxt = outTxt + '\\caption[HRRPUA of %s using scaling model %s]\n'%(s.replace('_',' '), sending)
                        outTxt = outTxt + '{Comparison of predicted and measured heat release rate per unit area using scaling-based approach for cone calorimeter experiments.}\n'
                        outTxt = outTxt + '\\label{%s_HRR_%s}\n'%(s, mc)
                        outTxt = outTxt + '\\end{figure}\n\n'
                        outTxt = outTxt + '\\begin{figure}[p]\n'
                        outTxt = outTxt + '\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}r}\n'
                        ending = '&'
                        outTxt = outTxt + '\\includegraphics[height=2.10in]{SCRIPT_FIGURES/Scaling_Pyrolysis/%s.pdf} %s\n'%(fname, ending)
                        counter = 1
            if counter > 0:
                outTxt = outTxt + '\\end{tabular*}\n'
                outTxt = outTxt + '\\caption[HRRPUA of %s using scaling model %s]\n'%(s.replace('_',' '), sending)
                outTxt = outTxt + '{Comparison of predicted and measured heat release rate per unit area using scaling-based approach for cone calorimeter experiments.}\n'
                outTxt = outTxt + '\\label{%s_HRR_%s}\n'%(s, mc)
                outTxt = outTxt + '\\end{figure}\n\n'
                outTxt = outTxt + '\\begin{figure}[p]\n'
                outTxt = outTxt + '\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}r}\n'
    
    with open('material_figures.tex', 'w') as f:
        f.write(outTxt)
    
    
    
    
    #from materials.scripts.algorithms import get_ignition_temperature
    #cases = spec_file_dict[material]['cases']
    #baseDir = os.path.join(systemPath,'materials','input_files')
    #get_ignition_temperature(material, baseDir, cases)
    
    material_count = pd.DataFrame(material_count)
    material_output_data = pd.DataFrame(material_output_data)
    material_output_data.to_csv('material_output_data.csv', float_format='%.1f')
    
    
    
    
    
    
    
    
    files = glob.glob(os.path.join(systemPath, 'materials','input_files','*','*_git.txt'))
    for file in files:
        ignores = ['-29mm-15.csv','-29mm-20.csv','-29mm-25.csv','-29mm-35.csv','-29mm-50.csv','-29mm-75.csv']
        if os.path.basename(file) not in ignores:
            shutil.copy(file, os.path.join(systemPath,'..','..','..','..','..','out','Scaling_Pyrolysis',os.path.basename(file)))
    
    files = glob.glob(os.path.join(systemPath, 'materials','input_files','*','*_devc.csv'))
    for file in files:
        ignores = ['-29mm-15.csv','-29mm-20.csv','-29mm-25.csv','-29mm-35.csv','-29mm-50.csv','-29mm-75.csv']
        if os.path.basename(file) not in ignores:
            shutil.copy(file, os.path.join(systemPath,'..','..','..','..','..','out','Scaling_Pyrolysis',os.path.basename(file)))
    
    
    

    
