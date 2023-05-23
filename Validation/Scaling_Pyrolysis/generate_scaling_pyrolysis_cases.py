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
import os
import subprocess
import pandas as pd
import shutil
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
    
    process = subprocess.Popen([fdscmd, outfile, ">&", "log.err"], cwd=r'%s'%(outdir), env=my_env, shell=False, stdout=subprocess.DEVNULL)
    
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

def interpolateExperimentalData(times, HRRs, targetDt=False, filterWidth=False):
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
    HRRs = np.interp(targetTimes, times, HRRs)
    
    return targetTimes, HRRs

def getRepresentativeHrrpua(HRRPUA, factor=0.5):
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
    representativeHRRPUA = HRRPUA[HRRPUA > HRRPUA.max()*factor].mean()
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
    hrrpua_ref = getRepresentativeHrrpua(hrrpua)
    qref = estimateExposureFlux(coneExposure, hrrpua_ref)
    
    tempOutput = '.TRUE.' if outputTemperature else '.FALSE.'
    DT_DEVC = devc_dt
    if calculateDevcDt:
        NFRAMES = 1200/1.
        DT_DEVC = tend/NFRAMES
        
    txt = "&HEAD CHID='%s', /\n"%(chid)
    txt = txt+"&TIME DT=1., T_END=%0.1f /\n"%(tend)
    txt = txt+"&DUMP DT_CTRL=%0.1f, DT_DEVC=%0.1f, DT_HRR=%0.1f, SIG_FIGS=4, SIG_FIGS_EXP=2, /\n"%(DT_DEVC, DT_DEVC, DT_DEVC)
    txt = txt+"&MISC SOLID_PHASE_ONLY=.TRUE., TMPA=27., /\n"
    txt = txt+"&MESH ID='MESH', IJK=3,3,3, XB=0.,0.3,0.,0.3,0.,0.3, /\n"
    txt = txt+"&REAC ID='PROPANE', FUEL='PROPANE', /\n"
    txt = txt+"&MATL ID='BACKING', CONDUCTIVITY=0.10, DENSITY=65., EMISSIVITY=0.9, SPECIFIC_HEAT=1.14, /\n"
    #txt = txt+"&MATL ID='BACKING', CONDUCTIVITY=0.2, DENSITY=585., EMISSIVITY=1., SPECIFIC_HEAT=0.8, /\n"
    txt = txt+"&MATL ID='SAMPLE', CONDUCTIVITY=%0.4f, DENSITY=%0.1f, EMISSIVITY=%0.4f, SPECIFIC_HEAT=%0.4f, /\n"%(k, rho, e, cp)
    
    for i in range(0, len(time)):
        txt = txt+"&RAMP ID='CONE-RAMP', T=%0.4f, F=%0.4f, /\n"%(time[i]-time[0], hrrpua[i])
    
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
    if np.var(np.log(y / x)) - (sigma_e**2) < 0:
        sigma_m = sigma_e
    else:
        sigma_m = (np.var(np.log(y / x)) - (sigma_e**2))**0.5
    #sigma_m2 = np.var(np.log(y / x)) / 2
    sigma_m = np.max([sigma_m, sigma_e])
    delta = np.exp(np.mean(np.log(y / x)) + (sigma_m**2)/2 - (sigma_e**2)/2)
    return delta, sigma_m, sigma_e, np.log(y/x)

def calculateUncertaintyBounds(flatx, flaty, flatFlux, split=False):
    d2 = pd.DataFrame(np.array([flatx, flaty, flatFlux]).T, columns=['exp','mod','flux'])
    d2[d2 == 0] = np.nan
    #d2[d2 < 0] = np.nan
    mask = np.logical_and(~np.isnan(np.array(d2.values[:,0], dtype=float)),
                          ~np.isnan(np.array(d2.values[:,1], dtype=float)))
    d = d2.loc[mask]
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
        num_points = x.shape[0]
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
    
    
    
    if diff is not None:
        cases = np.array(list(set(diff)))
        for j in range(0, len(f)):
            caseInd = np.where(cases == diff[j])[0][0]
            #c = 0 if diff[j] > 0 else 1
            ax.scatter(x[j], y[j], marker=markers[caseInd], s=100, color=colors2[caseInd])
        customMarkers = []
        for caseInd, case in enumerate(cases):
            if labelName is None:
                customMarkers.append(Line2D([0],[0],marker=markers[caseInd], color='w', markerfacecolor=colors2[caseInd], label=case, markersize=15))
            else:
                customMarkers.append(Line2D([0],[0],marker=markers[caseInd], color='w', markerfacecolor=colors2[caseInd], label=labelName[case], markersize=15))
        
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

def computeAnova(delta, sigma_m, num_points):
    means = dict()
    stds = dict()
    nums = dict()
    
    for key in list(sigma_m.keys()):
        means[key] = 1/delta[key]
        stds[key] = (sigma_m[key]**2)*(means[key]**2)
        nums[key] = num_points[key]
    
    overallMean = np.mean([float(x) for x in means.values()])
    
    # Calculate the between groups SSD across full set
    BSSD = 0
    for key in list(sigma_m.keys()):
        BSSD += nums[key]*(means[key] - overallMean)**2
    
    # Between groups degrees of freedom is one less than number of groups
    BDOF = len(means.keys())-1
    
    # Between goups mean square value is SSD divided by DOF
    BMSV = SSD/BDOF
    
    # Within groups SSD
    ISSD = 0
    for key in list(sigma_m.keys()):
        ISSD += nums[key]*(stds[key]**2)
    
    # In groups DOF is number of groups x (number of samples in group - 1)
    IDOF = dict()
    for key in list(sigma_m.keys()):
        IDOF[key] = len(means.keys())*(nums[key]-1)
    
    # In groups MSV is SSD / DOF
    IMSV = dict()
    for key in list(sigma_m.keys()):
        IMSV[key] = ISSD/IDOF[key]
        
    # F ratio
    Fratio = dict()
    for key in list(sigma_m.keys()):
        Fratio[key] = BMSV/IMSV[key]
    
    # F threshold
    Fthresh = dict()
    for key in list(sigma_m.keys()):
        Fthresh[key] = scipy.stats.f.ppf(0.05, IDOF[key], BDOF)
    
    return Fratio, Fthresh

def getNormalStats(delta, sigma_m):
    mu = 1/delta
    sig = ((sigma_m**2)*(mu**2))**0.5
    var = sig ** 2
    
    return mu, var

if __name__ == "__main__":
    
    fdsdir, fdscmd = findFds()
    
    specificationFile = "scaling_pyrolysis_cone_cases.csv"
    figoutdir = "figures"
    
    if figoutdir is not None:
        if os.path.exists(figoutdir) is not True: os.mkdir(figoutdir)
        import matplotlib.pyplot as plt
    
    specificationFile = pd.read_csv(specificationFile)
    material_uncertainty = defaultdict(bool)
    material_output_data = defaultdict(bool)
    uncertainty = defaultdict(bool)
    uncertainty[60] = defaultdict(bool)
    uncertainty[180] = defaultdict(bool)
    uncertainty[300] = defaultdict(bool)
    uncertainty['peak'] = defaultdict(bool)
    uncertainty['t_peak'] = defaultdict(bool)
    
    for key in list(uncertainty.keys()):
        uncertainty[key]['EXP'] = []
        uncertainty[key]['MOD'] = []
    uncertainty['flux'] = []
    uncertainty['delta'] = []
    uncertainty['case'] = []
    uncertainty['material'] = []
    uncertainty['series'] = []
    uncertainty['MaterialClass'] = []
    
    runSimulations = True
    showStats = False
    closePlots = False
    
    for i in range(1, specificationFile.shape[0]):
        
        # Check for run code
        code = specificationFile.iloc[i]['Code']
        num_id = specificationFile.iloc[i]['Number']
        if code == 's':
            print("Skipping row %d"%(i))
            continue
        elif code =='d':
            pass
        else:
            print("Unknown code %s in row %d"%(code, i))
            continue
        
        # Extract specification file data
        material = specificationFile.iloc[i]['Material']
        series = specificationFile.iloc[i]['Series']
        materialClass = specificationFile.iloc[i]['MaterialClass']
        coneExposure = specificationFile.iloc[i]['ReferenceExposure']
        conductivity = specificationFile.iloc[i]['Conductivity']
        specific_heat = specificationFile.iloc[i]['SpecificHeat']
        density = specificationFile.iloc[i]['Density']
        emissivity = specificationFile.iloc[i]['Emissivity']
        thickness = specificationFile.iloc[i]['Thickness']
        
        resultDir = specificationFile.iloc[i]['ResultDir'].replace('\\\\','\\').replace('"','')
        if os.path.exists(resultDir) is not True: os.mkdir(resultDir)
        workingDir = os.path.join(resultDir, 'tmp')
        if os.path.exists(workingDir) is not True: os.mkdir(workingDir)
        inputFileDir = specificationFile.iloc[i]['InputFileDir'].replace('\\\\','\\').replace('"','')
        expFileDir = specificationFile.iloc[i]['ExpFileDir'].replace('\\\\','\\').replace('"','')
        
        referenceTimeColumn = specificationFile.iloc[i]['ReferenceTime']
        referenceHrrpuaColumn = str(specificationFile.iloc[i]['ReferenceHRRPUA'])
        
        validationTimeColumns = specificationFile.iloc[i]['ValidationTimes'].split('|')
        validationHrrpuaColumns = specificationFile.iloc[i]['ValidationHrrpuaColumns'].split('|')
        validationFluxes = specificationFile.iloc[i]['ValidationFluxes'].split('|')
        fluxes = [float(f) for f in validationFluxes]
        
        ignitionTemperature = specificationFile.iloc[i]['IgnitionTemperature']
        if ignitionTemperature == 'Calculate':
            calculateIgnitionTemperature = True
            ignitionTemperatureBasis = specificationFile.iloc[i]['IgnitionTemperatureBasis'].split('|')
            ignitionTemperatureBasis = [float(x) for x in ignitionTemperatureBasis]
            Tign = 1000
        else:
            Tign = float(ignitionTemperature)
            calculateIgnitionTemperature = False
        dataFile = specificationFile.iloc[i]['DataFile'].replace('\\\\','\\').replace('"','')
        headerRows = specificationFile.iloc[i]['HeaderRows']
        if '|' in dataFile:
            dfs = dataFile.split('|')
            hrs = headerRows.split('|')
            exp_data = dict()
            for df, hr in zip(dfs, hrs):
                fname = df.split(os.sep)[-1]
                # Read data file, manually due to differing number of header rows
                with open(df, 'r') as f:
                    d = f.readlines()
                d = np.array([dd.replace('\n','').replace('/','_').split(',') for dd in d])
                hr = int(hr)
                for ii in range(hr, len(d)):
                    for j in range(0, len(d[ii])):
                        try:
                            d[ii,j] = float(d[ii,j])
                        except:
                            d[ii,j] = np.nan
                columns = [fname + '-' + str(c) for c in d[0]]
                for ii, c in enumerate(columns):
                    c2 = os.path.abspath(c).split(os.sep)[-1]
                    exp_data[c2] = np.array(d[hr:, ii], dtype=float)
            multipleFiles = True
        else:
            headerRows = int(headerRows)
            # Read data file, manually due to differing number of header rows
            with open(dataFile, 'r') as f:
                d = f.readlines()
            d = np.array([dd.replace('\n','').split(',') for dd in d])
            
            for ii in range(headerRows, len(d)):
                for j in range(0, len(d[ii])):
                    try:
                        d[ii,j] = float(d[ii,j])
                    except:
                        d[ii,j] = np.nan
            columns = [str(c) for c in d[0]]
            exp_data = pd.DataFrame(np.array(d[headerRows:, :], dtype=float), columns=columns)
            multipleFiles = False
        if multipleFiles:
            # Identify time to ignition
            HRRs = exp_data[referenceHrrpuaColumn]
            times = exp_data[referenceTimeColumn]
        else:
            # Identify time to ignition
            HRRs = exp_data.loc[~np.isnan(exp_data[referenceTimeColumn]),referenceHrrpuaColumn].values
            times = exp_data.loc[~np.isnan(exp_data[referenceTimeColumn]),referenceTimeColumn].values
        #targetTimes, HRRs_interp = interpolateExperimentalData(times, HRRs, targetDt=15, filterWidth=False)
        tign, times_trimmed, hrrs_trimmed = findLimits(times, HRRs)
            
        tigns = dict()
        for ii in range(0, len(validationTimeColumns)):
            timeColumn = validationTimeColumns[ii]
            hrrpuaColumn = validationHrrpuaColumns[ii]
            if multipleFiles:
                HRRs_v = exp_data[hrrpuaColumn]
                times_v = exp_data[timeColumn]
            else:
                HRRs_v = exp_data.loc[~np.isnan(exp_data[timeColumn]),hrrpuaColumn].values
                times_v = exp_data.loc[~np.isnan(exp_data[timeColumn]),timeColumn].values
            tign_v = findIgnitionTime(times_v, HRRs_v)
            tigns[fluxes[ii]] = tign_v
        
        # Calculate reference heat flux
        # If the trimmed HRR curve fails, use the full curve and print a warning
        try:
            hrrpua_ref = getRepresentativeHrrpua(hrrs_trimmed)
        except:
            hrrpua_ref = getRepresentativeHrrpua(HRRs)
            print("Warning: Failed to get representative HRRPUA for material %s on trimmed data. Using full HRR curve."%(material))
            hrrs_trimmed = HRRs
            times_trimmed = times
        qref = estimateExposureFlux(coneExposure, hrrpua_ref)
        
        # Set chid
        chid = series + '_' + material.replace(' ','_')+"_cone"
        if len(chid) > 50:
            chid = chid[:50]
        
        if calculateIgnitionTemperature and runSimulations:
            # Generate initial FDS input file to calculate ignition temperature
            tend = times.max() + 300 # Arbitrary, overwritten in actual calculation after finding Tign
            txt = buildFdsFile(chid, coneExposure, emissivity, conductivity, density, 
                                   specific_heat, 1000, thickness, times_trimmed, hrrs_trimmed,
                                   tend, fluxes, 15.0, ignitionMode='Temperature', outputTemperature=True,
                                   calculateDevcDt=False, devc_dt=0.1)
            '''
            (chid, coneExposure, e, k, rho, cp, Tign, d, time, hrrpua, tend,
                             HFs, front_h, hfs_tign=False,
                             ignitionMode='Temperature', outputTemperature=False,
                             calculateDevcDt=True, devc_dt=1.)
            '''
            with open("%s%s%s.fds"%(workingDir, os.sep, chid), 'w') as f:
                f.write(txt)
            
            if runSimulations:
                runModel(workingDir, chid+".fds", 1, fdsdir, fdscmd, printLiveOutput=False)
            data = load_csv(workingDir, chid)
            
            Tigns = []
            Tmaxs = []
            for flux in ignitionTemperatureBasis:
                tign = tigns[flux]
                Tign = data.loc[data['Time'] > tign, '"WALL TEMPERATURE-%02d"'%(flux)].values[0]
                Tmaxs.append(np.nanmax(data['"WALL TEMPERATURE-%02d"'%(flux)].values))
                Tigns.append(Tign)
            Tign = np.mean(Tigns)
            t = Tigns.copy()
            while Tign > np.min(Tmaxs):
                t = [x for x in t if x < np.max(t)]
                Tign = np.max(t)
            print(i, "Material %s ignition temperature %0.1f reference flux %0.1f"%(material, Tign, qref))
        else:
            tend = times.max()*1.5 + 300
            print(i, "Material %s ignition temperature %0.1f reference flux %0.1f"%(material, Tign, qref))
        
        try:
            tMaxes = [exp_data[key].max() for key in list(exp_data.keys()) if 'Time' in key]
            tend = max(tMaxes) * 1.5
        except:
            tend = times.max()*1.5
        # Generate fds input file with updated ignition temperature
        txt = buildFdsFile(chid, coneExposure, emissivity, conductivity, density, 
                               specific_heat, Tign, thickness, times_trimmed, hrrs_trimmed,
                               tend, fluxes, 0.0, hfs_tign=tigns, ignitionMode='Time',
                               outputTemperature=False, calculateDevcDt=True)
        
        with open("%s%s%s.fds"%(workingDir, os.sep, chid), 'w') as f:
            f.write(txt)
        if runSimulations:
            runModel(workingDir, chid+".fds", 1, fdsdir, fdscmd, printLiveOutput=False)
        data = load_csv(workingDir, chid)
        
        # Plot results
        if figoutdir is not None:
            fig = plotResults_exp(data, exp_data, fluxes, validationTimeColumns, validationHrrpuaColumns, lw=3, fs=16)
            fig.savefig(os.path.join(figoutdir, '%03d_'%(i)+chid+'.png'), dpi=300)
            if closePlots:
                plt.close()
        
        
        validationFluxes = list(fluxes)
        validationFluxes.remove(coneExposure)
        
        validationHrrpuaColumns.remove(referenceHrrpuaColumn)
        validationTimeColumns.remove(referenceTimeColumn)
        
        for flux, hc, tc in zip(validationFluxes, validationHrrpuaColumns, validationTimeColumns):
            time = data['Time'].values
            model_dt = np.median(time[1:]-time[:-1])
            
            scaling = data['"HRRPUA-%02d"'%(flux)].values
            if type(exp_data) is dict:
                expTime = exp_data[tc]
                expHRR = exp_data[hc]
            else:
                expTime = exp_data[tc].values
                expHRR = exp_data[hc].values
            
            uncertainty['peak']['EXP'].append(np.nanmax(expHRR))
            uncertainty['peak']['MOD'].append(np.nanmax(scaling))
            
            uncertainty['t_peak']['EXP'].append(expTime[np.argmax(expHRR)])
            uncertainty['t_peak']['MOD'].append(time[np.argmax(scaling)])
            uncertainty['flux'].append(flux)
            uncertainty['delta'].append(flux-coneExposure)
            uncertainty['case'].append(num_id)
            uncertainty['material'].append(material)
            uncertainty['series'].append(series)
            uncertainty['MaterialClass'].append(materialClass)
            for filterWidth in [60, 180, 300]:
                filterWidth_mod = int(filterWidth/model_dt)
                fil_mod = np.ones(filterWidth_mod)/filterWidth_mod
                scaling_filtered = np.convolve(scaling, fil_mod, mode='same')
                uncertainty[filterWidth]['MOD'].append(np.nanmax(scaling_filtered))
                
                exp_dt = np.nanmin(expTime[1:]-expTime[:-1])
                exp_times2 = np.linspace(0, np.nanmax(expTime), int(np.round(np.nanmax(expTime)/exp_dt)+1))
                exp_hrr2 = np.interp(exp_times2, expTime, expHRR)
                
                filterWidth_exp = int(filterWidth/exp_dt)
                fil_exp = np.ones(filterWidth_exp)/filterWidth_exp
                exp_filtered = np.convolve(exp_hrr2, fil_exp, mode='same')
                uncertainty[filterWidth]['EXP'].append(np.nanmax(exp_filtered))
        
        # Copy results to directories for building guide
        shutil.copy(os.path.join(workingDir, chid+"_devc.csv"), resultDir)
        shutil.copy(os.path.join(workingDir, chid+"_git.txt"), resultDir)
        shutil.copy(os.path.join(workingDir, chid+".fds"), inputFileDir)
        
        #if material_output_data[series] is False:
        #    material_output_data[series] = defaultdict(bool)
        material_output_data[(series,material)] = defaultdict(bool)
        material_output_data[(series,material)]['Conductivity\n($\mathrm{W/(m\cdot K)}$)'] = conductivity
        material_output_data[(series,material)]['Density\n($\mathrm{kg/m^{3}}$)'] = density
        material_output_data[(series,material)]['Emissivity\n(-)'] = emissivity
        material_output_data[(series,material)]['Specific Heat\n($\mathrm{kJ/(kg\cdot C}$)'] = specific_heat
        material_output_data[(series,material)]['Thickness\n($\mathrm{mm}$)'] = thickness
        material_output_data[(series,material)]['Ignition Temperaure\n($\mathrm{^{\circ}C}$)'] = Tign
        material_output_data[(series,material)]['Reference Heat Flux\n($\mathrm{kW/m^{2}}$)'] = qref
        
    material_output_data = pd.DataFrame(material_output_data)
    material_output_data.to_csv('material_output_data.csv')
    
    material_data = material_output_data.T
    
    
    if showStats:
        axmin, axmax = (1e1, 1e4)
        
        axmin = 0
        axmax = 3000
        
        cases = np.array(uncertainty['case'])
        uncertainties = [60, 180, 300, 't_peak', 'peak']
        labels = ['60s Avg', '180s Avg', '300s Avg', 'Time to Peak', 'Peak']
        diff2 = np.array([-1 if d < 0 else 1 for d in uncertainty['delta']])
        labelNames = {1 : '$\mathrm{q_{exp} > q_{ref}}$', -1 : '$\mathrm{q_{exp} < q_{ref}}$' }
        for loglog in [False, True]:
            if loglog:
                axmin = 5e0
                axmaxes = [5e3, 5e3, 5e3, 1e4, 5e3]
            else:
                axmin = 0
                axmaxes = [2500, 2000, 1500, 10000, 3000]
            for i in range(0, len(uncertainties)):
                label = labels[i]
                axmax = axmaxes[i]
                u = uncertainties[i]
                x = np.array(uncertainty[u]['EXP'])
                y = np.array(uncertainty[u]['MOD'])
                f = np.array(uncertainty['flux'])
                fig, sigma_m, delta = plotMaterialExtraction(x, y, f, label, diff=diff2, axmin=axmin, axmax=axmax, loglog=loglog, labelName=labelNames)
                if loglog:
                    fname = os.path.join(figoutdir, 'statistics_loglog_direction_%s.png'%(label))
                else:
                    fname = os.path.join(figoutdir, 'statistics_nolog_direction_%s.png'%(label))
                plt.savefig(fname, dpi=300)
                plt.close()
                
                
                
                print(label, "\tBias\tSigma_m")
                delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x, y, diff2, split=False)
                print("Total\t\t%0.4f\t%0.4f"%(delta, sigma_m))
                delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x, y, diff2, split=True)
                print("Lower\t\t%0.4f\t%0.4f"%(delta[-1], sigma_m[-1]))
                print("Upper\t\t%0.4f\t%0.4f"%(delta[1], sigma_m[1]))
        
        
        
        cases = np.array(uncertainty['case'])
        uncertainties = [60, 't_peak', 'peak']
        labels = ['60s Avg', 'Time to Peak', 'Peak']
        
        for loglog in [False, True]:
            if loglog:
                axmin = 5e0
                axmaxes = [5e3, 1e4, 5e3]
            else:
                axmin = 0
                axmaxes = [2500, 10000, 3000]
            for i in range(0, len(uncertainties)):
                label = labels[i]
                axmax = axmaxes[i]
                u = uncertainties[i]
                x = np.array(uncertainty[u]['EXP'])
                y = np.array(uncertainty[u]['MOD'])
                f = np.array(uncertainty['flux'])
                fig, sigma_m, delta = plotMaterialExtraction(x, y, f, label, diff=uncertainty['MaterialClass'], axmin=axmin, axmax=axmax, loglog=loglog)
                if loglog:
                    fname = os.path.join(figoutdir, 'statistics_loglog_matclass_%s.png'%(label))
                else:
                    fname = os.path.join(figoutdir, 'statistics_nolog_matclass_%s.png'%(label))
                plt.savefig(fname, dpi=300)
                plt.close()
                
        cases = np.array(uncertainty['case'])
        uncertainties = [60, 't_peak', 'peak']
        labels = ['60s Avg', 'Time to Peak', 'Peak']
        
        for loglog in [False, True]:
            if loglog:
                axmin = 5e0
                axmaxes = [5e3, 1e4, 5e3]
            else:
                axmin = 0
                axmaxes = [2500, 10000, 3000]
            for i in range(0, len(uncertainties)):
                label = labels[i]
                axmax = axmaxes[i]
                u = uncertainties[i]
                x = np.array(uncertainty[u]['EXP'])
                y = np.array(uncertainty[u]['MOD'])
                f = np.array(uncertainty['flux'])
                fig, sigma_m, delta = plotMaterialExtraction(x, y, f, label, diff=uncertainty['series2'], axmin=axmin, axmax=axmax, loglog=loglog)
                if loglog:
                    fname = os.path.join(figoutdir, 'statistics_loglog_series_%s.png'%(label))
                else:
                    fname = os.path.join(figoutdir, 'statistics_nolog_series_%s.png'%(label))
                plt.savefig(fname, dpi=300)
                plt.close()
        
        
        
        i = 3
        label = labels[i]
        axmax = axmaxes[i]
        u = uncertainties[i]
        x = np.array(uncertainty[u]['EXP'])
        y = np.array(uncertainty[u]['MOD'])
        f = np.array(uncertainty['flux'])
        
        delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x, y, diff2, split=False)
        outputs = dict()
        outputs['delta'] = dict()
        outputs['sigma_m'] = dict()
        outputs['count'] = dict()
        #outputs['delta']['all'] = delta
        #outputs['sigma_m']['all'] = sigma_m
        #outputs['count']['all'] = len(x)
        hrrpua_thresholds = np.logspace(1,4, num=100)
        for j in hrrpua_thresholds: #[50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]:
            x2 = x[x > j]
            y2 = y[x > j]
            diff3 = diff2[x > j]
            delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x2, y2, diff3, split=False)
            outputs['delta'][j] = delta
            outputs['sigma_m'][j] = sigma_m
            outputs['count'][j] = len(x2)
        outputs2 = pd.DataFrame(outputs)
        
        ind = np.where(outputs2['sigma_m'] == outputs2['sigma_m'].min())[0][0]
        hrrpua_at_min = hrrpua_thresholds[ind]
        sigma_m_at_max = outputs2['sigma_m'].max()
        sigma_m_at_min = outputs2['sigma_m'].min()
        slope = -1 * (sigma_m_at_max - sigma_m_at_min) / (hrrpua_at_min - 0) *1.05
        
        approx = hrrpua_thresholds*slope + sigma_m_at_max*1.05
        approx[approx < sigma_m_at_min] = sigma_m_at_min
        
        sigma_e = np.zeros_like(hrrpua_thresholds) + sigma_m_at_min
        
        fs=16
        colors = getColors()
        plt.figure(figsize=(10, 6))
        #plt.scatter(hrrpua_thresholds, outputs2['sigma_m'], label='$\mathrm{\sigma_{m}}$', s=30, color=colors[0])
        plt.plot(hrrpua_thresholds, outputs2['sigma_m'], label='$\mathrm{\sigma_{m}}$', color=colors[0], linewidth=3)
        app_label = "$\mathrm{\sigma_{m}}=%0.2f - (%0.1fE-4) q''_{peak}$"%(sigma_m_at_max*1.05, -1e4*slope)
        plt.plot(hrrpua_thresholds, sigma_e, '--', label='$\mathrm{\sigma_{e}}$', color=colors[2], linewidth=3)
        plt.plot(hrrpua_thresholds, approx, label=app_label, linewidth=3, color=colors[1])
        plt.xlabel('Peak HRRPUA ($\mathrm{kW/m^{2}}$)', fontsize=fs)
        plt.ylabel('$\mathrm{\sigma}$', fontsize=fs)
        plt.ylim(0.07, 0.18)
        plt.xlim(0, 1000)
        plt.tick_params(labelsize=fs)
        plt.legend(fontsize=fs)
        plt.grid()
        plt.savefig(os.path.join(figoutdir, 'statistics_sigma_m_with_hrrpua.png'), dpi=300)
        
        
        uncertainty['MaterialClass2'] = [x if (('Others' not in x) and ('Composites' not in x)) else 'Others & Composites' for x in uncertainty['MaterialClass']]
        uncertainty['series2'] = [x if (('JH' not in x)) else 'JH' for x in uncertainty['series']]
        
        
        print("Deviation by scaling up/down")
        cases = np.array(uncertainty['case'])
        uncertainties = [60, 'peak']
        labels = ['60s Avg', 'Peak']
        diff2 = np.array(['Lower' if d < 0 else 'Higher' for d in uncertainty['delta']])
        
        for i in range(0, len(uncertainties)):
            label = labels[i]
            axmax = axmaxes[i]
            u = uncertainties[i]
            x = np.array(uncertainty[u]['EXP'])
            y = np.array(uncertainty[u]['MOD'])
            
            print(label.ljust(30), "\tN\t\tBias\tSigma_m")
            tdelta, tsigma_m, tsigma_e, tnum_points, tpoints = calculateUncertaintyBounds(x, y, diff2, split=False)
            print("%s\t%0d\t\t%0.4f\t\t%0.4f"%('Total'.ljust(30), tnum_points, tdelta, tsigma_m))
            delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x, y, diff2, split=True)
            for s in np.unique(diff2):
                mn1, var1 = getNormalStats(delta[s], sigma_m[s])
                N1 = num_points[s]
                mn2, var2 = getNormalStats(tdelta, tsigma_m)
                N2 = tnum_points
                
                stat, p = scipy.stats.levene(points[s], tpoints, center='mean')
                #p = scipy.stats.f.cdf(var2/var1, N2, N1) if var1 > var2 else scipy.stats.f.cdf(var1/var2, N1, N2)
                
                print("%s\t%04d\t%0.4f\t\t%0.4f\t\t%0.4f"%(s.ljust(30), num_points[s], delta[s], sigma_m[s], p))
                
            print("\n")
            
            
        
        print("Deviation by series")
        uncertainties = [60, 'peak']
        labels = ['60s Avg', 'Peak']
        #diff2 = np.array([-1 if d < 0 else 1 for d in uncertainty['delta']])
        diff2 = np.array(uncertainty['series2'])
        
        for i in range(0, len(uncertainties)):
            label = labels[i]
            axmax = axmaxes[i]
            u = uncertainties[i]
            x = np.array(uncertainty[u]['EXP'])
            y = np.array(uncertainty[u]['MOD'])
            
            print(label.ljust(30), "\tN\tBias\tSigma_m")
            tdelta, tsigma_m, tsigma_e, tnum_points, tpoints = calculateUncertaintyBounds(x, y, diff2, split=False)
            print("%s\t%0d\t\t%0.4f\t\t%0.4f"%('Total'.ljust(30), tnum_points, tdelta, tsigma_m))
            delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x, y, diff2, split=True)
            
            for s in np.unique(diff2):
                mn1, var1 = getNormalStats(delta[s], sigma_m[s])
                N1 = num_points[s]
                mn2, var2 = getNormalStats(tdelta, tsigma_m)
                N2 = tnum_points
                
                stat, p = scipy.stats.levene(points[s], tpoints, center='mean')
                #p = scipy.stats.f.cdf(var2/var1, N2, N1) if var1 > var2 else scipy.stats.f.cdf(var1/var2, N1, N2)
                
                print("%s\t%04d\t%0.4f\t\t%0.4f\t\t%0.4f"%(s.ljust(30), num_points[s], delta[s], sigma_m[s], p))
                
                
            print("\n")
        
        print("Deviation by material class")
        uncertainties = [60, 'peak']
        labels = ['60s Avg', 'Peak']
        #diff2 = np.array([-1 if d < 0 else 1 for d in uncertainty['delta']])
        diff2 = np.array(uncertainty['MaterialClass'])
        
        for i in range(0, len(uncertainties)):
            label = labels[i]
            axmax = axmaxes[i]
            u = uncertainties[i]
            x = np.array(uncertainty[u]['EXP'])
            y = np.array(uncertainty[u]['MOD'])
            
            print(label.ljust(30), "\tN\t\tBias\t\tSigma_m\tP-value")
            tdelta, tsigma_m, tsigma_e, tnum_points, tpoints = calculateUncertaintyBounds(x, y, diff2, split=False)
            print("%s\t%04d\t%0.4f\t\t%0.4f"%('Total'.ljust(30), tnum_points, tdelta, tsigma_m))
            delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x, y, diff2, split=True)
            
            for s in np.unique(diff2):
                mn1, var1 = getNormalStats(delta[s], sigma_m[s])
                N1 = num_points[s]
                mn2, var2 = getNormalStats(tdelta, tsigma_m)
                N2 = tnum_points
                
                stat, p = scipy.stats.levene(points[s], tpoints, center='mean')
                #p = scipy.stats.f.cdf(var2/var1, N2, N1) if var1 > var2 else scipy.stats.f.cdf(var1/var2, N1, N2)
                
                print("%s\t%04d\t%0.4f\t\t%0.4f\t\t%0.4f"%(s.ljust(30), num_points[s], delta[s], sigma_m[s], p))
                
            print("\n")
    
        
        print("Deviation by series and material class")
        uncertainties = [60, 'peak']
        labels = ['60s Avg', 'Peak']
        diff2 = [u1 + u2 for u1, u2 in zip(uncertainty['series2'], uncertainty['MaterialClass'])]
        
        diff3 = [u1 for u1 in uncertainty['MaterialClass']]
        sep = 20
        for i in range(0, len(uncertainties)):
            label = labels[i]
            axmax = axmaxes[i]
            u = uncertainties[i]
            x = np.array(uncertainty[u]['EXP'])
            y = np.array(uncertainty[u]['MOD'])
            
            
            tdelta, tsigma_m, tsigma_e, tnum_points, tpoints = calculateUncertaintyBounds(x, y, diff3, split=True)
            #print("%s\t%0d\t\t%0.4f\t\t%0.4f"%('Total'.ljust(30), tnum_points, tdelta, tsigma_m))
            delta, sigma_m, sigma_e, num_points, points = calculateUncertaintyBounds(x, y, diff2, split=True)
            
            for j in np.unique(diff3):
                print(label.ljust(sep), "\tN\t\tBias\t\tSigma_m")
                mn2, var2 = getNormalStats(tdelta[j], tsigma_m[j])
                print("%s\t%04d\t%0.4f\t\t%0.4f"%(j.ljust(sep), tnum_points[j], tdelta[j], tsigma_m[j]))
                tpoint_test = tpoints[j]
                for s in np.unique(diff2):
                
                    if j in s:
                        l = s.replace(j,'')
                        mn2, var2 = getNormalStats(tdelta[j], tsigma_m[j])
                        N2 = tnum_points[j]
                        stat, p = scipy.stats.levene(points[s], tpoint_test, center='mean')
                        print("%s\t%04d\t%0.4f\t\t%0.4f\t\t%0.4f"%(l.ljust(sep), num_points[s], delta[s], sigma_m[s], p))
                        
                        
                    #mn1, var1 = getNormalStats(delta[s], sigma_m[s])
                    #N1 = num_points[s]
    
                
                #p = scipy.stats.f.cdf(var2/var1, N2, N1) if var1 > var2 else scipy.stats.f.cdf(var1/var2, N1, N2)
                
                
                
                print("\n")
    
    
    '''
    def buildLatexTable(series, material_data):
        txt = '\begin{table}[!h]\n'
        if 'material' in series.lower():
            txt = txt  + '\caption[Properties of %s]{Properties of %s ~\cite{Luo:FRA2019}.}'
        else:
            txt = txt  + '\caption[Properties of %s Materials]{Properties of %s ~\cite{Luo:FRA2019}.}'
    

    \centering
    \begin{tabular}{|l|p{1.2cm}|p{1.2cm}|p{1.2cm}|p{1.2cm}|p{1.2cm}|}
    \hline
    Property                                    &     Acrylic   &     PC Blend  &  Phenolic Resin FRP & PVC Blend  & Vinyl Ester Resin FRP \\ \hline
    Conductivity    ($\mathrm{W/(m\cdot K)}$)    &  0.24         &     0.89      &  0.35               &  0.29      &  0.54 \\ \hline
    Density        ($\mathrm{kg/m^{3}}$)        &  1178         &     1320      &  1846               &  1314      &  1600  \\ \hline
    Emissivity                                  &  0.88         &     0.93      &  0.89               &  0.88      &  0.89   \\ \hline
    Specific Heat  ($\mathrm{kJ/(kg\cdot C}$)   &  0.62         &     0.67      &  1.29               &  0.42      &  0.94 \\ \hline
    Thickness ($\mathrm{mm}$)                   &  4.5          &     4.5       &  3.3                &  3.3       &  4.5 \\ \hline
    Ignition Temperaure ($\mathrm{^{\circ}C}$)  &  390.0        &     453.7     &  520.0              &  499.5     &  579.6 \\ \hline
    Reference Heat Flux ($\mathrm{kW/m^{2}}$)   &  83.7         &     75.1      &  73.0               &  76.1      &  73.8 \\ \hline
    \end{tabular}
    \label{Properties_JH_FRA_Materials_polymers}
    \end{table}
    '''
    
