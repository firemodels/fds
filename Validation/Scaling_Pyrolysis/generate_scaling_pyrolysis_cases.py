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

def buildFdsFile(chid, coneExposure, e, k, rho, cp, Tign, d, time, hrrpua, tend, HFs):
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
    
    if Tign == 'Calculated':
        Tign = 1000
        DT_DEVC = 0.1 
        tempOutput = '.TRUE.'
    else:
        DT_DEVC = 10.
        tempOutput = '.FALSE.'
    
    txt = "&HEAD CHID='%s', /\n"%(chid)
    txt = txt+"&TIME DT=1., T_END=%0.1f /\n"%(tend)
    txt = txt+"&DUMP DT_CTRL=%0.1f, DT_DEVC=%0.1f, DT_HRR=%0.1f, SIG_FIGS=4, SIG_FIGS_EXP=2, /\n"%(DT_DEVC, DT_DEVC, DT_DEVC)
    txt = txt+"&MISC SOLID_PHASE_ONLY=.TRUE., TMPA=27., /\n"
    txt = txt+"&MESH ID='MESH', IJK=3,3,3, XB=0.,0.3,0.,0.3,0.,0.3, /\n"
    txt = txt+"&REAC ID='PROPANE', FUEL='PROPANE', /\n"
    txt = txt+"&MATL ID='BACKING', CONDUCTIVITY=0.2, DENSITY=585., EMISSIVITY=1., SPECIFIC_HEAT=0.8, /\n"
    txt = txt+"&MATL ID='SAMPLE', CONDUCTIVITY=%0.4f, DENSITY=%0.1f, EMISSIVITY=%0.4f, SPECIFIC_HEAT=%0.4f, /\n"%(k, rho, e, cp)
    
    for i in range(0, len(time)):
        txt = txt+"&RAMP ID='CONE-RAMP', T=%0.1f, F=%0.4f, /\n"%(time[i]-time[0], hrrpua[i])
    
    y = -0.05
    for i, hf in enumerate(HFs):
        hf_ign = estimateHrrpua(coneExposure, hrrpua_ref, hf)
        if i%3 == 0: y = y + 0.1
        XYZ = [((i % 3))*0.1+0.05, y, 0.0]
        XB = [XYZ[0]-0.05, XYZ[0]+0.05, XYZ[1]-0.05, XYZ[1]+0.05, 0.0,0.0]
        
        txt = txt+"&SURF ID='SAMPLE-%02d', EXTERNAL_FLUX=1., "%(hf)
        txt = txt+"HEAT_TRANSFER_COEFFICIENT=0., HEAT_TRANSFER_COEFFICIENT_BACK=15., "
        txt = txt+"HRRPUA=1., IGNITION_TEMPERATURE=%0.1f, MATL_ID(1:2,1)='SAMPLE','BACKING', "%(Tign)
        txt = txt+"RAMP_EF='IGNITION_RAMP-%02d', RAMP_Q='CONE-RAMP', "%(hf)
        txt = txt+"REFERENCE_HEAT_FLUX=%0.4f, REFERENCE_HEAT_FLUX_TIME_INTERVAL=1.,"%(qref)
        txt = txt+'THICKNESS(1:2)=%0.4f,%0.4f, /\n'%(d, 0.0254/2)
        
        txt = txt+"&RAMP ID='IGNITION_RAMP-%02d', T=%0.1f, F=%0.4f, DEVC_ID='IGNITION_DEVC-%02d', /\n"%(hf, 0.0, hf, hf)
        txt = txt+"&RAMP ID='IGNITION_RAMP-%02d', T=%0.1f, F=%0.4f, /\n"%(hf, 1.0, hf_ign)
        
        txt = txt+"&VENT ID='SAMPLE-%02d', SURF_ID='SAMPLE-%02d', XB="%(hf, hf)
        for x in XB:
            txt = txt+"%0.4f,"%(x)
        txt = txt+', /\n'
        
        txt = txt+"&DEVC ID='WALL TEMPERATURE-%02d', INITIAL_STATE=.FALSE., IOR=3, OUTPUT=%s, "%(hf, tempOutput)
        txt = txt+"QUANTITY='WALL TEMPERATURE', SETPOINT=%0.1f, XYZ=%0.4f,%0.4f,%0.4f, /\n"%(Tign, XYZ[0], XYZ[1], XYZ[2])
        
        txt = txt+"&CTRL ID='IGNITION-CTRL-%02d', FUNCTION_TYPE='ANY', INPUT_ID='WALL TEMPERATURE-%02d', /\n"%(hf, hf)
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
        

if __name__ == "__main__":
    
    fdsdir, fdscmd = findFds()
    
    specificationFile = "scaling_pyrolysis_cone_cases.csv"
    figoutdir = "figures"
    
    if figoutdir is not None:
        if os.path.exists(figoutdir) is not True: os.mkdir(figoutdir)
        import matplotlib.pyplot as plt
    
    specificationFile = pd.read_csv(specificationFile)
    material_output_data = dict()
    for i in range(0, specificationFile.shape[0]):
        
        # Check for run code
        code = specificationFile.iloc[i]['Code']
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
        targetTimes, HRRs_interp = interpolateExperimentalData(times, HRRs, targetDt=15, filterWidth=False)
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
            hrrpua_ref = getRepresentativeHrrpua(hrrs_trimmed, factor=0.5)
        except:
            hrrpua_ref = getRepresentativeHrrpua(HRRs, factor=0.5)
            print("Warning: Failed to get representative HRRPUA for material %s on trimmed data. Using full HRR curve."%(material))
            hrrs_trimmed = HRRs
            times_trimmed = times
        qref = estimateExposureFlux(coneExposure, hrrpua_ref)
        
        # Set chid
        chid = material+"_cone"
        if len(chid) > 50:
            chid = chid[:50]
        
        if calculateIgnitionTemperature:
            # Generate initial FDS input file to calculate ignition temperature
            tend = times.max() + 300 # Arbitrary, overwritten in actual calculation after finding Tign
            txt = buildFdsFile(chid, coneExposure, emissivity, conductivity, density, 
                                   specific_heat, 'Calculated', thickness, times_trimmed, hrrs_trimmed,
                                   tend, fluxes)
            
            with open("%s%s%s.fds"%(workingDir, os.sep, chid), 'w') as f:
                f.write(txt)
            
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
            print("Material %s ignition temperature %0.1f reference flux %0.1f"%(material, Tign, qref))
        else:
            tend = times.max()*1.5 + 300
            print("Material %s ignition temperature %0.1f reference flux %0.1f"%(material, Tign, qref))
        
        try:
            tMaxes = [exp_data[key].max() for key in list(exp_data.keys()) if 'Time' in key]
            tend = max(tMaxes) * 1.5
        except:
            tend = times.max()*1.5
        # Generate fds input file with updated ignition temperature
        txt = buildFdsFile(chid, coneExposure, emissivity, conductivity, density, 
                               specific_heat, Tign, thickness, times_trimmed, hrrs_trimmed,
                               tend, fluxes)
        with open("%s%s%s.fds"%(workingDir, os.sep, chid), 'w') as f:
            f.write(txt)
        runModel(workingDir, chid+".fds", 1, fdsdir, fdscmd, printLiveOutput=False)
        data = load_csv(workingDir, chid)
        
        # Plot results
        if figoutdir is not None:
            fig = plotResults_exp(data, exp_data, fluxes, validationTimeColumns, validationHrrpuaColumns, lw=3, fs=16)
            fig.savefig(os.path.join(figoutdir, chid+'.png'), dpi=300)
            plt.close()
        
        # Copy results to directories for building guide
        shutil.copy(os.path.join(workingDir, chid+"_devc.csv"), resultDir)
        shutil.copy(os.path.join(workingDir, chid+"_git.txt"), resultDir)
        shutil.copy(os.path.join(workingDir, chid+".fds"), inputFileDir)
        
        material_output_data[material] = dict()
        material_output_data[material]['conductivity'] = conductivity
        material_output_data[material]['density'] = density
        material_output_data[material]['emissivity'] = emissivity
        material_output_data[material]['specific_heat'] = specific_heat
        material_output_data[material]['thickness'] = thickness
        material_output_data[material]['Tign'] = Tign
        material_output_data[material]['qref'] = qref
        
    material_output_data = pd.DataFrame(material_output_data)
    material_output_data.to_csv('material_output_data.csv')
