import pyfdstools as fds
import os, glob, sys, argparse
import numpy as np

if __name__ == '__main__':
    
    args = sys.argv
    systemPath = os.path.dirname(os.path.abspath(__file__))
    
    filename = 'export_to_vtk.py'
    
    availableArguments = "--help\n"
    availableArguments = availableArguments + "--all Outputs all FDS outputs to VTK\n"
    availableArguments = availableArguments + "--binary Selects whether VTK files are ASCIi or Binary\n"
    availableArguments = availableArguments + "--bnde Exports BNDE files to VTK, default False\n"
    availableArguments = availableArguments + "--bndf Exports BNDF files to VTK, default False\n"
    availableArguments = availableArguments + "--chid CHID of the FDS case, defaults to basename of the first .fds file in dir\n"
    availableArguments = availableArguments + "--clean Deletes VTK and related files from dir prior to run\n"
    availableArguments = availableArguments + "--indir Directory of the FDS output files, defaults to current working directory\n"
    availableArguments = availableArguments + "--outdir Directory of the VTK files, defaults to current working directory\n"
    availableArguments = availableArguments + '--times Establishes the time interval to export data [start time, end time, NFRAMES]\n'
    availableArguments = availableArguments + "--prt5 Exports PRT5 files to VTK, default False\n"
    availableArguments = availableArguments + "--sl2d Exports SL2D files to VTK, default False\n"
    availableArguments = availableArguments + "--sl3d Exports SL3D files to VTK, default False\n"
    availableArguments = availableArguments + "--s3d Exports S3D files to VTK, default False\n"
    availableArguments = availableArguments + "--stl Generates STL file based on geometry, default True"
    
    parser = argparse.ArgumentParser()
    parser.add_argument('call')
    parser.add_argument('--all', action='store_true', default=False, help='Outputs all FDS outputs to VTK, default False')
    parser.add_argument('--binary', action='store_true', default=False, help='Selects whether VTK files are ASCII or Binary, default ASCII')
    parser.add_argument('--bnde', action='store_true', default=False, help='Exports BNDE files to VTK, default False')
    parser.add_argument('--bndf', action='store_true', default=False, help='Exports BNDF files to VTK, default False')
    parser.add_argument('--chid', nargs=1, default=[], help='CHID of the FDS case, defaults to basename of the first .fds file in dir')
    parser.add_argument('--clean', action='store_true', help='Deletes VTK and related files from dir prior to run')
    parser.add_argument('--indir', nargs=1, default=[], help='Directory of the FDS output files, defaults to current working directory')
    parser.add_argument('--outdir', nargs=1, default=[], help='Directory of the VTK output files, defaults to current working directory')
    parser.add_argument('--prt5', action='store_true', default=False, help='Exports PRT5 files to VTK, default False')
    parser.add_argument('--sl2d', action='store_true', default=False, help='Exports SL2D files to VTK, default False')
    parser.add_argument('--sl3d', action='store_true', default=False, help='Exports SL3D files to VTK, default False')
    parser.add_argument('--s3d', action='store_true', default=False, help='Exports S3D files to VTK, default False')
    parser.add_argument('--stl', action='store_true', default=False, help='Disables export of STL file, default True')
    parser.add_argument('--times', nargs=3, default=[], help='Establishes the time interval to export data [start time, end time, NFRAMES]')
    #parser.add_argument('--help', action='store_true', default=False)
    
    if len(args) == 1:
        sys.exit("Available arguments\n%s"%(availableArguments))
    
    cmdargs = parser.parse_args(args)
    
    if len(cmdargs.indir) == 0:
        print("Warning, no input directory was provided, assuming the script directory.")
        rdir = os.getcwd() + os.sep
    else:
        if os.path.isabs(cmdargs.indir[0]):
            rdir = cmdargs.indir[0] + os.sep
        else:
            rdir = os.path.abspath(cmdargs.indir[0]) + os.sep
    if len(cmdargs.outdir) == 0:
        print("Warning, no input directory was provided, assuming the script directory.")
        odir = os.getcwd() + os.sep
    else:
        if os.spath.isabs(cmdargs.outdir[0]):
            odir = cmdargs.outdir[0] + os.sep
        else:
            odir = os.path.abspath(cmdargs.outdir[0]) + os.sep
    
    if len(cmdargs.chid) == 0:
        print("Warning, no chid was provided, assuming the namespace of the fds file.")
        fdsFiles = fds.getFileList(rdir, '', 'fds')
        if len(fdsFiles) == 0:
            sys.exit('No fds files found in %s.'%(rdir))
        if len(fdsFiles) > 1:
            print("Warning, multiple fds files found, assuming %s is correct."%(fdsFiles[0]))
        chid = os.path.basename(fdsFiles[0]).replace('.fds','')
        print('CHID = %s'%(chid))
    else:
        chid = cmdargs.chid[0]
    
    if cmdargs.clean:
        for f in glob.glob(odir + '*.vtr'):
            os.remove(f)
        for f in glob.glob(odir + '*.pvtr'):
            os.remove(f)
        for f in glob.glob(odir + '*.vti'):
            os.remove(f)
        for f in glob.glob(odir + '*.pvti'):
            os.remove(f)
        for f in glob.glob(odir + '*.vtp'):
            os.remove(f)
        for f in glob.glob(odir + '*.pvtp'):
            os.remove(f)
    
    if cmdargs.all:
        cmdargs.sl2d = True
        cmdargs.sl3d = True
        cmdargs.s3d = True
        cmdargs.bndf = True
        cmdargs.bnde = True
        cmdargs.prt5 = True
    
    if len(cmdargs.times) == 0:
        outTimes = None
    else:
        outTimes = np.linspace(float(cmdargs.times[0]), float(cmdargs.times[1]), int(cmdargs.times[2]))
    
    if not cmdargs.stl:
        fds.obstToStl(rdir, chid)

    if cmdargs.sl2d:
        print("Converting SL2D to VTK...")
        fds.exportSl2dDataToVtk(chid, rdir, outtimes=outTimes, binary=cmdargs.binary)
    
    if cmdargs.sl3d:
        print("Converting SL3D to VTK...")
        fds.exportSl3dDataToVtk(chid, rdir, outtimes=outTimes, binary=cmdargs.binary)
        
    if cmdargs.s3d:
        print("Converting S3D to VTK...")
        fds.exportS3dDataToVtk(chid, rdir, outtimes=outTimes, binary=cmdargs.binary)
    
    if cmdargs.bndf:
        print("Converting BNDF to VTK...")
        fds.exportBndfDataToVtk(chid, rdir, outtimes=outTimes, binary=cmdargs.binary)
    
    if cmdargs.bnde:
        print("Converting BNDE to VTK...")
        fds.exportBndeDataToVtk(chid, rdir, outtimes=outTimes, binary=cmdargs.binary)
        
    if cmdargs.prt5:
        print("Converting PRT5 to VTK...")
        fds.exportPrt5DataToVtk(chid, rdir, outtimes=outTimes, binary=cmdargs.binary)
        

