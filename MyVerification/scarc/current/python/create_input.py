# -*- coding: utf-8 -*-
import sys, math

def print_head(fout, case):
   fout.write("&HEAD CHID=\'%s\', TITLE=\'2D-case for comparison of different pressure solvers\' /\n" %(case))

def print_meshes(fout, top, cell):

   topx = int(top[0])
   topy = 1
   topz = int(top[2])

   cellx = int(cell[0])
   celly = 1
   cellz = int(cell[2])

   xs      =  0.0
   xf      =  0.8

   ys      =  -0.01
   yf      =   0.01

   zs      =  0.0
   zf      =  0.8

   nx   =  cellx/topx
   ny   =  1
   nz   =  cellz/topz

   stepx = (xf-xs)/topx
   stepy = (yf-ys)/topy
   stepz = (zf-zs)/topz

   dx = (xf-xs)/topx
   dy = (yf-ys)/topy
   dz = (zf-zs)/topz

   sizex   = topx-1
   sizey   = topy-1
   sizez   = topz-1

   fout.write("&MESH IJK=%d,%d,%d, XB=%.2f,%.2f,%.2f,%.2f,%.2f,%.2f, MULT_ID=\'mesh\' /\n"  %(nx, ny, nz, xs, xs+stepx, ys, ys+stepy, zs, zs+stepz))
   fout.write("&MULT ID=\'mesh\', DX=%.2f, DY=%.2f, DZ=%.2f, I_UPPER=%d, J_UPPER=%d, K_UPPER=%d /\n\n" %(dx, dy, dz, sizex, sizey, sizez ))

def print_vel_tolerance(fout, exponent):
   if exponent > 0:
      tol = 10**(-exponent)
      fout.write("      VELOCITY_TOLERANCE      = %.0e \n" %tol)

def print_max_iterations(fout, nite):
   fout.write("      MAX_PRESSURE_ITERATIONS = %d /\n" %nite)

def dump_xyz(ptype):
   if ptype:
      fout.write("&DUMP WRITE_XYZ    = .TRUE.     \n")
   else:
      fout.write("&DUMP WRITE_XYZ    = .FALSE.    \n")
   
def dump_velocity_error(ptype):
   if ptype:
      fout.write("&DUMP VELOCTY_ERROR_FILE = .TRUE.     \n")
   else:
      fout.write("&DUMP VELOCTY_ERROR_FILE = .FALSE.    \n")
   
def dump_pl3d(q1,q,q3,q4,q5,dt_pl3d):
   fout.write("&DUMP DT_PL3D=%.2f, PLOT3D_QUANTITY= \'%s\',\'%s\',\'%s\',\'%s\',\'%s\'    \n" %(dt_pl3d,q1,q2,q3,q4,q5))

def print_discretization(fout, dtype):
   if dtype == 'structured':
      fout.write("&PRES SCARC_DISCRETIZATION    = \'STRUCTURED\'    \n")
   else:
      fout.write("&PRES SCARC_DISCRETIZATION    = \'UNSTRUCTURED\'    \n")

def print_cg_cluster(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'KRYLOV\'    \n")
   fout.write("      SCARC_KRYLOV            = \'CG\'        \n")
   fout.write("      SCARC_KRYLOV_ACCURACY   = %s            \n" %tol)
   fout.write("      SCARC_PRECON            = \'CLUSTER\'   \n")
   
def print_cg_pardiso(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'KRYLOV\'    \n")
   fout.write("      SCARC_KRYLOV            = \'CG\'        \n")
   fout.write("      SCARC_KRYLOV_ACCURACY   = %s            \n" %tol)
   fout.write("      SCARC_PRECON            = \'PARDISO\'   \n")
   
def print_cg_ssor(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'KRYLOV\'    \n")
   fout.write("      SCARC_KRYLOV            = \'CG\'        \n")
   fout.write("      SCARC_KRYLOV_ACCURACY   = %s            \n" %tol)
   fout.write("      SCARC_PRECON            = \'SSOR\'      \n")
   
def print_cg_fft(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'KRYLOV\'    \n")
   fout.write("      SCARC_KRYLOV            = \'CG\'        \n")
   fout.write("      SCARC_KRYLOV_ACCURACY   = %s            \n" %tol)
   fout.write("      SCARC_PRECON            = \'SSOR\'      \n")
   
def print_fft_cg_fft(fout, tol):
   fout.write("      SCARC_FFT               = .TRUE.        \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'KRYLOV\'    \n")
   fout.write("      SCARC_KRYLOV            = \'CG\'        \n")
   fout.write("      SCARC_KRYLOV_ACCURACY   = %s            \n" %tol)
   fout.write("      SCARC_PRECON            = \'FFT\'       \n")
   
def print_fft_cg_ssor(fout, tol):
   fout.write("      SCARC_FFT               = .TRUE.        \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'KRYLOV\'    \n")
   fout.write("      SCARC_KRYLOV            = \'CG\'        \n")
   fout.write("      SCARC_KRYLOV_ACCURACY   = %s            \n" %tol)
   fout.write("      SCARC_PRECON            = \'SSOR\'      \n")
   
def print_cluster(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'MKL\'       \n")
   fout.write("      SCARC_MKL               = \'GLOBAL\'    \n")
   fout.write("      SCARC_MKL_MTYPE         = \'SYMMETRIC\' \n")
   
def print_gmg_ssor_dir2(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'MULTIGRID\'    \n")
   fout.write("      SCARC_MULTIGRID         = \'GEOMETRIC\'    \n")
   fout.write("      SCARC_MULTIGRID_ACCURACY= %s               \n" %tol)
   fout.write("      SCARC_MULTIGRID_CYCLE   = \'V\'            \n")
   fout.write("      SCARC_MULTIGRID_LEVEL   =  2               \n")
   fout.write("      SCARC_COARSE            = \'DIRECT\'       \n")
   fout.write("      SCARC_SMOOTH            = \'SSOR\'      \n")
   fout.write("      SCARC_SMOOTH_OMEGA      =  0.8E+0          \n")
   fout.write("      SCARC_SMOOTH_ITERATIONS =   5              \n")
   
def print_gmg_ssor_dir3(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'MULTIGRID\'    \n")
   fout.write("      SCARC_MULTIGRID         = \'GEOMETRIC\'    \n")
   fout.write("      SCARC_MULTIGRID_ACCURACY= %s               \n" %tol)
   fout.write("      SCARC_MULTIGRID_CYCLE   = \'V\'            \n")
   fout.write("      SCARC_MULTIGRID_LEVEL   =  3               \n")
   fout.write("      SCARC_COARSE            = \'DIRECT\'       \n")
   fout.write("      SCARC_SMOOTH            = \'SSOR\'      \n")
   fout.write("      SCARC_SMOOTH_OMEGA      =  0.8E+0          \n")
   fout.write("      SCARC_SMOOTH_ITERATIONS =   5              \n")
   
def print_gmg_ssor(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'MULTIGRID\'    \n")
   fout.write("      SCARC_MULTIGRID         = \'GEOMETRIC\'    \n")
   fout.write("      SCARC_MULTIGRID_ACCURACY= %s               \n" %tol)
   fout.write("      SCARC_MULTIGRID_CYCLE   = \'V\'            \n")
   fout.write("      SCARC_COARSE            = \'ITERATIVE\'     \n")
   fout.write("      SCARC_SMOOTH            = \'SSOR\'      \n")
   fout.write("      SCARC_SMOOTH_OMEGA      =  0.8E+0          \n")
   fout.write("      SCARC_SMOOTH_ITERATIONS =   5              \n")

def print_gmg_ssor_ite(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'MULTIGRID\'    \n")
   fout.write("      SCARC_MULTIGRID         = \'GEOMETRIC\'    \n")
   fout.write("      SCARC_MULTIGRID_ACCURACY= %s               \n" %tol)
   fout.write("      SCARC_MULTIGRID_CYCLE   = \'V\'            \n")
   fout.write("      SCARC_MULTIGRID_LEVEL   =  3               \n")
   fout.write("      SCARC_COARSE            = \'ITERATIVE\'     \n")
   fout.write("      SCARC_SMOOTH            = \'SSOR\'      \n")
   fout.write("      SCARC_SMOOTH_OMEGA      =  0.8E+0          \n")
   fout.write("      SCARC_SMOOTH_ITERATIONS =   5              \n")
   
def print_gmg_pardiso_dir2(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'MULTIGRID\'    \n")
   fout.write("      SCARC_MULTIGRID         = \'GEOMETRIC\'    \n")
   fout.write("      SCARC_MULTIGRID_ACCURACY= %s               \n" %tol)
   fout.write("      SCARC_MULTIGRID_CYCLE   = \'V\'            \n")
   fout.write("      SCARC_MULTIGRID_LEVEL   =  2               \n")
   fout.write("      SCARC_COARSE            = \'DIRECT\'       \n")
   fout.write("      SCARC_SMOOTH            = \'PARDISO\'      \n")
   fout.write("      SCARC_SMOOTH_OMEGA      =  0.8E+0          \n")
   fout.write("      SCARC_SMOOTH_ITERATIONS =   5              \n")
   
def print_gmg_pardiso_dir3(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'MULTIGRID\'    \n")
   fout.write("      SCARC_MULTIGRID         = \'GEOMETRIC\'    \n")
   fout.write("      SCARC_MULTIGRID_ACCURACY= %s               \n" %tol)
   fout.write("      SCARC_MULTIGRID_CYCLE   = \'V\'            \n")
   fout.write("      SCARC_MULTIGRID_LEVEL   =  3               \n")
   fout.write("      SCARC_COARSE            = \'DIRECT\'       \n")
   fout.write("      SCARC_SMOOTH            = \'PARDISO\'      \n")
   fout.write("      SCARC_SMOOTH_OMEGA      =  0.8E+0          \n")
   fout.write("      SCARC_SMOOTH_ITERATIONS =   5              \n")
   
def print_gmg_pardiso_ite(fout, tol):
   fout.write("      SCARC_FFT               = .FALSE.       \n")
   fout.write("      SCARC_ACCURACY          = \'ABSOLUTE\'  \n")
   fout.write("      SCARC_METHOD            = \'MULTIGRID\'    \n")
   fout.write("      SCARC_MULTIGRID         = \'GEOMETRIC\'    \n")
   fout.write("      SCARC_MULTIGRID_ACCURACY= %s               \n" %tol)
   fout.write("      SCARC_MULTIGRID_CYCLE   = \'V\'            \n")
   fout.write("      SCARC_MULTIGRID_LEVEL   =  3               \n")
   fout.write("      SCARC_COARSE            = \'ITERATIVE\'     \n")
   fout.write("      SCARC_SMOOTH            = \'PARDISO\'      \n")
   fout.write("      SCARC_SMOOTH_OMEGA      =  0.8E+0          \n")
   fout.write("      SCARC_SMOOTH_ITERATIONS =   5              \n")
   
def print_line_devices(fout, cell, point, quantities):
   for quan in quantities:
      fout.write("&DEVC XB=%.2f,%.2f,%.2f,%.2f,%.2f,%.sf, QUANTITY = \'%s\'" %(p[0],p[1],p[2],p[3],p[4],p[5],quan))

   
def create_input(base_input, script, discret, top, cell, solver, obsttype, tolerances):
   
   nmeshes = int(top[0])*int(top[1])*int(top[2])
   if discret == 'structured':
      for tol in tolerances:
         case='s%d_%s_tol%d_%s' %(nmeshes, obsttype, tol, solver)
         print_fds_file(base_input, case, discret, tol)
         add_to_script (script, case, nmeshes)
   elif discret == 'unstructured':
      case='u%d_%s_tol%d_%s' %(nmeshes, obsttype, 0, solver)
      print_fds_file(base_input, case, discret, 0)
      add_to_script (script, case, nmeshes)
   else:
      sys.exit("Wrong discretization type %s"%discret)


def print_fds_file(base_input, case, discret, tol):
   case_file = "input/%s.fds" %(case)
   print "creating %s" % case_file
   fds  = open(case_file ,'w')
   tol_scarc=1E-08
   for line in base_input:
      line, null = line.split ("\n")
      if 'HEAD' in line:
         print_head(fds, case)
      elif 'MESH' in line:
         print_meshes(fds, top, cell)
         print_discretization(fds, discret)
         if solver == "cluster":
            print_cluster(fds, tol_scarc)
         elif solver == "cg_ssor":
            print_cg_ssor(fds, tol_scarc)
         elif solver == "cg_fft":
            print_cg_fft(fds, tol_scarc)
         elif solver == "fft_cg_fft":
            print_fft_cg_fft(fds, tol_scarc)
         elif solver == "fft_cg_ssor":
            print_fft_cg_ssor(fds, tol_scarc)
         elif solver == "cg_cluster":
            print_cg_cluster(fds, tol_scarc)
         elif solver == "cg_pardiso":
            print_cg_pardiso(fds, tol_scarc)
         elif solver == "gmg_ssor":
            print 'hallo'
            print_gmg_ssor(fds, tol_scarc)
         elif solver == "gmg_pardiso_dir2":
            print_gmg_pardiso_dir3(fds, tol_scarc)
         elif solver == "gmg_pardiso_dir3":
            print_gmg_pardiso_dir2(fds, tol_scarc)
         elif solver == "gmg_pardiso_ite":
            print_gmg_pardiso_ite(fds, tol_scarc)
         elif solver == "gmg_ssor_dir2":
            print_gmg_ssor_dir2(fds, tol_scarc)
         elif solver == "gmg_ssor_dir3":
            print_gmg_ssor_dir3(fds, tol_scarc)
         elif solver == "gmg_ssor_ite":
            print_gmg_ssor_ite(fds, tol_scarc)
         if discret == "structured":
            print_vel_tolerance(fds, tol)
            if tol == 0:
               print_max_iterations(fds, 10)
            else:
               print_max_iterations(fds, 1000)
         else:
            print_max_iterations(fds, 10)
      #elif 'TAIL' in line:
      #   print_line_devices(fds, cell)
      #   fds.write("%s\n" % line)
      else:
         fds.write("%s\n" % line)
   fds.close()
   
   
def add_to_script (script, case, nmeshes):
   script.write("echo \"------> Running %s\"\n" %case)
   script.write("mkdir cases/%s \n" %case)
   script.write("cd cases/%s \n" %case)
   script.write("cp ../../input/%s.fds . \n" %case)
   script.write("mpirun -np %d $BUD/fds_mpi_intel_osx_64_db %s.fds  \n" %(nmeshes, case))
   script.write("cd ../.. \n\n")
   

obsttype = sys.argv[1]

#discretizations = ['structured','unstructured']
discretizations = ['structured']
topologies      = [[1,1,1],[2,1,2]]
cells           = [[16,1,16]]
#solvers         = ['fft','cg_ssor','cg_fft','gmg_ssor','fft_cg_ssor']
solvers         = ['gmg_ssor']
tolerances      = [2,6,10]

base_name   = "base/s1_%s_fft.fds" %obsttype
script_name = "scripts/do_%s.sh"  %obsttype
script  = open(script_name ,'w')
   
f = open(base_name ,'r')
base_input  = f.readlines()
f.close()

for discret in discretizations:
   for top in topologies[:]:
      for cell in cells[:]:
         for solver in solvers[:]:
            if discret == 'unstructured' and 'fft' in solver: continue
            create_input (base_input, script, discret, top, cell, solver, obsttype, tolerances) 

script.close()
