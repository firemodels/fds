#!/bin/bash 
# For use in spark.el.nist.gov:
rm -f loop3d_gnu loop3d_ifort loop3d_ifx *.txt
# GNU - gfortran:
module load gcc/13.1.1
echo "Compiling loop3d_gnu.."
gfortran -m64 -O3 -cpp  main.f90 -o loop3d_gnu
echo "Submitting loop3d_gnu job.."
sbatch submit_gnu.slog

#Intel:
module unload gcc/13.1.1
module load tbb/latest compiler-rt/latest compiler/latest
# ifort:
echo "Compiling loop3d_ifort.."
ifort -m64 -O2 -ipo -fpp -DIFORT -diag-disable=10448 main.f90 -o loop3d_ifort
echo "Submitting loop3d_ifort job.."
sbatch submit_ifort.slog

# ifx:
echo "Compiling loop3d_ifx.."
ifx -m64 -O2 -ipo -fpp -DIFX main.f90 -o loop3d_ifx
echo "Submitting loop3d_ifx job.."
sbatch submit_ifx.slog

