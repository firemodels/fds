## Detecting OpenMP Thread Race Conditions

To check for OpenMP race conditions, add `-fsanitize=thread` to the arguments (`FFLAGS`) for `impi_intel_linux_openmp_db` and add the lines
```
export OMP_TOOL_LIBRARIES='libarcher.so'
export TSAN_OPTIONS='ignore_noninstrumented_modules=1'
```
to the bash script containing the PBS/Torque or Slurm directives.

Then run a short, small job and race errors should be written to the `.err` file.

