# Intel Optimization Tools

This folder contains a few configuration files that work with the Intel Inspector. This tool is mainly useful on Linux, with multiple processes and threads.

## Intel Inspector

The [Intel Inspector](https://software.intel.com/en-us/node/622387) can help detect improperly coded OpenMP directives.

### Collection

1. Compile the `inspect` version of FDS using the script `make_fds.sh` in this directory. The relevant compiler options are listed [here](https://software.intel.com/en-us/inspector-user-guide-linux-building-applications).
2. Run fds with 'qfds' using -x [result directory] - this runs in Intel's ti2 form with a stack depth of 32. You can read more [here](https://software.intel.com/en-us/inspector-user-guide-linux-threading-error-analysis-types)
3. From firemodels/fds/Verification/Thread_Check, you can get inspect_report.sh. This script can handle the multiple files generated in any Inspector run with more than one process. 

### Analysis

From firemodels/fds/Verification/Thread_Check, you can get inspect_report.sh. This script can handle the multiple files generated in any Inspector run with more than one process.

Alternatively, you can download the gui version of [Intel Inspector](https://software.intel.com/en-us/inspector), and use it on the inspxe files in the created directories.

## Intel Cluster Checker

The Intel Cluster Checker is a diagnostic tool that checks the overall health of your compute cluster. To use it, first consult the [User's Guide](https://software.intel.com/en-us/cluster-checker-user-guide-2019-beta). In brief, do the following:

   1. Install Intel Cluster Checker and run `source /opt/intel19/clck_latest/clckvars.sh`. Your path might be slightly different.
   2. Create a `nodefile`, which is just a text file with a list of the cluster node names, one per line.
   3. Run `clck -f nodefile`

