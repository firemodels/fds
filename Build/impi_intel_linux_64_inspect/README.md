## Intel Inspector

This folder contains a few configuration files that work with the [Intel Inspector](https://software.intel.com/en-us/node/622387), a tool that can help detect improperly coded OpenMP directives.

### Required Setup

1. Run the setup script `inspxe-vars.sh` that is located in the installation folder for the Intel compiler. For example:
```
source /opt/intel19/inspector_2019/inspxe-vars.sh
```
2. Compile a specially instrumented FDS executable using the script `make_fds.sh` in this directory. The relevant compiler options are listed [here](https://software.intel.com/en-us/inspector-user-guide-linux-building-applications). Consult the `makefile` to see what options are used for FDS.

3. For threading error analysis, you must use at least 2 OpenMP threads.

### Collection Procedure

Intel Inspector can be run from the command line or via a graphical user interface (GUI). To invoke the GUI, first add X11 forwarding to your PuTTY (or equivalent) session if logging into a remote linux cluster. Then type `inspxe-gui` to open the GUI. Inside the GUI, use the new project option to create a project for FDS, selecting `fds_impi_intel_linux_64_inspect` as the target, the current directory for binary files, and `firemodels/fds/Source` for source files. To use the project, place collection results, obtained later, into the generated project folder.

To run Inspector from the command line, type the following:
```
mpiexec -np 1 inspxe-cl -collect ti2 -- $HOME/firemodels/fds/Build/impi_intel_linux_64_inspect/fds_impi_intel_linux_64_inspect simple_test.fds
```
The analysis-type, `ti2`, generates a reasonably balanced analysis for threading errors, identifying more errors than `ti1` in less time than `ti3`. If you need to use another analysis type, you can read more [here](https://software.intel.com/en-us/inspector-user-guide-linux-threading-error-analysis-types). 

The one option we use is `-knob stack-depth=32`. This ensures that, regardless of how many calls are made by a function, even deep errors can be found. You can find more options [here](https://software.intel.com/en-us/inspector-user-guide-linux-inspxe-cl-actions-options-and-arguments).

Rather than issuing the `mpiexec` call at the command line, you can invoke `qfds.sh -x [result_directory] ...` which will run a `ti2` analysis with a stack depth of 32, automatically. The script `inspection.sh [case_name.fds]` in `fds/Verification/Thread_Check` is another automated script which reports successes and failures to stdout. You can run this script in the background by adding &, and redirect the output to a file. 

### Analysis Procedure using the Command Line option

To analyze the results using command line formm of Inspector, type:
```
inspxe-cl -report <report_type> -r <result_dir>
```
This will display results according to `report_type` to standard out. `problems` is a suitable `report_type`, but you can find them all [here](https://software.intel.com/en-us/inspector-user-guide-linux-report).

Alternatively, `fds/Verification/Thread_Check` contains a script called `inspect_report.sh`. This script can handle the multiple files generated in any Inspector run with more than one process.
```
./inspect_report.sh -n <directory_prefix>
```

### Analysis Procedure using the Graphical User Interface

Alternatively, you can open the GUI via `inspxe-gui`. It can be run directly on a `inspxe` file, or, if a project was created, the result directories can be generated/placed in the project's folder, and viewed by clicking the names of the results on the left.

