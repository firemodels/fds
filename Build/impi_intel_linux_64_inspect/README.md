## Intel Inspector

This folder contains a few configuration files that work with the [Intel Inspector](https://software.intel.com/en-us/node/622387), a tool that can help detect improperly coded OpenMP directives.

### Required setup
1. Source inspxe-vars.sh, from the installation folder for Inspector. For example:
```
source /opt/intel19/inspector_2019/inspxe-vars.sh
```
2. Compile the `inspect` version of FDS using the script `make_fds.sh` in this directory. The relevant compiler options are listed [here](https://software.intel.com/en-us/inspector-user-guide-linux-building-applications).

3. For threading error analysis, you must use at least 2 OpenMP threads.
### Recommended Steps

For useful analysis capabilities, using inspxe-gui on the platform used for collection is ideal. By doing so, Inspector can have direct access to relevant source files, and thus give more useful results that do not require investigation with a secondary code viewer. 

To setup, use inspxe-gui before collecting data. X11 forwarding is necessary if logging in to a remote cluster. Inside the GUI, use the new project option to create a project for FDS, selecting the 'inspect' version's executable as the target, the 'inspect' version's folder for binary files, and firemodels/fds/Source for source files. To use the project, place collection results, obtained later, into the generated project folder.

### Collection

The base command used on one's platform to run FDS is a proper starting point. `mpiexec fds [test case]` is a common input, and thus we'll use it for the model here.

Before the fds executable in the run command, such as between mpiexec and the executable, place `inspxe-cl -collect <analysis-type> (optional options) -- `. Most frequently, analysis-type can be ti2. ti2 is a balanced analysis for threading errors, identifying more errors than ti1 in less time than ti3. If you need to use another analysis type, you can read more [here](https://software.intel.com/en-us/inspector-user-guide-linux-threading-error-analysis-types)

Thus, a possible input could be:
```
mpiexec -np 1 inspxe-cl -collect ti2 -- $HOME/firemodels/fds/Build/impi_intel_linux_64_inspect/fds_impi_intel_linux_64_inspect simple_test.fds
```

One optional knob we use is: `-knob stack-depth=32` This ensures that, regardless of how many calls are made by a function, even deep errors can be found. You can find more options [here](https://software.intel.com/en-us/inspector-user-guide-linux-inspxe-cl-actions-options-and-arguments)

Alternatively, using `qfds.sh -x [result_directory]` will run ti2 with a stack depth of 32, automatically. 

### Analysis

#### Command Line

To gather data from the command line, use:
```
inspxe-cl -report <report_type> -r <result_dir>
```

This will display results according to report_type to standard out. `problems` tends to be a suitable method, but you can find them all [here](https://software.intel.com/en-us/inspector-user-guide-linux-report)

Alternatively firemodels/fds/Verification/Thread_Check, you can get inspect_report.sh. This script can handle the multiple files generated in any Inspector run with more than one process.
```
./inspect_report.sh -n <directory_prefix>
```

#### Graphical User Interface

Alternatively, you can use the recommended inspxe-gui. It can be run directly on a inspxe file, or, if a project was created, the result directories can be generated/placed in the project's folder, and viewed by clicking the names of the results on the left.

