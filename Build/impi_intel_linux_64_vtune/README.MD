## Intel VTune Amplifier

This folder contains the optimal build for Intel VTune Amplifier, which can locate hotspots and the performance of both parallel and serial code. This helps for targeting optimizations.

### Required Setup

1. Source amplxe-vars.sh, from the installation folder for Inspector. For example:
```
source /opt/intel19/vtune_amplifier_2019/amplxe-vars.sh
```
2. Compile the `vtune` version of FDS using the script `make_fds.sh` in this directory. The relevant compiler options are listed [here](https://software.intel.com/en-us/vtune-amplifier-help-compiler-switches-for-performance-analysis-on-linux-targets).

### Recommended Steps

For useful analysis capabilities, using amplxe-gui on the platform used for collection is ideal, or at least having access to the used source and executable on the machine where Amplifier is installed.

To setup, use `amplxe-gui` before collecting data. X11 forwarding is necessary if logging in to a remote cluster. Inside the GUI, use the new project option to create a project for FDS, selecting the 'vtune' version's executable as the target. When prompted in analysis, select the 'vtune' version's folder for binary files, and firemodels/fds/Source for source files. To use the project, place collection results, obtained later, into the generated project folder.

### Collection

The base command used on one's platform to run FDS is a proper starting point. `mpiexec fds [test case]` is a common input, and thus we'll use it for the model here.

Before the fds executable in the run command, such as between mpiexec and the executable, place `amplxe-cl -collect <analysis-type> [Optional actions] -- `. Most frequently, analysis-type can be hpc-performance. This focuses on serial and parallel performance in a high-performance computing environment, i.e. clusters. You can look at other analysis types [here](https://software.intel.com/en-us/vtune-amplifier-help-running-command-line-analysis)

Thus, a possible input could be:
```
mpiexec -np 1 amplxe-cl -collect hpc-performance -- $HOME/firemodels/fds/Build/impi_intel_linux_64_vtune/fds_impi_intel_linux_64_vtune simple_test.fds
```

Alternatively, using `qfds.sh -a [result_directory]` will run hpc-performance, automatically. 

### Analysis

#### Command Line

To gather data from the command line, use:
```
amplxe -cl -report <report_type> -r <result_dir>
```
This will display results according to report_type to standard out. `summary` tends to be a good place to start, and you can find them all [here](https://software.intel.com/en-us/vtune-amplifier-help-generating-command-line-reports)


#### Graphical User Interface (Recommended)

Due to the large amount of data collected in any Amplifier run, using the `amplxe-gui` is highly recommended With results placed/generated in the project folder, Amplifier can display the results in a convenient GUI, instead of necessitating parsing the output on the command line.

You can access the GUI on a cluster that has it installed by running:
```
amplxe-gui &
```
This runs it without exclusivity, so you can continue using your shell.

In the GUI, there are several display options under 'Bottom-up' in a results screen. When optimizing existing OpenMP directives, focus on using options that highlight OpenMP regions. Otherwise, focus on function or source file views. Hovering over data labels normally gives a cursory explanation of what is represented.
