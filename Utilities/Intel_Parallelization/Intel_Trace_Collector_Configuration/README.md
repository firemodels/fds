# Intel Trace Collector

The Intel Trace Collector is essentially a compiler option that instructs the FDS executable to record a timeline of each MPI process as it makes calls to subroutines and functions. For details on Intel Trace Collector, read the [manual](https://software.intel.com/sites/default/files/intel-trace-collector-2018-user-and-reference-guide.pdf). 

## Creating the Configuration File

1. Compile a special version of FDS in `Build/impi_intel_linux_64_inspect`. This is essentially a debug compilation with the additional compiler option `-tcollect`. This option will force FDS to create a trace (`.stf`) file after running a job.

2. Run a very short version of the FDS job that you want to trace. Just a few time steps is sufficient.

3. The trace (`.stf`) file that is created in the same directory as the FDS job output has a history of every function and subroutine call, which makes it very clumsy to work with.

4. To streamline the trace analysis, you want to create a configuration file that reduces the number of subroutines and functions to trace. To do this, run the Configuration Assistant at the command line:
```
itcconfig <trace_file.stf>
```

5. A graphical window will pop up, called Trace Collector Configurator. Click on `Filters`, and then select the `Functions` tab in the center panel. Right click in the empty box beneath `Function Name Pattern`, select the functions and subroutines you want to trace. You might consider turning everything off on the first line, and then add back in those routines that you want to trace.

6. Save the configuration (`.conf`) file.

6. Add `export VT_CONFIG=<configuration_file.conf>` to your run script, created with `qfds.sh`.

7. Run your FDS job again, looking for the new trace (`.stf`) file.

8. Run the Intel Trace Analyzer:
```
traceanalyzer <trace_file.stf>
```
