## Intel Trace Collector and Analyzer

 The Intel Trace Collector and Analyzer are two separate programs with a single purpose---to enable you to visualize the work flow of each MPI process of an FDS simulation. The Trace Collector is essentially built into the FDS executable via a compiler option (`-tcollect`), and it outputs a trace (`.stf`) file at the end of the FDS job. The Trace Analyzer is a visualization tool that reads the trace file and displays its contents graphically.

 For details on the Intel Trace Collector, read the [manual](https://software.intel.com/sites/default/files/intel-trace-collector-2018-user-and-reference-guide.pdf). For details on the Intel Trace Analyzer, read the [manual](https://software.intel.com/en-us/ita-user-and-reference-guide).
 
 To enable Trace Collector and Analyzer, add the following line to your `.bashrc` file:
 ```
 source /opt/intel19/parallel_studio_xe_2019/psxevars.sh
 ```
The main consideration in tracing FDS is that the trace file can become enormous if you run a long job and trace each and every function and subroutine call. To prevent this, there is a configuration file called `fds_trace.conf` in this directory that contains a list of the main subroutines called in FDS. Only these subroutines are traced, keeping the trace file to a reasonable size and enabling you to more easily visualize the work flow. 

 To output trace information, use the -r flag when invoking qfds.sh.

 To use the configuration file, add the `-c <filepath>/<configfilename>.conf` flag to qfds.sh. If using a custom script, add
```
export VT_CONFIG=<Full path to FDS repo>/Build/impi_intel_linux_64_trace/fds_trace.conf
```
to your PBS script.

Make sure that the job only runs a handful of time steps, as there's no need to make the trace file bigger than it already is. The default trace file name is `fds_trace.stf`, and once the job is finished, you can start the Trace Analyzer:
```
traceanalyzer fds_trace.stf &
```
If you get an error message, it is probably because you cannot open an X11 window on your computer. To do so, download the free program called `Xming`, start it by double-clicking, then add to your PuTTY session `Connection --> SSH --> X11 --> Enable X11 forwarding`. When you open up a linux command prompt, type `xeyes &` and if a pair of mouse-following googly eyes pop up, you're in business.

The most important graphic in the Trace Analyzer is the timeline. Get this from the `Charts` menu, `Event Timeline`. You will first see the entire timeline, but you can click and drag over shorter time intervals to see details. You will also notice that the first time you use the Trace Analyzer, everything is either colored red (MPI) or blue (Application). Go to the chart in the lower left corner and right click on the `Groups`, and choose to ungroup them. You should see the modules and subroutines you've chosen to trace. Keep ungrouping until you get down to the subroutine level. If you right-click again, you can choose to color the various routines, making it much easier to visualize. Your chosen color scheme will be saved in a file called `.itarc` in your home directory.


 ### Creating a new configuration file

 If you do not want to use the configuration file `fds_trace.conf` that is included above, you can create your own by following these steps.

 1. Compile a special version of FDS in `Build/impi_intel_linux_64_trace`. This is essentially a debug compilation with the additional compiler option `-tcollect`. 

 2. Run a very short version of the FDS job that you want to trace. Just a few time steps is sufficient. You should not use a configuration file for this run. You want to collect everything.

 3. The trace (`.stf`) file that is created in the same directory as the FDS job output has a history of every function and subroutine call.

 4. To streamline the trace analysis, create a configuration file that reduces the number of subroutines and functions to trace. To do this, run the Configuration Assistant at the command line:
```
itcconfig <trace_file.stf>
```

 5. A graphical window will pop up, called Trace Collector Configurator. Click on `Filters`, and then select the `Functions` tab in the center panel. Right click in the empty box beneath `Function Name Pattern`, select the functions and subroutines you want to trace. You might consider turning `Application` off on the first line, and then add back in those routines that you want to trace. The term `Application` means everything besides MPI calls. It is your application, i.e. FDS. 

 6. Save the configuration (`.conf`) file and exit the Configurator.

 6. Add `export VT_CONFIG=<configuration_file.conf>` to your run script, created with `qfds.sh`, or just run `qfds.sh` with the `-c <configuration_file.conf> option.

 7. Run your FDS job again, looking for the new trace (`.stf`) file.

 8. Run the Intel Trace Analyzer:
```
traceanalyzer <trace_file.stf>
```
