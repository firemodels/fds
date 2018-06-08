# Intel Trace Collector Configuration

For details on Intel Trace Collector, view the manual [here](https://software.intel.com/sites/default/files/intel-trace-collector-2018-user-and-reference-guide.pdf). The example configuration file in this directory was made and used according to section 2.4.1 of the manual. This example filters collection of the functions that create the most trace data for a small (.1 second duration) simulation. For details on creating a .conf file, see section 2.5 of the manual.

**Quickstart Guide**

1. Install Intel Trace Analyzer and Collector.
2. Create a configuration file manually or by launching Trace Analyzer with '''traceanalyzer''', loading a file, and then using the configuration assistant.
3. Run '''export VT_CONFIG=/<full file path>/myconfiguration.conf'''.
4. Trace application normally using Intel Trace Analyzer.
