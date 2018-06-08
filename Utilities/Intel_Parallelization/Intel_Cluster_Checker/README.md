# Intel Cluster Checker

In order to use Intel Cluster Checker, follow the user guide [here](https://software.intel.com/en-us/cluster-checker-user-guide-2019-beta). Included in this directory is an example nodefile. Note that nodes that were not operation at time of writing are commented out in the file, meaning that the nodefile should be updated to reflect which nodes are operational when the tool is used.

**Quickstart Guide**

1. Install Intel Cluster Checker and run '''source /opt/intel/clck/201n/bin/clckvars.sh''' after install (note that this path will be different if the standard install path is not used for the tool).
2. Create a nodefile based on the current status of the cluster.
3. Run '''clck -f nodefile'''.
