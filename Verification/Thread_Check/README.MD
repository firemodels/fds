## Thread Checker

This folder contains tests and scripts that ensure the health of FDS parallel code. 

### inspection.sh

Takes an fds case that starts with the word "inspector" and a number of processes, and runs the Inspect build of FDS, with Inspector running automatically. Returns basic progress notifications, and a notification that identifies whether Inspector reported clean code or found errors. The "inspector" prefix is required for automatic output analysis. Example:

```
inspection.sh -p 6 inspector_test.fds
```

### inspect_report.sh

Takes a prefix of the directory titles of an Inspector result set as an argument - i.e., the set including insp.blaze01 and insp.blaze02 could be called on with 
```
inspect_report.sh insp.blaze
```
Performs Inspector result 'problem' analysis, reporting to stdout.