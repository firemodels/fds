# PR checks

- Listed PR authors may merge and start builds before the virus scan finishes.
- Other authors must pass the virus scan before merging or starting builds.
- Merging never cancels or skips the scan or selected builds.

`pr-security.yml` checks the author and runs the virus scan.
`pr-builds.yml` runs the existing build/test jobs.
The author list is `trusted-developers.json`.

Run the regression tests from the repository root:

```sh
python3 -B -m unittest discover -s .github/security -p 'test_*.py'
```

See the [Configuration Management Plan](../../Manuals/FDS_Config_Management_Plan/FDS_Config_Management_Plan.tex) for configuration documentation.
