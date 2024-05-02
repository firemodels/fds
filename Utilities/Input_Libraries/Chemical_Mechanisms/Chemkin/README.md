## Generating Cantera input files

First, make sure Cantera is installed.  (Activate the `fds_python_env`.)

```
$ cd <path-to-fds-repo>/fds/.github/
$ source fds_python_env/bin/activate
(fds_python_env) $
```

## To generate the Cantera .yaml file, run ck2yaml

Example:

```
(fds_python_env) $ ck2yaml --input=grimech30_kinetics.dat --thermo=grimech30_thermo.dat --transport=grimech30_transport.dat --output=../Cantera/grimech30.yaml
 --output=../Cantera/grimech30.yaml
Wrote YAML mechanism file to '../Cantera/grimech30.yaml'.
Mechanism contains 53 species and 325 reactions.
Validating mechanism...
PASSED
```