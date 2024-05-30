## MATL library

_NOTE: The values provided here are obtained from references or best judgement. They are intended primarily as guidance and for model testing and should not be used blindly._


This directory serves as a library of input text for specific materials - particularly those which include detailed pyrolysis models. These files can be directly appended to FDS inputs using the `&CATF` functionality. This is particularly useful for materials which are repeatedly used in verification and validation cases as changes made to material properties here will be applied to any referencing case. Note that these files primarily contain `&MATL` namelists but may include other inputs which are relevant to the material in question.

**Pine Wood:**

`pine_wood_1C_MATL` uses a single component with two n-th order reactions (pyrolysis and oxidative pyrolysis). Both form char which also undergoes an n-th order oxidation reaction.

| MATL    | Parameter        | Value                    | Ref. | 
|:--------|:-----------------|:-------------------------|-----:|
| `PINE`  | `A(1:2)`         | 4.70E+06, 1.45E+10       | [^1] |
| `PINE`  | `E(1:2)`         | 1.05E+05, 1.27E+05 J/mol | [^1] |
| `PINE`  | `N_S(1:2)`       | 0.87,0.63                | [^1] |
| `PINE`  | `N_O2(1:2)`      | 0.0,0.72                 | [^1] |
| `PINE`  | `N_MATL(1,1:2)`  | 0.31,0.31                | [^1] |
| `CHAR`  | `A`              | 7.55E+07                 | [^1] |
| `CHAR`  | `E`              | 1.24E+05 J/mol           | [^1] |
| `CHAR`  | `N_S`            | 0.56                     | [^1] |
| `CHAR`  | `N_O2`           | 0.68                     | [^1] |

`pine_wood_3C_MATL` uses a three components, each with two n-th order reactions (pyrolysis and oxidative pyrolysis). All initial reactions form char which also undergoes an n-th order oxidation reaction.

| MATL      | Parameter        | Value                    | Ref. | 
|:----------|:-----------------|:-------------------------|-----:|
| `PINE 1`  | `A(1:2)`         | 4.9E+9,8.9E+9            | [^1] |
| `PINE 1`  | `E(1:2)`         | 1.46E+05, 1.66E+05 J/mol | [^1] |
| `PINE 1`  | `N_S(1:2)`       | 0.56,0.3                 | [^1] |
| `PINE 1`  | `N_O2(1:2)`      | 0.0,0.61                 | [^1] |
| `PINE 1`  | `N_MATL(1,1:2)`  | 0.25,0.25                | [^1] |
| `PINE 2`  | `A(1:2)`         | 5.0E+10,2.0E+4           | [^1] |
| `PINE 2`  | `E(1:2)`         | 1.44E+05, 7.5E+04  J/mol | [^1] |
| `PINE 2`  | `N_S(1:2)`       | 1.0,1.0                  | [^1] |
| `PINE 2`  | `N_O2(1:2)`      | 0.0,0.49                 | [^1] |
| `PINE 2`  | `N_MATL(1,1:2)`  | 0.25,0.25                | [^1] |
| `PINE 3`  | `A(1:2)`         | 2.9E+11,2.6              | [^1] |
| `PINE 3`  | `E(1:2)`         | 1.64E+05, 1.64E+05 J/mol | [^1] |
| `PINE 3`  | `N_S(1:2)`       | 1.25,5.67                | [^1] |
| `PINE 3`  | `N_O2(1:2)`      | 0.0,0.66                 | [^1] |
| `PINE 3`  | `N_MATL(1,1:2)`  | 0.25,0.25                | [^1] |
| `CHAR`    | `A`              | 7.55E+07                 | [^1] |
| `CHAR`    | `E`              | 1.24E+05 J/mol           | [^1] |
| `CHAR`    | `N_S`            | 0.56                     | [^1] |
| `CHAR`    | `N_O2`           | 0.68                     | [^1] |


Following [^1], the the mass fractions of `PINE 1`, `PINE 2`, and `PINE 3` are 0.55, 0.10, and 0.35, respectively.


### Vegetation

Vegetative materials are based on the general set of inputs specified in `generic_vegetation_1C_MATL.fds`. References for the values used by this material file can be found in the FDS user guide. Unique vegetation types are listed below, along with the parameters which have been modified from the generic model and specific references.

**Little Bluestem Grass:**

`little_bluestem_1C_MATL.fds` uses a single component with a first-order reaction to model pyrolysis.

| MATL              | Parameter   | Value          | Ref. | 
|:------------------|:------------|:---------------|-----:|
| `LITTLE_BLUESTEM` | `A`         | 2.455E+04 1/s  | [^2] |
| `LITTLE_BLUESTEM` | `E`         | 5.82E+04 J/mol | [^2] |
| `LITTLE_BLUESTEM` | `NU_MATL`   | 0.23           | [^3] |
| `LITTLE_BLUESTEM` | `NU_SPEC`   | 0.77           | [^3] |

**Douglas Fir:**

`douglas_fir_3C_MATL.fds` uses three components with first-order reactions to model pyrolysis.

| MATL          | Parameter   | Value          | Ref. | 
|:--------------|:------------|:---------------|-----:|
| `COMPONENT 1` | `A`         | 3.35E+04 1/s   | [^4] |
| `COMPONENT 1` | `E`         | 7.07E+04 J/mol | [^4] |
| `COMPONENT 1` | `NU_MATL`   | 0.23           | [^4] |
| `COMPONENT 1` | `NU_SPEC`   | 0.77           | [^4] |
| `COMPONENT 2` | `A`         | 1.45E+07 1/s   | [^4] |
| `COMPONENT 2` | `E`         | 1.09E+05 J/mol | [^4] |
| `COMPONENT 2` | `NU_MATL`   | 0.23           | [^4] |
| `COMPONENT 2` | `NU_SPEC`   | 0.77           | [^4] |
| `COMPONENT 3` | `A`         | 2.13E+01 1/s   | [^4] |
| `COMPONENT 3` | `E`         | 5.17E+04 J/mol | [^4] |
| `COMPONENT 3` | `NU_MATL`   | 0.23           | [^4] |
| `COMPONENT 3` | `NU_SPEC`   | 0.77           | [^4] |


[^1]:[Anca-Couce et al., 2012](https://doi:10.1016/j.combustflame.2011.11.015)
[^2]:[Amini et al., 2021](https://doi.org/10.1016/j.jaap.2021.105167)
[^3]:[Amini et al., 2019](https://doi.org/10.1016/j.fuel.2018.08.112)
[^4]:[Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762)
