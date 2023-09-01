## MATL library

This directory serves as a library of input text for specific materials - particularly those which include detailed pyrolysis models. These files can be directly appended to FDS inputs using the `&CATF` functionality. This is particularly useful for materials which are repeatedly used in verification and validation cases as changes made to material properties here will be applied to any referencing case. Note that these files primarily contain `&MATL` namelists but may include other inputs which are relevant to the material in question.

### Vegetation

Vegetative materials are based on the general set of inputs specified in `generic_vegetation_1R_MATL.fds`. References for the values used by this material file can be found in the FDS user guide. Unique vegetation types are listed below, along with the parameters which have been modified from the generic model and specific references.

**Little Bluestem Grass:**

`little_bluestem_1R_MATL.fds` uses a single 1-step reaction to model pyrolysis.

|Namelist | ID                | Parameter                 | Value          | Reference                 |
|:--------|:------------------|:--------------------------|:---------------|:-------------------------:|
| `MATL`  | `LITTLE BLUESTEM` | `A`                       | 2.455E+04 1/s  | [Amini et al., 2021](https://doi.org/10.1016/j.jaap.2021.105167) |
| `MATL`  | `LITTLE BLUESTEM` | `E`                       | 5.82E+04 J/mol | [Amini et al., 2021](https://doi.org/10.1016/j.jaap.2021.105167) |
| `MATL`  | `LITTLE BLUESTEM` | `NU_MATL`, `NU_SPEC`      | 0.23, 0.77     | [Amini et al., 2019](https://doi.org/10.1016/j.fuel.2018.08.112) |


**Douglas Fir:**

`douglas_fir_3R_MATL.fds` uses three 1-step reactions to model pyrolysis.

|Namelist | ID                | Parameter                 | Value          | Reference                 |
|:--------|:------------------|:--------------------------|:---------------|:-------------------------:|
| `MATL`  | `COMPONENT 1`     | `A`                       | 3.35E+04 1/s   | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |
| `MATL`  | `COMPONENT 1`     | `E`                       | 7.07E+04 J/mol | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |
| `MATL`  | `COMPONENT 1`     | `NU_MATL`, `NU_SPEC`      | 0.23, 0.77     | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |
| `MATL`  | `COMPONENT 2`     | `A`                       | 1.45E+07 1/s   | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |
| `MATL`  | `COMPONENT 2`     | `E`                       | 1.09E+05 J/mol | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |
| `MATL`  | `COMPONENT 2`     | `NU_MATL`, `NU_SPEC`      | 0.23, 0.77     | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |
| `MATL`  | `COMPONENT 3`     | `A`                       | 2.13E+01 1/s   | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |
| `MATL`  | `COMPONENT 3`     | `E`                       | 5.17E+04 J/mol | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |
| `MATL`  | `COMPONENT 3`     | `NU_MATL`, `NU_SPEC`      | 0.23, 0.77     | [Leventon et al., 2022](https://doi.org/10.1016/j.firesaf.2023.103762) |