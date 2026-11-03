# BRE spray input generator

Requires Python 3.9 or newer. Uses only the Python standard library; Excel, macros,
openpyxl and other packages are not needed.

## Run

    python generate_bre_spray.py

Writes 24 `.fds` inputs and `BRE_Spray_summary.csv` to `generated` in the current
working directory. Existing outputs require `--overwrite`:

    python generate_bre_spray.py --overwrite

Edit `CONFIG` near the top of `generate_bre_spray.py` directly. Its `nozzles`
section contains `speed_F`, `speed_n`, `angle_deg`, and `diameter_1bar_um`
for each nozzle A, B and D. The experimental data are in `CONFIG["experiments"]`.
There is no separate configuration file.

An optional external FDS template can be selected with
`--template spray_template.fds.in`; otherwise the embedded template is used.

## Workbook mapping

Source: BRE_Spray_input_creator.xlsm supplied 2 November 2026.
The spray template is FDS_A_Input!A6:I91, rather than FDS_Input, which contains
an unrelated wood-heating model. Experimental and nozzle data are from
Worksheet_B!A7:O30. Names follow FDS_A_Input!B4: BRE_Spray_A_1.fds through
BRE_Spray_D_8.fds. Case order is A 1–8, B 1–8, D 1–8.

The original numerical formulas are recomputed in Python:

- K = flow_l_min / sqrt(pressure_bar)
- speed = F * (2 * pressure_bar * 100000 / water_density_kg_m3)**n
- diameter = diameter_1bar_um * (1 / pressure_bar)**(1/3)

Nozzle parameters:

| Nozzle | Speed F | Speed n | Angle (degrees) | Diameter at 1 bar (um) |
| --- | --- | --- | --- | --- |
| A | 0.3 | 0.5 | 50 | 350 |
| B | 0.5 | 0.6 | 70 | 500 |
| D | 0.05 | 0.8 | 75 | 600 |

Water density is 1000 kg/m3. The empirical exponent n is preserved exactly,
including its application to the complete 2P/rho expression. The measured speed
and diameter are reference data; FDS uses the generated values. K is calculated
individually from each measured flow, rather than using one constant per nozzle.


The CSV follows columns A–O of Worksheet_B, omitting the obsolete column I,
with descriptive headers and two additional fields, chid and input_file.
F and n are repeated per case instead of leaving blanks. The separator is a comma; encoding is UTF-8 with BOM
for Excel. Calculations retain full floating-point precision; FDS values use
15 significant digits. Experimental attenuation is in percent, not fraction.
The workbook does not explicitly label experimental speed/diameter units;
the CSV interprets them as m/s and micrometres, consistent with the generated
speed header and FDS diameter convention.

## Template preservation and validation

The active 40 x 20 x 30 mesh, geometry, timing, radiation, PDPA settings and
other model choices are retained. Inactive alternative mesh/slice lines and
the Geometry heading are prefixed with `!`. The two SPRAY_ANGLE values are
separated by a comma. Optional template edits can be made in the script or
in spray_template.fds.in when passing --template.

All 24 cases were generated. 288 comparisons of the 12 applicable columns
against cached workbook values passed (numeric tolerance 1e-12). Checks also
covered case names, substituted FDS parameters, one active mesh per input,
24 CSV rows, and protection against accidental output overwrites.
FDS was not available in the generation environment, so these inputs have
not been run through FDS; the original model has not been updated for a
particular FDS release. No simulated results or calculated attenuation have
been added to the CSV.
