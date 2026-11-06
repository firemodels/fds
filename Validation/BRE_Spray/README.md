# BRE spray input generator

Requires Python 3.9 or newer. 

## Run

    python generate_bre_spray.py

Writes 24 `.fds` inputs and `BRE_Spray_summary.csv` to `generated` in the current
working directory. Existing outputs require `--overwrite`:

    python generate_bre_spray.py --overwrite

Edit `CONFIG` near the top of `generate_bre_spray.py` directly. Its `nozzles`
section contains `speed_F`, `speed_n`, `angle_deg`, and `diameter_1bar_um`
for each nozzle A, B and D. The experimental data are in `CONFIG["experiments"]`.

An optional external FDS template can be selected with
`--template spray_template.fds.in`; otherwise the embedded template is used.

The following formulas are computed 

- K = flow_l_min / sqrt(pressure_bar)
- speed = F * (2 * pressure_bar * 100000 / water_density_kg_m3)**n
- diameter = diameter_1bar_um * (1 / pressure_bar)**(1/3)

Nozzle parameters as per November 1, 2026:

| Nozzle | Speed F | Speed n | Angle (degrees) | Diameter at 1 bar (um) |
| --- | --- | --- | --- | --- |
| A | 0.3 | 0.5 | 50 | 350 |
| B | 0.5 | 0.6 | 70 | 500 |
| D | 0.05 | 0.8 | 75 | 600 |

