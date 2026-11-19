#!/usr/bin/env python3
"""Generate the 24 BRE spray validation inputs. Python 3.9+, standard library only.
Source: BRE_Spray_input_creator.xlsm, Worksheet_B and FDS_A_Input.
Edit CONFIG below directly. No Excel or separate configuration file needed.
"""
import argparse
import csv
import math
import re
from pathlib import Path

# USER SETTINGS: edit nozzle speed_F, speed_n, angle_deg and diameter_1bar_um here.
CONFIG = {'project': 'BRE_Spray',
 'water_density_kg_m3': 1000,
 'nozzles': {'A': {'speed_F': 0.3,
                   'speed_n': 0.5,
                   'angle_deg': 50,
                   'diameter_1bar_um': 350},
             'B': {'speed_F': 0.5,
                   'speed_n': 0.6,
                   'angle_deg': 70,
                   'diameter_1bar_um': 500},
             'D': {'speed_F': 0.05,
                   'speed_n': 0.8,
                   'angle_deg': 75,
                   'diameter_1bar_um': 600}},
 'experiments': [{'case': 1,
                  'nozzle': 'A',
                  'pressure_bar': 1,
                  'attenuation_pct': 2.8,
                  'flow_l_min': 0.35,
                  'measured_speed_m_s': 1.2,
                  'measured_diameter_um': 268},
                 {'case': 2,
                  'nozzle': 'A',
                  'pressure_bar': 2,
                  'attenuation_pct': 5.2,
                  'flow_l_min': 0.55,
                  'measured_speed_m_s': 0.9,
                  'measured_diameter_um': 175},
                 {'case': 3,
                  'nozzle': 'A',
                  'pressure_bar': 3,
                  'attenuation_pct': 6.9,
                  'flow_l_min': 0.625,
                  'measured_speed_m_s': 1.7,
                  'measured_diameter_um': 110},
                 {'case': 4,
                  'nozzle': 'A',
                  'pressure_bar': 4,
                  'attenuation_pct': 9.2,
                  'flow_l_min': 0.7,
                  'measured_speed_m_s': 1.5,
                  'measured_diameter_um': 104},
                 {'case': 5,
                  'nozzle': 'A',
                  'pressure_bar': 5,
                  'attenuation_pct': 10.6,
                  'flow_l_min': 0.75,
                  'measured_speed_m_s': 2,
                  'measured_diameter_um': 102},
                 {'case': 6,
                  'nozzle': 'A',
                  'pressure_bar': 6,
                  'attenuation_pct': 12.4,
                  'flow_l_min': 0.875,
                  'measured_speed_m_s': 1.8,
                  'measured_diameter_um': 102},
                 {'case': 7,
                  'nozzle': 'A',
                  'pressure_bar': 7,
                  'attenuation_pct': 11.1,
                  'flow_l_min': 0.95,
                  'measured_speed_m_s': 2.3,
                  'measured_diameter_um': 93},
                 {'case': 8,
                  'nozzle': 'A',
                  'pressure_bar': 8,
                  'attenuation_pct': 13.9,
                  'flow_l_min': 1,
                  'measured_speed_m_s': 1.8,
                  'measured_diameter_um': 126},
                 {'case': 9,
                  'nozzle': 'B',
                  'pressure_bar': 1,
                  'attenuation_pct': 3.6,
                  'flow_l_min': 1.4,
                  'measured_speed_m_s': 1.9,
                  'measured_diameter_um': 392},
                 {'case': 10,
                  'nozzle': 'B',
                  'pressure_bar': 2,
                  'attenuation_pct': 6.3,
                  'flow_l_min': 1.825,
                  'measured_speed_m_s': 2.2,
                  'measured_diameter_um': 266},
                 {'case': 11,
                  'nozzle': 'B',
                  'pressure_bar': 3,
                  'attenuation_pct': 7.7,
                  'flow_l_min': 2,
                  'measured_speed_m_s': 2.9,
                  'measured_diameter_um': 167},
                 {'case': 12,
                  'nozzle': 'B',
                  'pressure_bar': 4,
                  'attenuation_pct': 9.5,
                  'flow_l_min': 2.25,
                  'measured_speed_m_s': 3.5,
                  'measured_diameter_um': 162},
                 {'case': 13,
                  'nozzle': 'B',
                  'pressure_bar': 5,
                  'attenuation_pct': 11.5,
                  'flow_l_min': 2.5,
                  'measured_speed_m_s': 4.6,
                  'measured_diameter_um': 115},
                 {'case': 14,
                  'nozzle': 'B',
                  'pressure_bar': 6,
                  'attenuation_pct': 12.2,
                  'flow_l_min': 2.75,
                  'measured_speed_m_s': 5,
                  'measured_diameter_um': 126},
                 {'case': 15,
                  'nozzle': 'B',
                  'pressure_bar': 7,
                  'attenuation_pct': 13.6,
                  'flow_l_min': 3,
                  'measured_speed_m_s': 5.5,
                  'measured_diameter_um': 156},
                 {'case': 16,
                  'nozzle': 'B',
                  'pressure_bar': 8,
                  'attenuation_pct': 14.9,
                  'flow_l_min': 3.25,
                  'measured_speed_m_s': 3.6,
                  'measured_diameter_um': 148},
                 {'case': 17,
                  'nozzle': 'D',
                  'pressure_bar': 1,
                  'attenuation_pct': 6.1,
                  'flow_l_min': 2.6,
                  'measured_speed_m_s': 1,
                  'measured_diameter_um': 691},
                 {'case': 18,
                  'nozzle': 'D',
                  'pressure_bar': 2,
                  'attenuation_pct': 11,
                  'flow_l_min': 3.75,
                  'measured_speed_m_s': 1.3,
                  'measured_diameter_um': 753},
                 {'case': 19,
                  'nozzle': 'D',
                  'pressure_bar': 3,
                  'attenuation_pct': 16,
                  'flow_l_min': 4.5,
                  'measured_speed_m_s': 1.6,
                  'measured_diameter_um': 794},
                 {'case': 20,
                  'nozzle': 'D',
                  'pressure_bar': 4,
                  'attenuation_pct': 20.5,
                  'flow_l_min': 5,
                  'measured_speed_m_s': 2.4,
                  'measured_diameter_um': 638},
                 {'case': 21,
                  'nozzle': 'D',
                  'pressure_bar': 5,
                  'attenuation_pct': 24.7,
                  'flow_l_min': 5.75,
                  'measured_speed_m_s': 2.9,
                  'measured_diameter_um': 550},
                 {'case': 22,
                  'nozzle': 'D',
                  'pressure_bar': 6,
                  'attenuation_pct': 27.3,
                  'flow_l_min': 6,
                  'measured_speed_m_s': 3.7,
                  'measured_diameter_um': 182},
                 {'case': 23,
                  'nozzle': 'D',
                  'pressure_bar': 7,
                  'attenuation_pct': 35.4,
                  'flow_l_min': 6.75,
                  'measured_speed_m_s': 4.2,
                  'measured_diameter_um': 178},
                 {'case': 24,
                  'nozzle': 'D',
                  'pressure_bar': 8,
                  'attenuation_pct': 35.5,
                  'flow_l_min': 7.5,
                  'measured_speed_m_s': 4.8,
                  'measured_diameter_um': 160}]}

FDS_TEMPLATE = """&HEAD CHID = '{chid}' TITLE='Radiation atten. test (Murrel et al. 1995)' /

! IJK = 80, 40, 60 XB = -2,2,-1,1,0,3 /
&MESH IJK = 40,20,30, XB = -2,2,-1,1,0,3 /


&TIME T_END=15.,DT = 0.05 /
&SPEC ID = 'WATER VAPOR' /
&DUMP DT_PL3D = 30.0 DT_DEVC = 0.1 /

&DEVC  XYZ = 0.0, 0.0,2.24,
PROP_ID =  'Nozzle'
QUANTITY = 'TIME'
SETPOINT = 4.0 /

&PROP  ID = 'Nozzle'
OFFSET  = 0.1
OPERATING_PRESSURE  = {pressure_bar}
K_FACTOR  = {k_factor_l_min_bar_half}
PARTICLES_PER_SECOND = 10000
PART_ID  = 'water drops'
PARTICLE_VELOCITY = {generated_speed_m_s}
SPRAY_ANGLE = 0 , {angle_deg}
SPRAY_PATTERN_SHAPE = 'UNIFORM'
SMOKEVIEW_ID='nozzle' /


&PART  ID='water drops',
SPEC_ID = 'WATER VAPOR'
DIAMETER = {generated_diameter_um} /

&RADI NUMBER_RADIATION_ANGLES = 1000 TIME_STEP_INCREMENT = 2 /


&DEVC XYZ = 2.0,0.00,1.50,  QUANTITY = 'GAUGE HEAT FLUX', ID = 'flux' IOR = -1 /
&DEVC XYZ = 0.0,0.00,1.20, QUANTITY = 'PDPA', PROP_ID = 'pdpa_d43' ID = 'd43' /
&DEVC XYZ = 0.0,0.00,1.20, QUANTITY = 'PDPA', PROP_ID = 'pdpa_w' ID = 'w' /

&PROP ID = 'pdpa_d43'
QUANTITY = 'DIAMETER'
PDPA_START = 2
PDPA_M = 4
PDPA_N = 3
PDPA_RADIUS = 0.05 /

&PROP ID = 'pdpa_w'
PDPA_START = 2
QUANTITY = 'W-VELOCITY'
PDPA_M = 0
PDPA_N = 0
PDPA_RADIUS = 0.05 /


! SLCF PBY=0.0,QUANTITY='TEMPERATURE'/
&SLCF PBY=0.0,QUANTITY='INTEGRATED INTENSITY' /
&SLCF PBY=0.0,QUANTITY='MPUV' PART_ID = 'water drops' /
! SLCF PBY=0.0,QUANTITY='ADD' PART_ID = 'water drops' /
! SLCF PBY=0.0,QUANTITY='WATER VAPOR' /

! Geometry
&VENT XB = -2.0,-2.0,-0.5,0.5,2.0,3.0, SURF_ID = 'OPEN' /
&VENT XB = -2.0,-2.0,-0.5,0.5,0.0,1.0, SURF_ID = 'OPEN'  /
&VENT XB = -2.0,-2.0,-2.0,-0.5,0.0,3.0, SURF_ID = 'OPEN'  /
&VENT XB = -2.0,-2.0, 0.5,2.0,0.0,3.0, SURF_ID = 'OPEN'  /

&VENT XB = 2.0,2.0,-0.5,0.5,2.0,3.0, SURF_ID = 'OPEN' /
&VENT XB = 2.0,2.0,-0.5,0.5,0.0,1.0, SURF_ID = 'OPEN'  /
&VENT XB = 2.0,2.0,-2.0,-0.5,0.0,3.0, SURF_ID = 'OPEN'  /
&VENT XB = 2.0,2.0, 0.5,2.0,0.0,3.0, SURF_ID = 'OPEN'  /
&VENT MB = 'YMIN', SURF_ID = 'OPEN' /
&VENT MB = 'YMAX', SURF_ID = 'OPEN' /
&VENT MB = 'ZMAX', SURF_ID = 'OPEN' /
&VENT MB = 'ZMIN', SURF_ID = 'OPEN' /

&VENT XB = -2.0,-2.0,-0.5,0.5,1.0,2.0, SURF_ID = 'HOT'  /
&SURF ID='HOT',TMP_FRONT = 900.0, COLOR = 'RED' HEAT_TRANSFER_COEFFICIENT = 0. RAMP_T = 'HOT' /
&RAMP ID = 'HOT' T = 0 F = 1. /
&RAMP ID = 'HOT' T = 4.0 F = 1. /
&RAMP ID = 'HOT' T = 4.1 F = 0. /
&RAMP ID = 'HOT' T = 9. F = 0. /
&RAMP ID = 'HOT' T = 9.1 F = 1. /

&VENT XB = 2.0,2.0,-0.5,0.5,1.0,2.0, SURF_ID = 'COOL',  /
&SURF ID='COOL',TMP_FRONT = 20.0, COLOR = 'BLUE' /

&TAIL /
"""

COLUMNS = [
    "case", "nozzle", "pressure_bar", "attenuation_pct", "flow_l_min",
    "measured_speed_m_s", "measured_diameter_um", "k_factor_l_min_bar_half",
    "speed_F", "speed_n", "generated_speed_m_s", "angle_deg",
    "diameter_1bar_um", "generated_diameter_um", "chid", "input_file",
]


def calculate_rows(config):
    """Recompute H, L and O of Worksheet_B, without cached Excel formula values."""
    rho = float(config["water_density_kg_m3"])
    if not math.isfinite(rho) or rho <= 0:
        raise ValueError("Water density must be positive and finite")
    project = config["project"]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", project):
        raise ValueError("Project must contain only letters, digits, underscores or hyphens")
    rows = []
    seen = set()
    for exp in config["experiments"]:
        row = dict(exp)
        nozzle = row["nozzle"]
        p = float(row["pressure_bar"])
        if nozzle not in {"A", "B", "D"} or p not in range(1, 9):
            raise ValueError("Expected nozzle A/B/D and pressure 1 through 8 bar")
        if (nozzle, p) in seen:
            raise ValueError("Duplicate nozzle/pressure combination")
        seen.add((nozzle, p))
        params = config["nozzles"][nozzle]
        row.update(params)
        row["k_factor_l_min_bar_half"] = float(row["flow_l_min"]) / math.sqrt(p)
        row["generated_speed_m_s"] = params["speed_F"] * (2*p*100000/rho)**params["speed_n"]
        row["generated_diameter_um"] = params["diameter_1bar_um"] * (1/p)**(1/3)
        row["chid"] = f"{project}_{nozzle}_{int(p)}"
        row["input_file"] = row["chid"] + ".fds"
        for key in COLUMNS:
            value = row[key]
            if isinstance(value, (int, float)) and not math.isfinite(value):
                raise ValueError(f"Non-finite {key}")
        if row["flow_l_min"] <= 0 or row["generated_speed_m_s"] <= 0 or row["generated_diameter_um"] <= 0:
            raise ValueError("Flow, generated speed and diameter must be positive")
        if not 0 <= row["angle_deg"] <= 180:
            raise ValueError("Angle must be between 0 and 180 degrees")
        rows.append(row)
    if seen != {(n, p) for n in "ABD" for p in range(1, 9)}:
        raise ValueError("Exactly 24 nozzle/pressure combinations are required")
    if [r["case"] for r in rows] != list(range(1, 25)):
        raise ValueError("Case numbers must be 1 through 24 in order")
    return rows


def generate(config, output_dir, template=FDS_TEMPLATE, overwrite=False):
    rows = calculate_rows(config)
    output_dir = Path(output_dir)
    targets = [output_dir / r["input_file"] for r in rows] + [output_dir / "BRE_Spray_summary.csv"]
    if not overwrite:
        existing = [str(p) for p in targets if p.exists()]
        if existing:
            raise FileExistsError("Outputs already exist; use --overwrite: " + ", ".join(existing))
    # Render everything before writing any file. Keep full calculation precision.
    texts = []
    for row in rows:
        values = {k: format(v, ".15g") if isinstance(v, (int, float)) else v for k, v in row.items()}
        texts.append(template.format_map(values))
    output_dir.mkdir(parents=True, exist_ok=True)
    for path, text in zip(targets, texts):
        path.write_text(text, encoding="utf-8")
    # UTF-8 BOM makes the CSV directly readable in Excel; dot decimals, comma delimiter.
    with targets[-1].open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("generated"))
    parser.add_argument("--template", type=Path, help="Optional text template using {column_name} placeholders")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    config = CONFIG
    template = args.template.read_text(encoding="utf-8") if args.template else FDS_TEMPLATE
    try:
        rows = generate(config, args.output_dir, template, args.overwrite)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"Error: {error}\n")
    print(f"Generated {len(rows)} FDS files and BRE_Spray_summary.csv in {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
