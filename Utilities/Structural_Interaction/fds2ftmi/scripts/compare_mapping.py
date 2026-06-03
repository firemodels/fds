import argparse
import re
import os
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Getting SVN number 
GIT="test"

# Defining fonts to be used in plots (Matching generate_plots.py)
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.family':'Times New Roman'})
rcParams.update({'font.size':16})
font1 = {'weight' : 'normal', 'size' : 18}
font2 = {'weight' : 'normal', 'size' : 12}
font3 = {'weight' : 'normal', 'size' : 10}

def parse_ansys_file(file_path):
    data = {}
    times = {} 
    pattern = re.compile(r"\*set,([AH])(\d+)\((\d+),(\d+)\),\s*([-\d.E+e]+)", re.IGNORECASE)

    with open(file_path, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                arr_type = match.group(1).upper()
                node_id = match.group(2)
                step = int(match.group(3))
                index = int(match.group(4))
                
                try:
                    val = float(match.group(5))
                except ValueError:
                    continue

                if index == 0:
                    times[(node_id, step)] = val
                elif index == 1:
                    time = times.get((node_id, step))
                    if time is not None:
                        key = (node_id, time)
                        if key not in data:
                            data[key] = {'T_AST': 0.0, 'h': 0.0}
                        
                        if arr_type == 'A':
                            data[key]['T_AST'] = val
                        elif arr_type == 'H':
                            data[key]['h'] = val
    return data

def compare_mappings(obst_file_path, geom_file_path):
    print("--- Starting Boundary Condition Comparison ---")
    
    try:
        obst_data = parse_ansys_file(obst_file_path)
        print(f"Successfully parsed obst file: {obst_file_path}. Records found: {len(obst_data)}")
    except FileNotFoundError:
        print(f"ERROR: obst file not found at {obst_file_path}. Exiting.")
        return None
    
    try:
        geom_data = parse_ansys_file(geom_file_path)
        print(f"Successfully parsed GEOM file: {geom_file_path}. Records found: {len(geom_data)}")
    except FileNotFoundError:
        print(f"ERROR: GEOM file not found at {geom_file_path}. Exiting.")
        return None

    all_keys = set(obst_data.keys()) | set(geom_data.keys())
    
    total_nodes_compared = 0
    t_ast_leg = []
    t_ast_geom = []
    h_leg = []
    h_geom = []
    
    max_t_ast_diff = 0.0
    max_h_diff = 0.0

    for node_id, time in sorted(list(all_keys)):
        obst_record = obst_data.get((node_id, time))
        geom_record = geom_data.get((node_id, time))
        
        is_missing_in_geom = obst_record and not geom_record
        is_missing_in_obst = geom_record and not obst_record

        if not is_missing_in_geom and not is_missing_in_obst:
            total_nodes_compared += 1
            
            t_ast_leg.append(obst_record['T_AST'])
            t_ast_geom.append(geom_record['T_AST'])
            h_leg.append(obst_record['h'])
            h_geom.append(geom_record['h'])
            
            t_ast_diff = abs(geom_record['T_AST'] - obst_record['T_AST'])
            h_diff = abs(geom_record['h'] - obst_record['h'])
            
            if t_ast_diff > max_t_ast_diff:
                max_t_ast_diff = t_ast_diff
            if h_diff > max_h_diff:
                max_h_diff = h_diff

    print(f"\n✅ Validated {total_nodes_compared} records. Max T_AST diff: {max_t_ast_diff:.4e}, Max h diff: {max_h_diff:.4e}")

    # --- Plotting the Differences Using Parity Plots ---
    os.makedirs('SCRIPT_FIGURES', exist_ok=True)

    if t_ast_leg:
        print("Generating diagnostic parity plots...")
        
        # Plot 1: T_AST Parity
        fig = plt.figure()
        plt.plot(t_ast_leg, t_ast_geom, 'ko', mfc='none', mec='black', label='Mapped Nodes', mew=1.15)
        
        # Add 1:1 Reference Line
        t_min = min(min(t_ast_leg), min(t_ast_geom))
        t_max = max(max(t_ast_leg), max(t_ast_geom))
        pad_t = (t_max - t_min) * 0.05 if t_max != t_min else 1.0
        
        plt.plot([t_min - pad_t, t_max + pad_t], [t_min - pad_t, t_max + pad_t], 'r--', label='1:1 Perfect Match', linewidth=1.5)
        
        plt.xlim([t_min - pad_t, t_max + pad_t])
        plt.ylim([t_min - pad_t, t_max + pad_t])
        plt.xlabel('OBST $T_{AST}$ ($^{\circ}$C)', fontdict=font1)
        plt.ylabel('GEOM $T_{AST}$ ($^{\circ}$C)', fontdict=font1)
        
        plt.legend(numpoints=1, loc=2, prop=font2) # Upper Left
        plt.text(0.05, 0.95, '$T_{AST}$ Mapping Parity (GEOM vs obst)', transform=plt.gca().transAxes, fontsize=14, va='top')
        plt.text(0.98, 1.02, GIT, fontdict=font3, transform=plt.gca().transAxes, ha='right')

        plt.savefig('SCRIPT_FIGURES/t_ast_mapping_parity.pdf', format='pdf')
        plt.close()

        # Plot 2: h Parity
        fig = plt.figure()
        plt.plot(h_leg, h_geom, 'bo', mfc='none', mec='blue', label='Mapped Nodes', mew=1.15)
        
        # Add 1:1 Reference Line
        h_min = min(min(h_leg), min(h_geom))
        h_max = max(max(h_leg), max(h_geom))
        pad_h = (h_max - h_min) * 0.05 if h_max != h_min else 1.0
        
        plt.plot([h_min - pad_h, h_max + pad_h], [h_min - pad_h, h_max + pad_h], 'r--', label='1:1 Perfect Match', linewidth=1.5)
        
        plt.xlim([h_min - pad_h, h_max + pad_h])
        plt.ylim([h_min - pad_h, h_max + pad_h])
        plt.xlabel('OBST $h$ ($W/m^2K$)', fontdict=font1)
        plt.ylabel('GEOM $h$ ($W/m^2K$)', fontdict=font1)
        
        plt.legend(numpoints=1, loc=2, prop=font2) # Upper Left
        plt.text(0.05, 0.95, '$h$ Mapping Parity (GEOM vs OBST)', transform=plt.gca().transAxes, fontsize=14, va='top')
        plt.text(0.98, 1.02, GIT, fontdict=font3, transform=plt.gca().transAxes, ha='right')

        plt.savefig('SCRIPT_FIGURES/h_mapping_parity.pdf', format='pdf')
        plt.close()
        
        print("Successfully saved diagnostic parity plots to SCRIPT_FIGURES/")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compares ANSYS boundary condition files (obst vs. GEOM mapping).")
    parser.add_argument("obst_file", help="Path to the obst boundary condition file.")
    parser.add_argument("geom_file", help="Path to the new GEOM mapping boundary condition file.")

    args = parser.parse_args()
    compare_mappings(args.obst_file, args.geom_file)