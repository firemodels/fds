import argparse
import pandas as pd
import numpy as np
import os
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

def clean_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Strips leading/trailing whitespace from all column names."""
    df.columns = [col.strip() for col in df.columns]
    return df

def plot_thermal_history(legacy_df: pd.DataFrame, geom_df: pd.DataFrame):
    print("Generating Thermal History plot...")
    
    fig = plt.figure()
    
    # OBST (Solid Lines)
    plt.plot(legacy_df['Time']/60, legacy_df['hightemp_A'], 'r-', mfc='none', label='OBST: Hightemp A', linewidth=1.5)
    plt.plot(legacy_df['Time']/60, legacy_df['lowtemp_A'], 'm-', mfc='none', label='OBST: Lowtemp A', linewidth=1.5)
    plt.plot(legacy_df['Time']/60, legacy_df['hightemp_B'], 'b-', mfc='none', label='OBST: Hightemp B', linewidth=1.5)
    plt.plot(legacy_df['Time']/60, legacy_df['lowtemp_B'], 'c-', mfc='none', label='OBST: Lowtemp B', linewidth=1.5)
    
    # GEOM (Dashed Lines)
    plt.plot(geom_df['Time']/60, geom_df['hightemp_A'], 'r--', mfc='none', label='GEOM: Hightemp A', linewidth=1.5)
    plt.plot(geom_df['Time']/60, geom_df['lowtemp_A'], 'm--', mfc='none', label='GEOM: Lowtemp A', linewidth=1.5)
    plt.plot(geom_df['Time']/60, geom_df['hightemp_B'], 'b--', mfc='none', label='GEOM: Hightemp B', linewidth=1.5)
    plt.plot(geom_df['Time']/60, geom_df['lowtemp_B'], 'c--', mfc='none', label='GEOM: Lowtemp B', linewidth=1.5)

    max_time = max(legacy_df['Time']/60)
    
    plt.xlim([0, max_time])
    plt.ylim([0, 700])
    plt.xlabel('Time (min)', fontdict=font1)
    plt.ylabel('Temperature ($^{\circ}$C)', fontdict=font1)
    
    plt.legend(numpoints=1, loc=4, prop=font2) # Lower Right
    
    # Text in plot units to match generate_plots.py style
    plt.text(3, 650, 'Temperature evolution (h_profile vs. h_profile_geom)')
    plt.text(45, 710, GIT, fontdict=font3)
    
    plt.savefig('SCRIPT_FIGURES/h_profile_geom_thermal.pdf', format='pdf')
    plt.close()


def plot_displacement_tracking(legacy_df: pd.DataFrame, geom_df: pd.DataFrame):
    print("Generating Displacement Tracking plots...")
    
    max_time = max(legacy_df['Time']/60)

    # --- X-Displacements ---
    fig = plt.figure()
    
    plt.plot(legacy_df['Time']/60, legacy_df['dx_A']*100, 'r-', mfc='none', label='OBST: dx_A', linewidth=1.5)
    plt.plot(legacy_df['Time']/60, legacy_df['dx_B']*100, 'b-', mfc='none', label='OBST: dx_B', linewidth=1.5)
    plt.plot(legacy_df['Time']/60, legacy_df['dx_C']*100, 'g-', mfc='none', label='OBST: dx_C', linewidth=1.5)
    
    plt.plot(geom_df['Time']/60, geom_df['dx_A']*100, 'r--', mfc='none', label='GEOM: dx_A', linewidth=1.5)
    plt.plot(geom_df['Time']/60, geom_df['dx_B']*100, 'b--', mfc='none', label='GEOM: dx_B', linewidth=1.5)
    plt.plot(geom_df['Time']/60, geom_df['dx_C']*100, 'g--', mfc='none', label='GEOM: dx_C', linewidth=1.5)

    plt.xlim([0, max_time])
    plt.ylim([0, 2])
    plt.xlabel('Time (min)', fontdict=font1)
    plt.ylabel('Displacement (cm)', fontdict=font1)
    
    plt.legend(numpoints=1, loc=4, prop=font2) # Lower Right
    
    # Text string and position in plot units
    plt.text(2, 1.83, 'Displacement in x axis (h_profile vs. h_profile_geom)')
    plt.text(45, 2.02, GIT, fontdict=font3)

    plt.savefig('SCRIPT_FIGURES/h_profile_geom_disp_x.pdf', format='pdf')
    plt.close()

    # --- Y-Displacements ---
    fig = plt.figure()
    
    plt.plot(legacy_df['Time']/60, legacy_df['dy_A']*100, 'r-', mfc='none', label='OBST: dy_A', linewidth=1.5)
    plt.plot(legacy_df['Time']/60, legacy_df['dy_C']*100, 'g-', mfc='none', label='OBST: dy_C', linewidth=1.5)
    
    plt.plot(geom_df['Time']/60, geom_df['dy_A']*100, 'r--', mfc='none', label='GEOM: dy_A', linewidth=1.5)
    plt.plot(geom_df['Time']/60, geom_df['dy_C']*100, 'g--', mfc='none', label='GEOM: dy_C', linewidth=1.5)

    plt.xlim([0, max_time])
    plt.ylim([-0.3, 0.3])
    plt.xlabel('Time (min)', fontdict=font1)
    plt.ylabel('Displacement (cm)', fontdict=font1)
    
    plt.legend(numpoints=1, loc=7, prop=font2) # Center Right
    
    # Text string and position in plot units (-0.25 to place it near the bottom)
    plt.text(2, 0.25, 'Displacement in y axis (h_profile vs. h_profile_geom)')
    plt.text(45, 0.308, GIT, fontdict=font3)

    plt.savefig('SCRIPT_FIGURES/h_profile_geom_disp_y.pdf', format='pdf')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Compares thermal and displacement results between legacy and GEOM models.")
    parser.add_argument("legacy_csv", help="Path to the legacy CSV file.")
    parser.add_argument("geom_csv", help="Path to the GEOM CSV file.")

    args = parser.parse_args()
    
    os.makedirs('SCRIPT_FIGURES', exist_ok=True)

    try:
        legacy_df = clean_column_names(pd.read_csv(args.legacy_csv))
        geom_df = clean_column_names(pd.read_csv(args.geom_csv))

        plot_thermal_history(legacy_df, geom_df)
        plot_displacement_tracking(legacy_df, geom_df)

    except FileNotFoundError as e:
        print(f"ERROR: One or more input CSV files not found. Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred during plotting: {e}")

if __name__ == "__main__":
    main()