import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os
import pandas as pd

# include FDS plot styles, etc.
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)


def reaction_rate(T, Y, dTdt, R0, E, A, residue):
    """
    Reaction rate function with safety checks and normalization of Y.
    """
    if T <= 0:
        raise ValueError("Temperature T must be positive.")

    # Ensure mass fractions are non-negative
    Y = np.maximum(Y, 0.0)

    r = np.zeros(len(Y))
    r[0] = -A[0] * Y[0] * np.exp(-E[0] / (R0 * T)) / dTdt
    r[1] = -A[1] * Y[1] * np.exp(-E[1] / (R0 * T)) / dTdt
    r[2] = -residue[1] * r[1]

    # Normalize Y to ensure the total remains valid
    Y_sum = np.sum(Y)
    if Y_sum > 0:
        Y /= Y_sum

    return r

# Get plot style parameters
plot_style = fdsplotlib.get_plot_style("fds")
plt.rcParams["font.family"] = plot_style["Font_Name"]
plt.rcParams["font.size"] = plot_style["Label_Font_Size"]
# print(plot_style)

# Define the base path for files
base_path = "../../Verification/Pyrolysis/"  # Update this to the relative path of your files

# Define parameters
dTdt = 5.0 / 60.0
R0 = 8314.3
show_fds = True

# Set up plots
for i_plot in range(2):
    fig, ax1 = plt.subplots(figsize=(plot_style["Paper_Width"], plot_style["Paper_Height"]))
    ax2 = ax1.twinx()

    # Define parameters for each plot
    if i_plot == 0:
        n_components = 3
        T_p = np.array([100.0 + 273.0, 300.0 + 273.0])
        Y_0 = np.array([0.00, 1.00, 0.00])
        delta_T = np.array([10.0, -80.0])
        r_p = np.array([0.0, 0.002])
        residue = np.array([0.0, 0.0])
    else:
        n_components = 3
        T_p = np.array([100.0 + 273.0, 300.0 + 273.0])
        Y_0 = np.array([0.1, 0.9, 0.0])
        delta_T = np.array([10.0, 80.0])
        r_p = np.zeros(2)
        residue = np.array([0.0, 0.2])

    # Calculate A and E
    A = np.zeros(n_components - 1)
    E = np.zeros(n_components - 1)
    for i in range(n_components - 1):
        if delta_T[i] > 0:
            r_p[i] = 2 * dTdt * (1.0 - residue[i]) / delta_T[i]
        E[i] = np.exp(1.0) * r_p[i] * R0 * T_p[i] ** 2 / dTdt
        A[i] = np.exp(1.0) * r_p[i] * np.exp(E[i] / (R0 * T_p[i]))

    # Solve ODE
    def odes(T, Y):
        return reaction_rate(T, Y, dTdt, R0, E, A, residue)

    sol = solve_ivp(odes, [273, 673], Y_0, method='BDF', rtol=1e-12, atol=1e-10)

    T = sol.t
    Y = sol.y.T
    m = np.sum(Y, axis=1)

    # Compute dm/dT
    dmdT = np.gradient(m, T)

    # Plot the solution
    color1="tab:blue"
    color2="tab:orange"
    ax1.plot(T - 273, m, label="Normalized Mass", color=color1)
    ax2.plot(T - 273, -dmdT * dTdt * 1000, label="Mass Loss Rate", color=color2)

    # Load FDS solution
    if show_fds:
        if i_plot == 0:
            fds_file = os.path.join(base_path, 'pyrolysis_1_devc.csv')
            if not os.path.exists(fds_file):
                print(f"Error: File {fds_file} does not exist. Skipping case.")
                continue
        else:
            fds_file = os.path.join(base_path, 'pyrolysis_2_devc.csv')
            if not os.path.exists(fds_file):
                print(f"Error: File {fds_file} does not exist. Skipping case.")
                continue

        # Read the FDS data using pandas
        fds_data = pd.read_csv(fds_file, skiprows=2, header=None)

        # Convert to a NumPy array if required for further processing
        fds_data = fds_data.to_numpy()

    # Plot attributes
    ax1.set_xlabel("Temperature (°C)",fontdict={"fontname": plot_style["Font_Name"], "fontsize": plot_style["Label_Font_Size"]})
    ax1.set_ylabel("Normalized Mass",color=color1,fontdict={"fontname": plot_style["Font_Name"], "fontsize": plot_style["Label_Font_Size"]})
    ax2.set_ylabel("Normalized Mass Loss Rate × 1000 (1/s)",color=color2,fontdict={"fontname": plot_style["Font_Name"], "fontsize": plot_style["Label_Font_Size"]})
    ax1.tick_params(axis="y", colors=color1)
    ax2.tick_params(axis="y", colors=color2)
    ax1.set_ylim([0, 1.1])
    ax2.set_ylim([0, 2.2])

    # ax1.legend(loc="upper left")
    # ax2.legend(loc="upper right")

    # Add vertical lines for temperature peaks
    for i in range(n_components - 1):
        ax2.axvline(T_p[i] - 273, color="black", linestyle="-", linewidth=1)

    # Add version sting
    chid = 'pyrolysis_1' if i_plot == 0 else 'pyrolysis_2'
    git_file = os.path.join(base_path, f"{chid}_git.txt")
    fdsplotlib.add_version_string(ax1, git_file, plot_type='linear')

    plt.tight_layout()

    # Save the figure
    plt.savefig(f"../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/pyrolysis_{i_plot + 1}.pdf")
    plt.close()
