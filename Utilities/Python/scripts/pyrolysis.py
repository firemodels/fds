import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
import sys

# Add path to Verification/Pyrolysis
sys.path.append('../../Verification/Pyrolysis')

# Define global variables
dTdt = None
R0 = None
E = None
A = None
residue = None

def plot_style():
    # Define plot style settings (placeholder)
    plt.style.use('seaborn-darkgrid')
    global Plot_Units, Plot_X, Plot_Y, Plot_Width, Plot_Height
    global Font_Name, Label_Font_Size, Font_Interpreter
    global Figure_Visibility, Paper_Units, Paper_Width, Paper_Height
    Plot_Units = 'normalized'
    Plot_X = 0.1
    Plot_Y = 0.1
    Plot_Width = 0.8
    Plot_Height = 0.8
    Font_Name = 'Arial'
    Label_Font_Size = 12
    Font_Interpreter = 'latex'
    Figure_Visibility = 'on'
    Paper_Units = 'inches'
    Paper_Width = 6
    Paper_Height = 4

def addverstr(ax, filename, mode):
    # Placeholder for addverstr function
    # In Python, you might read the version string from the file and add text to the plot
    if os.path.exists(filename):
        with open(filename, 'r') as file:
            ver_str = file.read().strip()
            ax.text(0.05, 0.95, ver_str, transform=ax.transAxes, fontsize=8,
                    verticalalignment='top')

def reaction_rate(T, Y):
    """
    Defines the ODE system dY/dt = -A * Y * exp(-E / (R * T))
    """
    global A, E, R0
    return -A * Y * np.exp(-E / (R0 * T))

def main():
    global dTdt, R0, E, A, residue

    plt.close('all')
    plt.figure()

    plot_style()

    show_fds = True

    for i_plot in range(1, 3):
        fig, ax1 = plt.subplots()
        fig.set_size_inches(Paper_Width, Paper_Height)
        ax1.set_position([Plot_X, Plot_Y, Plot_Width, Plot_Height])

        dTdt = 5.0 / 60.0  # degrees C per second
        R0 = 8314.3  # J/(kmol*K)

        if i_plot == 1:
            n_components = 3
            T_p = np.array([100.0 + 273.0, 300.0 + 273.0])
            Y_0 = np.array([0.00, 1.00, 0.00])
            delta_T = np.array([10.0, -80.0])
            r_p = np.zeros(n_components - 1)
            residue = np.array([0.0, 0.0])
        else:
            n_components = 3
            T_p = np.array([100.0 + 273.0, 300.0 + 273.0])
            Y_0 = np.array([0.1, 0.9, 0.0])
            delta_T = np.array([10.0, 80.0])
            residue = np.array([0.0, 0.2])
            r_p = np.zeros(n_components - 1)

        # Calculate A and E
        E = np.zeros(n_components - 1)
        A = np.zeros(n_components - 1)
        for i in range(n_components - 1):
            if delta_T[i] > 0:
                r_p[i] = 2.0 * dTdt * (1.0 - residue[i]) / delta_T[i]
            E[i] = np.exp(1.0) * r_p[i] * R0 * T_p[i] ** 2 / dTdt
            A[i] = np.exp(1.0) * r_p[i] * np.exp(E[i] / (R0 * T_p[i]))

        # Solve the ODE dY/dT = f(T, Y)
        # Since dY/dt = f(T, Y) and dT/dt is constant, dY/dT = f(T, Y) / dTdt
        def dydT(T, Y):
            return reaction_rate(T, Y) / dTdt

        sol = solve_ivp(dydT, [273, 673], Y_0, method='BDF',
                       rtol=1e-10, atol=1e-8)
        T = sol.t
        Y = sol.y.T

        # Compute dm/dT
        m = np.sum(Y, axis=1)
        dmdT = np.gradient(m, T)
        dmdT[0] = 0.0
        dmdT[-1] = 0.0

        # Plot the solution
        TC = T - 273.0
        ax2 = ax1.twinx()
        line1, = ax1.plot(TC, m, label='Normalized Mass')
        line2, = ax2.plot(TC, -dmdT * dTdt * 1000.0, label='Normalized Mass Loss Rate × 10³ (s⁻¹)', color='orange')
        ax1.set_xlabel('Temperature (°C)', fontsize=Label_Font_Size, fontname=Font_Name, 
                       fontfamily=Font_Interpreter)
        ax1.set_ylabel('Normalized Mass', fontsize=Label_Font_Size, fontname=Font_Name, 
                       fontfamily=Font_Interpreter)
        ax2.set_ylabel('Normalized Mass Loss Rate × 10³ (s⁻¹)', fontsize=Label_Font_Size, fontname=Font_Name,
                       fontfamily='LaTeX')
        ax1.set_ylim([0, 1.1])
        ax2.set_ylim([0, 2.2])
        ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax2.set_yticks([0, 0.4, 0.8, 1.2, 1.6, 2.0])
        line1.set_linestyle('-')
        line2.set_linestyle('-')

        # Add the FDS solution
        if show_fds:
            filename = f'pyrolysis_{i_plot}_devc.csv'
            if not os.path.exists(filename):
                print(f'Error: File {filename} does not exist. Skipping case.')
                continue
            FDS = np.genfromtxt(filename, delimiter=',', skip_header=2)
            ax1.plot(FDS[:, 3], FDS[:, 1], 'b--', label='FDS Mass')
            ax2.plot(FDS[:, 3], FDS[:, 2] * 500.0, 'r--', label='FDS Mass Loss Rate × 500')

        # Add vertical lines to indicate the temperature peaks
        for i in range(n_components - 1):
            ax2.axvline(x=T_p[i] - 273.0, color='black', linewidth=1)

        # Add version string if file is available
        chid = f'pyrolysis_{i_plot}'
        git_filename = f'{chid}_git.txt'
        addverstr(ax1, git_filename, 'linear')

        # Create the PDF files
        fig.set_visible(True if Figure_Visibility == 'on' else False)
        fig.set_size_inches(Paper_Width * 1.1, Paper_Height)
        pdf_filename = f'../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/{chid}.pdf'
        plt.savefig(pdf_filename, format='pdf')

    plt.show()

if __name__ == "__main__":
    main()

