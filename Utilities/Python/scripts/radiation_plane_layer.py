# radiation_plane_layer
import os
import pandas as pd

dir = '../../Verification/Radiation/'

infile = [['radiation_plane_layer_{}_{}_devc.csv'.format(i, j) for j in range(1, 6)] for i in range(1, 7)]

kappa = ['0.01', '0.1', '0.5', '1.0', '10', r'$\infty$']

exact = [2.8970, 24.9403, 82.9457, 116.2891, 148.9698, 148.9709]

flux = pd.DataFrame(index=range(6), columns=range(5))

for j in range(5):
    for i in range(6):
        filepath = os.path.join(dir, infile[i][j])
        if not os.path.exists(filepath):
            print(f'Error: File {filepath} does not exist. Skipping case.')
            continue
        df = pd.read_csv(filepath, skiprows=3, header=None)
        t = df.iloc[0, 0]
        flux.iloc[i, j] = df.iloc[0, 1]

filename = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/radiation_plane_layer.tex'
with open(filename, 'wt') as fid:
    fid.write('\\begin{center}\n')
    fid.write('\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline\n')
    fid.write('$\\tau$ & $S(\\tau)$ & \\multicolumn{2}{|c|}{FDS (I=20,J=20)} &\n')
    fid.write('\\multicolumn{2}{|c|}{FDS (I=20,J=1)} & FDS (I=150) \\\\ \\cline{3-7}\n')
    fid.write(' & (kW/m$^2$) & 1 band & 6 bands & 1 band & 6 bands & 1 band \\\\ \\hline\\hline\n')
    
    for i in range(6):
        fid.write(f'{kappa[i]} & {exact[i]:9.4f} & ')
        fid.write(' & '.join(f'{f:9.4f}' for f in flux.iloc[i, :5]))
        fid.write(' \\\\\n')
    
    fid.write('\\hline\n')
    fid.write('\\end{tabular}\n')
    fid.write('\\end{center}\n')
