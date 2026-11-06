#https://www.thermalradiation.net/sectionb/B-70.html
#note this is for infinite length rectangle so a little error in the exact values
import math
import numpy as np
import pandas as pd

phi = math.pi/4.
cosphi = math.cos(phi)
w = 2.85
TH = 1000
sigma = 5.67E-11
flux = sigma*TH**4

x = np.zeros(40)
x1 = np.zeros(40)
VF = np.zeros(40)
qnet = np.zeros(40)

for i in range(len(x)):
   x[i] = 0.025 + 0.05 * i
   x1[i] = x[i]*math.sqrt(2)
   bigx = x1[i]/w
   VF[i] = 0.5 + (cosphi - bigx)/(2*math.sqrt(1.+bigx**2-2*bigx*cosphi))
   qnet[i] = -flux*(1-VF[i])

df = pd.DataFrame({"x":x,"x1":x1,"VF":VF,"q_net":qnet})

outfile = "net_rad_1000_exact.csv"

# Create a dictionary mapping column names to their formatting styles
formats = {
    'x': '{:.3f}',
    'x1': '{:.3f}',
    'VF': '{:.3f}',
    'q_net': '{:.2f}'
}

# Apply formatting across all designated columns at once
df_export = df.assign(**{col: df[col].map(fmt.format) for col, fmt in formats.items()})

with open(outfile, 'w', encoding='utf-8', newline='') as f:
    # Match the number of commas to your data columns (4 columns = 3 commas)
    f.write("m,m,,kW/m2\n")

    # 3. Append the DataFrame to the same open file
    df_export.to_csv(f, index=False)


