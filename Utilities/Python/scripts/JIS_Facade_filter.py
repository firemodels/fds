
import pandas as pd

expdir = '../../../exp/Submodules/macfp-db/Wall_Fires/JIS_Facade/Experimental_Data/'

df = pd.read_csv(expdir + 'Sun_FAM_2024_mean_temperature.csv')

# Remove rows where the first column does not have a zero
df_filtered = df[df.iloc[:, 0] == 0]

# Save to a new CSV file with header
df_filtered.to_csv(expdir + 'Sun_FAM_2024_mean_temperature_surface.csv', index=False)

