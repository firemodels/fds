import pandas as pd
from pathlib import Path

def generate_csv_from_source(folder_name: str = "../../../out/SNL_Walls/"):
    # 1. Define paths
    folder = Path(folder_name)
    input_path = folder / "concrete_wall_prof_1.csv"
    output_path = folder / "concrete_wall.csv"

    # 2. Open "my_file.csv" skipping the first 3 rows
    # header=None is used to ensure we can access columns by integer index (0, 1, 2...)
    df = pd.read_csv(input_path, skiprows=3, header=None)

    # 3. Extract the data for the new columns
    # Interpretation: (Column Index, Row Range Inclusive)
    # Note: .iloc[start:end] in pandas is exclusive of the end, 
    # so we use +1 to make the range inclusive (e.g., 3:53 becomes 3:54).

    Npoints = df.iloc[0,1]
    
    data = {
        "Depth": df.iloc[ 0, 2:2+Npoints ].values,
        "0.25":  df.iloc[ 9, 2+Npoints:2+2*Npoints].values,
        "0.50":  df.iloc[18, 2+Npoints:2+2*Npoints].values,
        "0.75":  df.iloc[27, 2+Npoints:2+2*Npoints].values,
        "1.00":  df.iloc[36, 2+Npoints:2+2*Npoints].values,
        "1.25":  df.iloc[45, 2+Npoints:2+2*Npoints].values,
        "1.50":  df.iloc[54, 2+Npoints:2+2*Npoints].values
    }

    # 4. Create the new DataFrame with the specific header
    new_df = pd.DataFrame(data)

    # 5. Save to "new_file.csv" in the same directory
    new_df.to_csv(output_path, index=False)
    print(f"Successfully created: {output_path}")

if __name__ == "__main__":
    generate_csv_from_source()

