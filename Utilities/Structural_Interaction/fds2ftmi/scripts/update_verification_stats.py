import pandas as pd
import argparse
import numpy as np
import math

def calculate_stats(df_expected: pd.DataFrame, df_predicted: pd.DataFrame, metric_name: str, expected_col: str, predicted_col: str, error_type: str, tolerance: float) -> tuple[float, float, str]:
    
    # 1. Find absolute maximum values
    expected_max = df_expected[expected_col].abs().max()
    predicted_max = df_predicted[predicted_col].abs().max()
    
    # 2. Calculate Error
    error = 0.0  
    if error_type.lower() == 'relative':
        if expected_max == 0:
            error = 0.0
        else:
            error = abs(expected_max - predicted_max) / expected_max
    elif error_type.lower() == 'absolute':
        error = abs(expected_max - predicted_max)

    # 3. Evaluation and Formatting
    within_tolerance = "Yes" if error <= tolerance else "No"

    expected_formatted = f"{expected_max:.3e}"
    predicted_formatted = f"{predicted_max:.3e}"
    error_formatted = f"{error:.3e}"

    # Construct the LaTeX row (Added escaping for the underscores in the case name to prevent LaTeX errors)
    latex_row = f"h\\_profile\\_geom & ${expected_formatted}$ & ${predicted_formatted}$ & {error_type} & {error_formatted} & {tolerance} & {within_tolerance} \\\\"
    
    return expected_formatted, predicted_formatted, latex_row

def main():
    """
    Main function to calculate verification statistics and append LaTeX table rows.
    """
    parser = argparse.ArgumentParser(
        description="Calculates verification statistics and appends LaTeX table rows for comparison.",
        epilog="Usage: python update_verification_stats.py <legacy_csv> <geom_csv>"
    )
    parser.add_argument(
        "legacy_csv", 
        help="Path to the legacy CSV file (Expected Metric source)."
    )
    parser.add_argument(
        "geom_csv", 
        help="Path to the GEOM CSV file (Predicted Metric source)."
    )

    args = parser.parse_args()
    
    TOLERANCE = 0.01
    TOLERANCE_T = 0.025

    try:
        # Read data
        df_expected = pd.read_csv(args.legacy_csv)
        df_predicted = pd.read_csv(args.geom_csv)

        print("--- Calculating Verification Statistics ---")
        
        # Metric 1: Max Temperature (Relative Error)
        _, _, latex_temp = calculate_stats(
            df_expected, df_predicted, 
            "Max Temperature", 
            "hightemp_A", "hightemp_A", 
            "Relative", 
            TOLERANCE_T
        )

        # Metric 2: Max X-Displacement (Absolute Error)
        _, _, latex_dx = calculate_stats(
            df_expected, df_predicted, 
            "Max X-Displacement", 
            "dx_A", "dx_A", 
            "Absolute", 
            TOLERANCE
        )

        # Metric 3: Max Y-Displacement (Absolute Error)
        _, _, latex_dy = calculate_stats(
            df_expected, df_predicted, 
            "Max Y-Displacement", 
            "dy_A", "dy_A", 
            "Absolute", 
            TOLERANCE
        )

        # Append to verification_statistics.tex
        output_file = "verification_statistics.tex"
        with open(output_file, "a") as f:
            f.write(latex_temp + "\n")
            f.write(latex_dx + "\n")
            f.write(latex_dy + "\n")
            
        print(f"Successfully appended 3 metrics to {output_file}!")

    except FileNotFoundError as e:
        print(f"\nERROR: One or more input CSV files not found. Please check paths. Details: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred during statistics calculation: {e}")


if __name__ == "__main__":
    main()