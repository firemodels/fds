#!/usr/bin/env python3
"""
Extract RAW (pre-masking) Measured and Predicted values from an
fdsplotlib.dataplot pickle for a given Save_Quantity.

This reproduces the data that scatplot concatenates into
Measured_Values and Predicted_Values BEFORE any filtering.

Output columns:
    csv_rownum, Dataname, Predicted_Values, Measured_Values

Usage:
    python pkl2csv.py saved_data_validation.pkl "HGL Temperature; Natural Ventilation" raw_values.csv
"""

import argparse
import pickle
import csv
import sys
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract raw scatplot Measured/Predicted values "
            "for a given Save_Quantity from a dataplot pickle."
        )
    )
    parser.add_argument(
        "pkl_file",
        help="Pickle file produced by fdsplotlib.dataplot",
    )
    parser.add_argument(
        "quantity",
        help="Exact Scatter_Plot_Title / Save_Quantity to extract",
    )
    parser.add_argument(
        "output_csv",
        nargs="?",
        default="raw_scatplot_values.csv",
        help="Output CSV filename (default: raw_scatplot_values.csv)",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Load pickle
    # ------------------------------------------------------------
    try:
        with open(args.pkl_file, "rb") as f:
            saved_data, drange = pickle.load(f)
    except Exception as e:
        print(f"Error reading pickle file: {e}", file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------
    # Unpack saved_data (canonical dataplot order)
    # ------------------------------------------------------------
    (
        Save_Quantity,
        Save_Group_Style,
        Save_Fill_Color,
        Save_Group_Key_Label,
        Save_Measured_Metric,
        Save_Predicted_Metric,
        Save_Dataname,
        Save_Plot_Filename,
        Save_Dep_Title,
        Save_Error_Tolerance,
        Save_Metric_Type,
        Save_Measured_Quantity,
        Save_Predicted_Quantity,
        Save_csv_rownum,
    ) = saved_data

    # ------------------------------------------------------------
    # Match dataplot entries EXACTLY as scatplot does
    # ------------------------------------------------------------
    target = args.quantity.strip().lower()

    match_idx = [
        i for i, q in enumerate(Save_Quantity)
        if str(q).strip().lower() == target
    ]

    if not match_idx:
        print(
            f"No dataplot entries found for Save_Quantity = '{args.quantity}'",
            file=sys.stderr,
        )
        sys.exit(1)

    # ------------------------------------------------------------
    # Write CSV
    # ------------------------------------------------------------
    with open(args.output_csv, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "csv_rownum",
            "Dataname",
            "Predicted_Values",
            "Measured_Values",
        ])

        for idx in match_idx:
            rownum = Save_csv_rownum[idx]
            dataname = Save_Dataname[idx]

            # RAW values: exactly what scatplot starts from
            mvals = np.array(Save_Measured_Metric[idx], dtype=float).flatten()
            pvals = np.array(Save_Predicted_Metric[idx], dtype=float).flatten()

            n = min(len(mvals), len(pvals))

            for k in range(n):
                writer.writerow([
                    rownum,
                    dataname,
                    pvals[k],
                    mvals[k],
                ])

    print(
        f"Wrote {args.output_csv} "
        f"({len(match_idx)} dataplot entries, raw values)"
    )


if __name__ == "__main__":
    main()
