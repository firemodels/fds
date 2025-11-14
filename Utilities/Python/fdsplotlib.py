#!/usr/bin/env python
"""
fdsplotlib.py
by Randy McDermott
January 2025

Fire Dynamics Simulator (FDS) Plot Library

Collection of functions for plotting and analysis
"""

import os
import sys
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
# rc('text', usetex=True) # Enable TeX rendering
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd


def expand_ranges(items, df, header_rows=1):
    """
    Expand a list of row specifiers into pandas row indices.
    Supports:
      - int: row number (1-based, with headers counted)
      - "start:stop": inclusive range of row numbers
      - "start:": open-ended range (to the end)
      - "all": keep everything
      - plain string: match Dataname column (case-insensitive)
    """
    nrows = len(df)
    result = []

    for item in items:
        if isinstance(item, int):
            # single row number
            result.append(item - (header_rows + 1))

        elif isinstance(item, str):
            s = item.strip()
            low = s.lower()

            if low == "all":
                return list(range(nrows))  # everything

            if ":" in s:  # range form
                start, _, end = s.partition(":")
                start = int(start)

                if end == "":
                    # open-ended: go to the last row of df
                    end = nrows + header_rows  # Excel-style row count
                else:
                    end = int(end)

                # Convert to iloc positions (0-based, exclusive of end)
                start_pos = start - (header_rows + 1)
                end_pos = end - header_rows

                # Clamp so we don't go past the last row
                end_pos = min(end_pos, nrows)

                rng = range(start_pos, end_pos)
                result.extend(rng)

            else:
                # assume it's a Dataname match
                matches = df.index[df['Dataname'].str.lower() == low].tolist()
                if not matches:
                    raise ValueError(f"No match for Dataname '{item}'")
                result.extend(matches)

        else:
            raise TypeError(f"Unsupported plot_range element: {item}")

    return sorted(set(result))


def _compute_metrics_block(
    x, Y, metric, initial_value,
    comp_start, comp_end, dep_comp_start, dep_comp_end,
    variant_side="d1",
):
    """
    MATLAB dataplot.m metric logic, shape-safe Python equivalent.

    Returns:
      vals_flat : 1D np.array (metrics for each curve, or concatenated for 'all')
      titles    : list of metric labels
      per_curve_series : list of per-curve metric arrays (for Metric='all')
    """
    import numpy as np

    # --- normalize inputs ---
    x = np.asarray(x).reshape(-1)
    Y = np.asarray(Y)
    if Y.ndim == 1:
        Y = Y.reshape(-1, 1)
    elif Y.ndim > 2:
        Y = np.squeeze(Y)
        if Y.ndim == 1:
            Y = Y.reshape(-1, 1)

    N, ncols = Y.shape

    # --- comparison mask (like MATLAB) ---
    comp_mask = np.isfinite(x)
    if np.isfinite(comp_start):
        comp_mask &= (x >= comp_start)
    if np.isfinite(comp_end):
        comp_mask &= (x <= comp_end)

    if np.isfinite(dep_comp_start) or np.isfinite(dep_comp_end):
        y0 = Y[:, 0]
        dep_mask = np.isfinite(y0)
        if np.isfinite(dep_comp_start):
            dep_mask &= (y0 >= dep_comp_start)
        if np.isfinite(dep_comp_end):
            dep_mask &= (y0 <= dep_comp_end)
        comp_mask &= dep_mask

    if not np.any(comp_mask):
        return np.array([]), [], []

    x_sel = x[comp_mask]
    Y_sel = Y[comp_mask, :]

    # --- support patterns like mean_1_2, max_2_1, end_1_2 ---
    def _parse_stat_xy(m):
        for base in ("max", "mean", "end"):
            pref = base + "_"
            if m.startswith(pref):
                try:
                    a, b = m[len(pref):].split("_", 1)
                    return base, int(a), int(b)
                except Exception:
                    pass
        return m, None, None

    metric_str = str(metric).strip().lower()
    base, idx_first, idx_second = _parse_stat_xy(metric_str)

    vals = []
    titles = []
    per_curve_series = []

    # --- stat_x_y: only compute for the specified first index (MATLAB uses 1-based) ---
    if idx_first is not None:
        j = idx_first - 1
        if j < 0 or j >= ncols:
            return np.array([]), [], []
        yj = Y_sel[:, j].reshape(-1)

        if base == "max":
            out = np.nanmax(yj) - initial_value
        elif base == "mean":
            out = abs(np.nanmean(yj) - initial_value)
        elif base == "end":
            out = yj[-1] - initial_value
        else:
            out = np.nan

        if out == 0.0:
            out = 1e-12

        return np.array([out]), [f"curve{idx_first}"], []

    # --- metric='all': return all finite Y values (one per data point) ---
    if metric_str == "all":
        for j in range(ncols):
            yj = Y_sel[:, j].reshape(-1)
            mask = np.isfinite(yj)
            yj = yj[mask] - initial_value
            per_curve_series.append(yj)
            titles.extend([f"point{k+1}_curve{j+1}" for k in range(len(yj))])
        vals_flat = np.concatenate(per_curve_series) if per_curve_series else np.array([])
        return vals_flat, titles, per_curve_series

    # --- scalar per-curve metrics ---
    for j in range(ncols):
        yj = Y_sel[:, j].reshape(-1)
        if metric_str == "max":
            out = np.nanmax(yj) - initial_value
        elif metric_str == "min":
            out = initial_value - np.nanmin(yj)
        elif metric_str == "maxabs":
            out = np.nanmax(np.abs(yj - initial_value))
        elif metric_str == "slope":
            msk = np.isfinite(x_sel) & np.isfinite(yj)
            out = np.polyfit(x_sel[msk], yj[msk], 1)[0] if msk.sum() >= 2 else np.nan
        elif metric_str == "mean":
            out = abs(np.nanmean(yj) - initial_value)
        elif metric_str == "threshold":
            out = np.nanmin(yj) - initial_value
        elif metric_str == "tolerance":
            out = np.nanmax(np.abs(yj - initial_value))
        elif metric_str == "area":
            out = np.trapz(yj, x_sel) - initial_value
        elif metric_str == "end":
            out = yj[-1] - initial_value
        elif metric_str == "start":
            out = yj[0]
        elif metric_str == "ipct":
            out = 1e-12  # placeholder for parity
        else:
            out = 1e-12

        if out == 0.0:
            out = 1e-12

        vals.append(out)
        titles.append(f"curve{j+1}")

    return np.asarray(vals, dtype=float).reshape(-1), titles, []


def dataplot(config_filename, **kwargs):

    import os
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    from matplotlib import rc
    import logging

    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    _csv_cache = {}
    def read_csv_cached(path, **kwargs):
        if path not in _csv_cache:
            _csv_cache[path] = pd.read_csv(path, **kwargs)
        return _csv_cache[path].copy()

    configdir = kwargs.get('configdir', '')
    revision = kwargs.get('revision', '')
    expdir = kwargs.get('expdir', '')
    cmpdir = kwargs.get('cmpdir', '')
    pltdir = kwargs.get('pltdir', '')
    close_figs = kwargs.get('close_figs', False)
    verbose = kwargs.get('verbose', False)

    plot_list = kwargs.get('plot_list', ['all'])
    plot_range_in = kwargs.get('plot_range', None)
    header_rows = kwargs.get('header_rows', 1)

    drange = []
    Save_Quantity = []
    Save_Group_Style = []
    Save_Fill_Color = []
    Save_Group_Key_Label = []
    Save_Measured_Metric = []
    Save_Predicted_Metric = []
    Save_Dataname = []
    Save_Plot_Filename = []
    Save_Dep_Title = []
    Save_Error_Tolerance = []
    Save_Metric_Type = []
    Save_Measured_Quantity = []
    Save_Predicted_Quantity = []

    default_na = {
        '', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN',
        '1.#IND', '1.#QNAN', '<NA>', 'NA', 'NULL', 'NaN',
        'n/a', 'nan', 'null'
    }
    safe_na_values = default_na

    df = pd.read_csv(
        configdir + config_filename,
        sep=',',
        engine='python',
        quotechar='"',
        na_values=safe_na_values,
        keep_default_na=False
    )
    df["__orig_index__"] = df.index
    C = df.where(pd.notnull(df), None)

    quantity_filter = kwargs.get("quantity_filter", None)
    if quantity_filter:
        if isinstance(quantity_filter, str):
            quantity_filter = [quantity_filter]
        mask = False
        for q in quantity_filter:
            mask = mask | C["Quantity"].str.contains(q, case=False, na=False)
        C = C[mask]
        if C.empty:
            raise ValueError(f"[dataplot] No rows match Quantity filter(s): {quantity_filter}")
        if verbose:
            matched_rows = (C["__orig_index__"] + 2).tolist()
            print(f"[dataplot] Filtering by Quantity={quantity_filter}")

    if plot_range_in is not None:
        if isinstance(plot_range_in, (list, tuple)):
            adjusted = expand_ranges(plot_range_in, C, header_rows)
            if len(adjusted) < len(C):
                C = C.iloc[adjusted]
        elif isinstance(plot_range_in, range):
            start, end = plot_range_in.start, plot_range_in.stop
            adjusted = range(start - (header_rows + 1), end - header_rows)
            C = C.iloc[adjusted]
        else:
            raise TypeError("plot_range must be a list, tuple, or range")
    elif plot_list and 'all' not in [p.lower() for p in plot_list]:
        C = C[C['Dataname'].str.lower().isin([p.lower() for p in plot_list])]

    otest_active = any(str(C.iloc[j]['switch_id']).strip().lower() == 'o' for j in range(len(C)))
    f_Last = plt.figure()
    col_orig_idx = C.columns.get_loc("__orig_index__")
    col_idx = {col: i for i, col in enumerate(C.columns)}

    fast_mode = bool(kwargs.get('fast_mode', True))
    if verbose:
        print(f"[dataplot] {'Running in fast_mode' if fast_mode else 'Running in full mode'}")

    for pos, row in enumerate(C.itertuples(index=False, name=None)):

        csv_rownum = int(row[col_orig_idx]) + header_rows + 1
        pp = define_plot_parameters(C, pos, lightweight=fast_mode)

        if pp.switch_id == 's':
            continue
        if otest_active and pp.switch_id != 'o':
            continue

        gtest = (pp.switch_id == 'g')
        dtest = (pp.switch_id == 'd')
        ftest = (pp.switch_id == 'f')

        if not (dtest or ftest or gtest or pp.switch_id == 'o'):
            if verbose:
                print(f"[dataplot] Skipping unrecognized switch_id '{pp.switch_id}' on line {csv_rownum}")
            continue

        if not gtest:
            drange.append(int(row[col_orig_idx]) + 2)
            Save_Dataname.append(pp.Dataname)
            Save_Plot_Filename.append(pp.Plot_Filename)
            Save_Dep_Title.append(pp.Dep_Title)
            Save_Error_Tolerance.append(pp.Error_Tolerance)
            Save_Metric_Type.append(pp.Metric)
            Save_Group_Key_Label.append(pp.Group_Key_Label)
            Save_Quantity.append(pp.Quantity)
            Save_Group_Style.append(pp.Group_Style)
            Save_Fill_Color.append(pp.Fill_Color)
            Save_Measured_Metric.append(np.nan)
            Save_Predicted_Metric.append(np.nan)
            Save_Measured_Quantity.append(None)
            Save_Predicted_Quantity.append(None)

        # ---------------------- LOAD EXP ----------------------
        E = read_csv_cached(expdir + pp.d1_Filename,
                            header=int(pp.d1_Col_Name_Row - 1),
                            sep=',', engine='python', quotechar='"',
                            skip_blank_lines=True).dropna(how='all')
        E = E.loc[:E.dropna(how='all').last_valid_index()]
        E.columns = E.columns.str.strip()
        start_idx = int(pp.d1_Data_Row - pp.d1_Col_Name_Row - 1)
        x, _ = get_data(E, pp.d1_Ind_Col_Name, start_idx)
        y, _ = get_data(E, pp.d1_Dep_Col_Name, start_idx)

        flip_axis = str(pp.Flip_Axis).strip().lower() in ['yes', 'true', '1']
        x_scale = float(pp.Scale_Ind or 1.0)
        y_scale = float(pp.Scale_Dep or 1.0)
        x_scaled = np.asarray(x) / x_scale
        y_scaled = np.asarray(y) / y_scale

        if x_scaled.ndim == 2 and y_scaled.ndim == 2 and x_scaled.shape[1] == y_scaled.shape[1]:
            x_plot_list = [x_scaled[:, i] for i in range(x_scaled.shape[1])]
            y_plot_list = [y_scaled[:, i] for i in range(y_scaled.shape[1])]
        elif y_scaled.ndim == 2 and y_scaled.shape[1] > 1:
            x_plot_list = [x_scaled for _ in range(y_scaled.shape[1])]
            y_plot_list = [y_scaled[:, i] for i in range(y_scaled.shape[1])]
        else:
            x_plot_list = [np.ravel(x_scaled)]
            y_plot_list = [np.ravel(y_scaled)]

        for xi, yi in zip(x_plot_list, y_plot_list):
            if len(xi) != len(yi):
                print(f"[dataplot] Pair length mismatch in {pp.Dataname}: x={len(xi)}, y={len(yi)}")

        plot_type = str(pp.Plot_Type or '').strip().lower()
        if plot_type not in ['linear', 'loglog', 'semilogx', 'semilogy']:
            plot_type = 'linear'

        raw_styles = [c.strip() for c in (pp.d1_Style or '').split('|')] if pp.d1_Style else []
        styles = (raw_styles + [None] * len(y_plot_list))[:len(y_plot_list)]
        raw_keys = [c.strip() for c in (pp.d1_Key or '').split('|')] if pp.d1_Key else []
        key_labels = (raw_keys + [None] * len(y_plot_list))[:len(y_plot_list)]

        if dtest or gtest:
            if verbose:
                print(f"Generating plot {csv_rownum} {pltdir}{pp.Plot_Filename}...")
            if close_figs:
                plt.close('all')
            first_plot = True
        elif ftest:
            f = f_Last
            first_plot = False
        else:
            continue

        # --- Styles and keys (EXP d1) ---
        d1_raw_styles = [c.strip() for c in (pp.d1_Style or '').split('|')] if pp.d1_Style else []
        d1_styles = (d1_raw_styles + [None] * len(y_plot_list))[:len(y_plot_list)]
        d1_raw_keys = [c.strip() for c in (pp.d1_Key or '').split('|')] if pp.d1_Key else []
        d1_key_labels = (d1_raw_keys + [None] * len(y_plot_list))[:len(y_plot_list)]

        # --- Plot Exp curves ---
        for i, (x_i, y_i) in enumerate(zip(x_plot_list, y_plot_list)):
            f = plot_to_fig(
                x_data=y_i if flip_axis else x_i,
                y_data=x_i if flip_axis else y_i,
                figure_handle=None if (first_plot and i == 0) else f,
                data_label=d1_key_labels[i],
                x_label=pp.Dep_Title if flip_axis else pp.Ind_Title,
                y_label=pp.Ind_Title if flip_axis else pp.Dep_Title,
                marker_style=d1_styles[i],
                x_min=pp.Min_Dep if flip_axis else pp.Min_Ind,
                x_max=pp.Max_Dep if flip_axis else pp.Max_Ind,
                y_min=pp.Min_Ind if flip_axis else pp.Min_Dep,
                y_max=pp.Max_Ind if flip_axis else pp.Max_Dep,
                legend_location=matlab_legend_to_matplotlib(pp.Key_Position),
                plot_type=plot_type,
            )

        # --- Save measured (experimental) ---
        if not gtest:
            try:
                vals_meas_list = []
                qty_meas_list = []
                if y.ndim == 2 and x.ndim == 2 and y.shape[1] == x.shape[1]:
                    for j in range(y.shape[1]):
                        xj = np.ravel(x[:, j])
                        yj = np.ravel(y[:, j])
                        mask = np.isfinite(xj) & np.isfinite(yj)
                        xj, yj = xj[mask], yj[mask]
                        if len(xj) > 0 and len(yj) > 0:
                            vals_meas, qty_meas, _ = _compute_metrics_block(
                                x=xj, Y=yj, metric=pp.Metric,
                                initial_value=float(pp.d1_Initial_Value or 0.0),
                                comp_start=float(pp.d1_Comp_Start or np.nan),
                                comp_end=float(pp.d1_Comp_End or np.nan),
                                dep_comp_start=float(pp.d1_Dep_Comp_Start or np.nan),
                                dep_comp_end=float(pp.d1_Dep_Comp_End or np.nan),
                                variant_side="d1",
                            )
                            vals_meas_list.append(vals_meas)
                            qty_meas_list.append(qty_meas)
                else:
                    vals_meas, qty_meas, _ = _compute_metrics_block(
                        x=x, Y=y, metric=pp.Metric,
                        initial_value=float(pp.d1_Initial_Value or 0.0),
                        comp_start=float(pp.d1_Comp_Start or np.nan),
                        comp_end=float(pp.d1_Comp_End or np.nan),
                        dep_comp_start=float(pp.d1_Dep_Comp_Start or np.nan),
                        dep_comp_end=float(pp.d1_Dep_Comp_End or np.nan),
                        variant_side="d1",
                    )
                    vals_meas_list = [vals_meas]
                    qty_meas_list = [qty_meas]

                # Replace placeholder "curve#" labels with actual dependent column name
                qty_label = str(pp.d1_Dep_Col_Name).strip() or "Unknown"

                Save_Measured_Metric[-1] = np.array(vals_meas_list, dtype=object)
                Save_Measured_Quantity[-1] = np.array([qty_label] * len(vals_meas_list), dtype=object)

            except Exception as e:
                print(f"[dataplot] Error computing measured metric for {pp.Dataname}: {e}")
                Save_Measured_Metric[-1] = np.array([])
                Save_Measured_Quantity[-1] = []

        # ---------------------- LOAD MODEL ----------------------
        M = read_csv_cached(cmpdir + pp.d2_Filename,
                            header=int(pp.d2_Col_Name_Row - 1),
                            sep=',', engine='python', quotechar='"',
                            skip_blank_lines=True).dropna(how='all')
        M = M.loc[:M.dropna(how='all').last_valid_index()]
        M.columns = M.columns.str.strip()
        start_idx = int(pp.d2_Data_Row - pp.d2_Col_Name_Row - 1)

        version_string = revision
        if pp.VerStr_Filename:
            try:
                with open(cmpdir + pp.VerStr_Filename, "r") as fver:
                    Lines = fver.readlines()
                    if Lines:
                        version_string = Lines[0].strip()
            except Exception as e:
                print(f"[dataplot] Warning: could not read version string: {e}")

        x, _ = get_data(M, pp.d2_Ind_Col_Name, start_idx)
        y, _ = get_data(M, pp.d2_Dep_Col_Name, start_idx)

        x_scale = float(pp.Scale_Ind or 1.0)
        y_scale = float(pp.Scale_Dep or 1.0)
        x_scaled = np.asarray(x) / x_scale
        y_scaled = np.asarray(y) / y_scale

        if x_scaled.ndim == 2 and y_scaled.ndim == 2 and x_scaled.shape[1] == y_scaled.shape[1]:
            x_plot_list = [x_scaled[:, i] for i in range(x_scaled.shape[1])]
            y_plot_list = [y_scaled[:, i] for i in range(y_scaled.shape[1])]
        elif y_scaled.ndim == 2 and y_scaled.shape[1] > 1:
            x_plot_list = [x_scaled for _ in range(y_scaled.shape[1])]
            y_plot_list = [y_scaled[:, i] for i in range(y_scaled.shape[1])]
        else:
            x_plot_list = [np.ravel(x_scaled)]
            y_plot_list = [np.ravel(y_scaled)]

        for xi, yi in zip(x_plot_list, y_plot_list):
            if len(xi) != len(yi):
                print(f"[dataplot] Pair length mismatch in {pp.Dataname}: x={len(xi)}, y={len(yi)}")

        # --- Styles and keys (MODEL d2) ---
        d2_raw_styles = [c.strip() for c in (pp.d2_Style or '').split('|')] if pp.d2_Style else []
        d2_styles = (d2_raw_styles + [None] * len(y_plot_list))[:len(y_plot_list)]
        d2_raw_keys = [c.strip() for c in (pp.d2_Key or '').split('|')] if pp.d2_Key else []
        d2_key_labels = (d2_raw_keys + [None] * len(y_plot_list))[:len(y_plot_list)]

        # --- Plot model curves ---
        for i, (x_i, y_i) in enumerate(zip(x_plot_list, y_plot_list)):
            f = plot_to_fig(
                x_data=y_i if flip_axis else x_i,
                y_data=x_i if flip_axis else y_i,
                revision_label=version_string,
                figure_handle=f,                  # keep same figure
                data_label=d2_key_labels[i],
                line_style=d2_styles[i],
                x_label=pp.Dep_Title if flip_axis else pp.Ind_Title,
                y_label=pp.Ind_Title if flip_axis else pp.Dep_Title,
                x_min=pp.Min_Dep if flip_axis else pp.Min_Ind,
                x_max=pp.Max_Dep if flip_axis else pp.Max_Ind,
                y_min=pp.Min_Ind if flip_axis else pp.Min_Dep,
                y_max=pp.Max_Ind if flip_axis else pp.Max_Dep,
                legend_location=matlab_legend_to_matplotlib(pp.Key_Position),
                plot_type=plot_type,
                plot_title=pp.Plot_Title,
            )

        # --- Interpolated, metric-aware model logic ---
        if not gtest:
            try:
                metric_str = str(pp.Metric or '').strip().lower()
                meas_list, pred_list, qty_pred_list = [], [], []

                # Load experimental again for alignment (safe; cached)
                E = read_csv_cached(expdir + pp.d1_Filename,
                                    header=int(pp.d1_Col_Name_Row - 1),
                                    sep=',', engine='python', quotechar='"',
                                    skip_blank_lines=True).dropna(how='all')
                E.columns = E.columns.str.strip()
                start_idx_exp = int(pp.d1_Data_Row - pp.d1_Col_Name_Row - 1)
                x_exp, _ = get_data(E, pp.d1_Ind_Col_Name, start_idx_exp)
                y_exp, _ = get_data(E, pp.d1_Dep_Col_Name, start_idx_exp)

                # Normalize shapes to 2D (col-major semantics)
                x_exp = np.atleast_2d(x_exp)
                y_exp = np.atleast_2d(y_exp)
                x_mod = np.atleast_2d(x)
                y_mod = np.atleast_2d(y)

                ncols = min(y_exp.shape[1], y_mod.shape[1])

                for j in range(ncols):
                    # Cull NaNs per series keeping pairs aligned
                    xj_e = np.ravel(x_exp[:, j] if x_exp.shape[1] > 1 else x_exp)
                    yj_e = np.ravel(y_exp[:, j])
                    m_e = np.isfinite(xj_e) & np.isfinite(yj_e)
                    xj_e, yj_e = xj_e[m_e], yj_e[m_e]

                    xj_m = np.ravel(x_mod[:, j] if x_mod.shape[1] > 1 else x_mod)
                    yj_m = np.ravel(y_mod[:, j])
                    m_m = np.isfinite(xj_m) & np.isfinite(yj_m)
                    xj_m, yj_m = xj_m[m_m], yj_m[m_m]

                    if metric_str == 'all':
                        # align by interpolating model to exp x
                        if xj_m.size < 2 or xj_e.size == 0:
                            continue
                        yj_m_i = np.interp(xj_e, xj_m, yj_m, left=np.nan, right=np.nan)
                        mask_pair = np.isfinite(yj_m_i) & np.isfinite(yj_e)
                        if not np.any(mask_pair):
                            continue
                        x_use = xj_e[mask_pair]
                        y_exp_use = yj_e[mask_pair]
                        y_mod_use = yj_m_i[mask_pair]
                        # compute both on the same x grid
                        v_meas, _, _ = _compute_metrics_block(
                            x=x_use, Y=y_exp_use, metric=metric_str,
                            initial_value=float(pp.d1_Initial_Value or 0.0),
                            comp_start=float(pp.d1_Comp_Start or np.nan),
                            comp_end=float(pp.d1_Comp_End or np.nan),
                            dep_comp_start=float(pp.d1_Dep_Comp_Start or np.nan),
                            dep_comp_end=float(pp.d1_Dep_Comp_End or np.nan),
                            variant_side="d1",
                        )
                        v_pred, qty_pred, _ = _compute_metrics_block(
                            x=x_use, Y=y_mod_use, metric=metric_str,
                            initial_value=float(pp.d2_Initial_Value or 0.0),
                            comp_start=float(pp.d2_Comp_Start or np.nan),
                            comp_end=float(pp.d2_Comp_End or np.nan),
                            dep_comp_start=float(pp.d2_Dep_Comp_Start or np.nan),
                            dep_comp_end=float(pp.d2_Dep_Comp_End or np.nan),
                            variant_side="d2",
                        )
                    else:
                        # aggregate metrics: NO interpolation; operate independently
                        if yj_e.size == 0 or yj_m.size == 0:
                            continue
                        v_meas, _, _ = _compute_metrics_block(
                            x=xj_e, Y=yj_e, metric=metric_str,
                            initial_value=float(pp.d1_Initial_Value or 0.0),
                            comp_start=float(pp.d1_Comp_Start or np.nan),
                            comp_end=float(pp.d1_Comp_End or np.nan),
                            dep_comp_start=float(pp.d1_Dep_Comp_Start or np.nan),
                            dep_comp_end=float(pp.d1_Dep_Comp_End or np.nan),
                            variant_side="d1",
                        )
                        v_pred, qty_pred, _ = _compute_metrics_block(
                            x=xj_m, Y=yj_m, metric=metric_str,
                            initial_value=float(pp.d2_Initial_Value or 0.0),
                            comp_start=float(pp.d2_Comp_Start or np.nan),
                            comp_end=float(pp.d2_Comp_End or np.nan),
                            dep_comp_start=float(pp.d2_Dep_Comp_Start or np.nan),
                            dep_comp_end=float(pp.d2_Dep_Comp_End or np.nan),
                            variant_side="d2",
                        )

                    meas_list.append(np.atleast_1d(v_meas))
                    pred_list.append(np.atleast_1d(v_pred))
                    qty_pred_list.append(qty_pred)

                flat_meas = np.concatenate(meas_list) if meas_list else np.array([])
                flat_pred = np.concatenate(pred_list) if pred_list else np.array([])

                nmin = min(flat_meas.size, flat_pred.size)
                if nmin == 0:
                    print(f"[dataplot] Warning: no valid data pairs for {pp.Dataname}")
                else:
                    if flat_meas.size != flat_pred.size:
                        print(f"[dataplot] Truncated unequal vectors for {pp.Dataname}: "
                              f"Measured={flat_meas.size}, Predicted={flat_pred.size} → {nmin}")
                        # Truncate both sides to maintain one-to-one correspondence
                        flat_meas = flat_meas[:nmin]
                        flat_pred = flat_pred[:nmin]

                    # Save truncated paired arrays (but don’t overwrite earlier measured data)
                    Save_Measured_Metric[-1] = flat_meas
                    Save_Predicted_Metric[-1] = flat_pred

                    qty_label = str(pp.d2_Dep_Col_Name).strip() or "Unknown"
                    Save_Predicted_Quantity[-1] = np.array([qty_label] * len(flat_pred), dtype=object)

            except Exception as e:
                print(f"[dataplot] Error computing predicted metric for {pp.Dataname}: {e}")
                Save_Predicted_Metric[-1] = np.array([])
                Save_Predicted_Quantity[-1] = []

        plt.figure(f.number)
        os.makedirs(pltdir, exist_ok=True)
        plt.savefig(pltdir + pp.Plot_Filename + '.pdf', backend='pdf')
        f_Last = f

    saved_data = [
        Save_Quantity, Save_Group_Style, Save_Fill_Color, Save_Group_Key_Label,
        Save_Measured_Metric, Save_Predicted_Metric, Save_Dataname, Save_Plot_Filename,
        Save_Dep_Title, Save_Error_Tolerance, Save_Metric_Type,
        Save_Measured_Quantity, Save_Predicted_Quantity,
    ]

    for i, (m, p, name, qty) in enumerate(zip(
        Save_Measured_Metric, Save_Predicted_Metric, Save_Dataname, Save_Quantity
    )):
        len_m = np.size(m) if isinstance(m, np.ndarray) else 0
        len_p = np.size(p) if isinstance(p, np.ndarray) else 0
        csv_rownum = drange[i] if i < len(drange) else "?"
        if len_m != len_p:
            print(f"[dataplot] Length mismatch at CSV row {csv_rownum}: {name} | {qty} | Measured={len_m}, Predicted={len_p}")

    print("[dataplot] returning saved_data and drange")
    return saved_data, drange


def get_data(E, spec, start_idx):
    """
    Extract data columns from DataFrame E according to spec string.

    spec: "colA" | "colA|colB" | "colA+colB"
    start_idx: integer row index where numeric data starts

    Returns
    -------
    tuple (y, col_names)
        y : 2D numpy array (nrows x ncols)
        col_names : list of str, resolved column names (with '+' grouped)
    """
    names = [s.strip() for s in spec.split('|')]
    out = []
    col_names = []
    for name in names:
        if '+' in name:
            cols = [n.strip() for n in name.split('+')]
            series = E[cols].iloc[start_idx:].astype(float).sum(axis=1).values
            out.append(series)
            col_names.append('+'.join(cols))  # keep group name
        else:
            series = E[[name]].iloc[start_idx:].astype(float).values.ravel()
            out.append(series)
            col_names.append(name)
    y = np.column_stack(out) if len(out) > 1 else np.array(out[0]).reshape(-1, 1)
    return y, col_names


def plot_to_fig(x_data,y_data,**kwargs):
    """
    Create a simple x,y plot and return the fig handle
    """
    # # useful debug statements
    # print(x_data)
    # print(y_data)
    # for key, value in kwargs.items():
    #     print ("%s == %s" %(key, value))

    plot_style = get_plot_style("fds")

    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "pdf.use14corefonts": True,
        "text.usetex": False,

        # Text and math in Times New Roman
        "font.family": "serif",
        "font.serif": ["Times", "Times New Roman"],

        "mathtext.fontset": "custom",
        "mathtext.rm": "Times",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times:bold",
        "mathtext.cal": "Times New Roman:italic",
        "mathtext.tt": "Courier New",
        "mathtext.default": "it",

        "axes.unicode_minus": False,
        "pdf.compression": 9,
    })

    import logging
    # Suppress just the 'findfont' warnings from matplotlib's font manager
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    ##### default parameters ######
    default_figure_size = (plot_style["Paper_Width"],plot_style["Paper_Height"])
    default_plot_size = (plot_style["Plot_Width"],plot_style["Plot_Height"])
    default_plot_origin = (plot_style["Plot_X"],plot_style["Plot_Y"])
    default_ticklabel_fontsize = plot_style["Label_Font_Size"]
    default_axeslabel_fontsize = plot_style["Label_Font_Size"]
    default_legend_fontsize = plot_style["Key_Font_Size"]
    default_title_fontsize = plot_style["Title_Font_Size"]
    default_version_fontsize = 10
    default_legend_location = 'best'
    default_legend_framealpha = 1
    default_markevery = 1
    markerfacecolor = None
    markeredgecolor = 'black'
    marker = None
    linestyle = '-'
    color = 'black'
    ###############################

    figure_size=kwargs.get('figure_size',default_figure_size)
    plot_size=kwargs.get('plot_size',default_plot_size)
    plot_origin=kwargs.get('plot_origin',default_plot_origin)
    version_fontsize=kwargs.get('version_fontsize',default_version_fontsize)

    # if figure handle is passed, append to current figure, else generate a new figure
    if kwargs.get('figure_handle'):
        fig = kwargs.get('figure_handle')
        if fig.axes:
            ax = fig.axes[0]
        else:
            # Ensure we have at least one axes
            ax = fig.add_subplot(111)
        plt.figure(fig.number)
        using_existing_figure = True
    else:
        fig = plt.figure(figsize=figure_size)
        using_existing_figure = False
        # Convert to fractions of the figure size:
        ax_w = plot_size[0] / figure_size[0]
        ax_h = plot_size[1] / figure_size[1]
        left   = plot_origin[0] / figure_size[0]
        bottom = plot_origin[1] / figure_size[1]
        ax = fig.add_axes([left, bottom, ax_w, ax_h])

    # widen paper if legend is outside, keeping axes fixed in physical size
    if kwargs.get('legend_location') == 'outside':
        legend_expand = kwargs.get('legend_expand', 1.25)
        old_size = fig.get_size_inches()

        # Compute current axes rectangle in *inches*
        bbox = ax.get_position()
        ax_left_in = bbox.x0 * old_size[0]
        ax_bottom_in = bbox.y0 * old_size[1]
        ax_w_in = bbox.width * old_size[0]
        ax_h_in = bbox.height * old_size[1]

        # Widen the figure canvas
        new_width = old_size[0] * legend_expand
        fig.set_size_inches(new_width, old_size[1])

        # Recompute normalized coordinates to keep the axes fixed in size and location
        left_new = ax_left_in / new_width
        bottom_new = ax_bottom_in / old_size[1]
        ax_w_new = ax_w_in / new_width
        ax_h_new = ax_h_in / old_size[1]

        ax.set_position([left_new, bottom_new, ax_w_new, ax_h_new])


    # select plot type
    plot_type=kwargs.get('plot_type','linear')

    # convert matlab styles to matplotlib
    style = kwargs.get('marker_style','ko')
    color,marker,linestyle = parse_matlab_style(style)

    if kwargs.get('line_style'):
        style = kwargs.get('line_style')
        color,marker,linestyle = parse_matlab_style(style)

    marker_fill_color = kwargs.get('marker_fill_color',None)
    markerfacecolor = marker_fill_color

    error_fill_color = kwargs.get('error_fill_color',None)

    # adjust sizes if requested
    linewidth = kwargs.get('linewidth',1)
    markeredgewidth = kwargs.get('markeredgewidth',1)
    markersize = kwargs.get('markersize',5)

    # adjust ticks if required
    xnumticks = kwargs.get('xnumticks',None)
    ynumticks = kwargs.get('ynumticks',None)
    xticks = kwargs.get('xticks',None)
    yticks = kwargs.get('yticks',None)

    # other plot parameters
    markevery = kwargs.get('data_markevery',default_markevery)
    legend_location = kwargs.get('legend_location',default_legend_location)
    legend_framealpha = kwargs.get('legend_framealpha',default_legend_framealpha)

    # set dashes to default, or user requested
    # This set is the matplotlib default
    if linestyle == '': dashes = (None, None); linewidth = 0;
    if linestyle == '-': dashes = (None, None)
    if linestyle == '--': dashes = kwargs.get('line_dashes',(6, 6))
    if linestyle == '-.': dashes = kwargs.get('line_dashes',(6, 3, 1, 3))
    if linestyle == ':': dashes = kwargs.get('line_dashes',(1, 3))

    # This set is what we were using in Matlab
    # if linestyle == '': dashes = (None, None); linewidth = 0;
    # if linestyle == '-': dashes = (None, None)
    # if linestyle == '--': dashes = kwargs.get('line_dashes',(10, 6.2))
    # if linestyle == '-.': dashes = kwargs.get('line_dashes',(12, 7.4, 3, 7.4))
    # if linestyle == ':': dashes = kwargs.get('line_dashes',(1, 3))


    data_label = kwargs.get('data_label',None)

    # trap any data_labels set to blank (old matlab convention)
    if isinstance(data_label, str) and data_label.lower() == 'blank':
        data_label = None

    # generate the main x,y plot
    if plot_type=='linear':
        ax.plot(x_data,y_data,
            markevery=markevery,
            label=data_label,
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
            dashes=dashes)

    if plot_type=='loglog':
        ax.loglog(x_data,y_data,
            markevery=markevery,
            label=data_label,
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
            dashes=dashes)

    if plot_type=='semilogx':
        ax.semilogx(x_data,y_data,
            markevery=markevery,
            label=data_label,
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
            dashes=dashes)

    if plot_type=='semilogy':
        ax.semilogy(x_data,y_data,
            markevery=markevery,
            label=data_label,
            markerfacecolor=markerfacecolor,
            markeredgecolor=color,
            markeredgewidth=markeredgewidth,
            marker=marker,
            markersize=markersize,
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
            dashes=dashes)

    # if y fill range is passed, add it to the plot
    if kwargs.get('y_error_fill_absolute') and not kwargs.get('y_error_fill_relative'):
        if kwargs.get('y_error_fill_absolute')>0.:
            ax.fill_between(x_data,y_data-kwargs.get('y_fill_absolute'),y_data+kwargs.get('y_fill_absolute'),
                alpha=0.1,color=error_fill_color)

    if kwargs.get('y_error_fill_relative') and not kwargs.get('y_error_fill_absolute'):
        if kwargs.get('y_error_fill_relative')>0.:
            ax.fill_between(x_data,y_data*(1.-kwargs.get('y_error_fill_relative')),y_data*(1.+kwargs.get('y_error_fill_relative')),
                alpha=0.1,color=error_fill_color)

    if kwargs.get('y_error_fill_relative') and kwargs.get('y_error_fill_absolute'):
        if kwargs.get('y_error_fill_relative')>0.:
            ax.fill_between(x_data,y_data*(1.-kwargs.get('y_error_fill_relative'))-kwargs.get('y_error_fill_absolute'),
                                   y_data*(1.+kwargs.get('y_error_fill_relative'))+kwargs.get('y_error_fill_absolute'),
                alpha=0.1,color=error_fill_color)

    if kwargs.get('y_error_fill'):
        y_error_fill = kwargs.get('y_error_fill')
        if len(y_data)==len(y_error_fill):
            ax.fill_between(x_data,y_data-y_error_fill,y_data+y_error_fill,
                alpha=0.1,color=error_fill_color)
        else:
            raise ValueError(f"y_fill must the same length as y_data")

    xerr = kwargs.get('x_error', None)
    yerr = kwargs.get('y_error', None)
    err_linewidth=kwargs.get('error_linewidth', 1)
    if xerr is not None or yerr is not None:
        ax.errorbar(
            x_data, y_data,
            xerr=xerr,                               # can be scalar, array, or [lower, upper]
            yerr=yerr,                               # same flexibility
            fmt=style,                               # marker style for data points
            markeredgewidth=markeredgewidth,         # marker edge width
            markerfacecolor=markerfacecolor,         # make marker hollow
            markeredgecolor=color,                   # outline color
            elinewidth=err_linewidth,
            capsize=kwargs.get('error_capsize', 5),  # size of caps at ends
            capthick=linewidth,
        )


    ticklabel_fontsize=kwargs.get('ticklabel_fontsize',default_ticklabel_fontsize)
    plt.setp( ax.xaxis.get_majorticklabels(), rotation=0, fontsize=ticklabel_fontsize )
    plt.setp( ax.yaxis.get_majorticklabels(), rotation=0, fontsize=ticklabel_fontsize )

    axeslabel_fontsize=kwargs.get('axeslabel_fontsize',default_axeslabel_fontsize)
    if not using_existing_figure:
        plt.xlabel(kwargs.get('x_label'), fontsize=axeslabel_fontsize)
        plt.ylabel(kwargs.get('y_label'), fontsize=axeslabel_fontsize)

    legend_fontsize=kwargs.get('legend_fontsize',default_legend_fontsize)

    # --- always get current handles and labels ---
    handles, labels = ax.get_legend_handles_labels()

    # --- case 1 or 2: creating a new figure ---
    if not using_existing_figure:
        # record legend display properties on the Axes
        if legend_location == 'outside':
            ax._legend_loc = 'center left'
            ax._legend_bbox = (1.02, 0.5)
            ax.figure.subplots_adjust(right=0.8)
        else:
            ax._legend_loc = legend_location
            ax._legend_bbox = None

        ax._legend_fontsize = legend_fontsize
        ax._legend_framealpha = legend_framealpha

        # if we already have labeled data, draw legend now
        if labels:
            ax.legend(handles, labels,
                      loc=ax._legend_loc,
                      bbox_to_anchor=ax._legend_bbox,
                      fontsize=ax._legend_fontsize,
                      framealpha=ax._legend_framealpha)
        else:
            # optionally show empty placeholder (for consistent layout)
            ax.legend([], [],
                      loc=ax._legend_loc,
                      bbox_to_anchor=ax._legend_bbox,
                      fontsize=ax._legend_fontsize,
                      frameon=False)

    # --- case 3: existing figure, adding more data ---
    else:
        loc = getattr(ax, '_legend_loc', 'best')
        bbox = getattr(ax, '_legend_bbox', None)
        fontsize = getattr(ax, '_legend_fontsize', legend_fontsize)
        framealpha = getattr(ax, '_legend_framealpha', legend_framealpha)

        if labels:
            ax.legend(handles, labels,
                      loc=loc,
                      bbox_to_anchor=bbox,
                      fontsize=fontsize,
                      framealpha=framealpha)


    # plot title
    if kwargs.get('plot_title'):
        if kwargs.get('title_fontsize'):
            title_fontsize=kwargs.get('title_fontsize')
        else:
            title_fontsize=default_title_fontsize

        plt.text(0.05, 0.95, kwargs.get('plot_title'),
        transform=plt.gca().transAxes,
        fontsize=title_fontsize,
        verticalalignment='top',
        horizontalalignment='left')

    # set axes and tick properties
    xmin=kwargs.get('x_min')
    xmax=kwargs.get('x_max')
    ymin=kwargs.get('y_min')
    ymax=kwargs.get('y_max')

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    # set number of ticks if requested by the user
    if xnumticks != None:
        if plot_type in ('loglog', 'semilogx'):
            ax.set_xticks(np.logspace(xmin, xmax, xnumticks))
        else:
            ax.set_xticks(np.linspace(xmin, xmax, xnumticks))
    if ynumticks != None:
        if plot_type in ('loglog', 'semilogy'):
            ax.set_yticks(np.logspace(ymin, ymax, ynumticks))
        else:
            ax.set_yticks(np.linspace(ymin, ymax, ynumticks))

    # set ticks if requested by the user
    if xticks is not None: ax.set_xticks(xticks)
    if yticks is not None: ax.set_yticks(yticks)


    if plot_type in ('loglog', 'semilogy'):
        apply_global_exponent(ax, axis='y', fontsize=axeslabel_fontsize)
    if plot_type in ('loglog', 'semilogx'):
        apply_global_exponent(ax, axis='x', fontsize=axeslabel_fontsize)

    if kwargs.get('revision_label'):
        add_version_string(ax=ax, version_str=kwargs.get('revision_label'), plot_type=plot_type, font_size=version_fontsize)

    # fig.tight_layout() # this should not be needed if figure_size and plot_size are both specified

    set_ticks_like_matlab(fig)

    return fig


def apply_global_exponent(ax, axis='y', fontsize=10, minor_subs=None):
    import numpy as np
    from matplotlib.ticker import FuncFormatter, LogLocator

    if axis == 'y':
        ticks = ax.get_yticks()
        axis_obj = ax.yaxis
        lims = ax.get_ylim()
    else:
        ticks = ax.get_xticks()
        axis_obj = ax.xaxis
        lims = ax.get_xlim()

    # Keep only positive finite ticks (for log axes)
    ticks = np.array([t for t in ticks if t > 0 and np.isfinite(t)])
    if ticks.size == 0:
        return

    # Choose representative exponent
    exp = int(np.floor(np.log10(np.median(ticks))))
    scale = 10.0**exp

    # Major formatter: fixed-point decimals
    def fmt(val, pos):
        v = val / scale
        return "{:g}".format(v)

    axis_obj.set_major_formatter(FuncFormatter(fmt))
    axis_obj.get_offset_text().set_visible(False)

    # Decide what to do with minor tick labels
    span_decades = np.log10(lims[1]) - np.log10(max(lims[0], 1e-300))

    # Default subs = 2, 4, 6, 8
    if minor_subs is None:
        minor_subs = [2, 4, 6, 8]

    if span_decades <= 1.1:  # only ~1 decade
        axis_obj.set_minor_locator(LogLocator(base=10.0, subs=minor_subs, numticks=10))
        axis_obj.set_minor_formatter(FuncFormatter(lambda val, pos: "{:g}".format(val/scale)))
    else:
        ax.tick_params(axis=axis, which='minor', labelleft=False, labelbottom=False)

    # Force tick labels NOT to go through TeX
    if axis == 'y':
        for label in ax.get_yticklabels() + ax.get_yticklabels(minor=True):
            label.set_usetex(False)
    else:
        for label in ax.get_xticklabels() + ax.get_xticklabels(minor=True):
            label.set_usetex(False)

    # Place the ×10^exp text at the axis end
    if exp != 0:
        if axis == 'y':
            ax.text(0, 1.01, rf"$\times 10^{{{exp}}}$",
                    transform=ax.transAxes,
                    ha='left', va='bottom', fontsize=fontsize)
        else:
            ax.text(1.0, -0.1, rf"$\times 10^{{{exp}}}$",
                    transform=ax.transAxes,
                    ha='right', va='top', fontsize=fontsize)



def parse_matlab_style(style):
    color = ''
    marker = ''
    linestyle = ''

    # Check for the color (first character)
    if style[0] == 'k':
        color = 'black'
    elif style[0] == 'r':
        color = 'red'
    elif style[0] == 'g':
        color = 'green'
    elif style[0] == 'b':
        color = 'blue'
    elif style[0] == 'y':
        color = 'yellow'
    elif style[0] == 'm':
        color = 'magenta'
    elif style[0] == 'c':
        color = 'cyan'
    elif style[0] == 'w':
        color = 'white'
    else:
        raise ValueError(f"Unknown color code: {style[0]}")

    # Check for the marker style (rest of the string)
    for char in style[1:]:
        if char == 'o':
            marker = 'o'  # Circle
        elif char == 's':
            marker = 's'  # Square
        elif char == 'd':
            marker = 'd'  # Diamond
        elif char == '^':
            marker = '^'  # Triangle up
        elif char == 'v':
            marker = 'v'  # Triangle down
        elif char == '>':
            marker = '>'  # Triangle right
        elif char == '<':
            marker = '<'  # Triangle left
        elif char == '*':
            marker = '*'  # Star
        elif char == '+':
            marker = '+'  # Plus
        elif char == 'x':
            marker = 'x'  # X
        # Ignore unknown characters (to handle linestyle separately)

    # Check for the line style
    if '--' in style:
        linestyle = '--'  # Dashed
    elif '-.' in style:
        linestyle = '-.'  # Dash-dot
    elif '-' in style:
        linestyle = '-'   # Solid
    elif ':' in style:
        linestyle = ':'   # Dotted
    else:
        linestyle = ''    # No line style

    return color, marker, linestyle


def get_version_string(filename):
    file1 = open(filename,"r")
    Lines = file1.readlines()
    version_str = Lines[0].strip()
    file1.close()
    return version_str


def add_version_string(ax, version_str, plot_type='linear', scale_x=1.00, scale_y=1.02,
                       font_name='Times', font_size=10):
    """
    Adds a version string to a matplotlib plot.

    Parameters:
    ax (matplotlib.axes.Axes): The axes to add the version string to.
    filename (str): Path to the version string file.
    plot_type (str): Type of plot ('loglog', 'semilogx', 'semilogy', or 'linear').
    scale_x (float): Scaling factor for X-position.
    scale_y (float): Scaling factor for Y-position.
    font_name (str): Font name for the text.
    font_interpreter (str): Interpreter type (not applicable in matplotlib, kept for compatibility).
    font_size (int): Font size for the text.
    """
    import logging
    # Suppress just the 'findfont' warnings from matplotlib's font manager
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    if (version_str):

        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()

        eps = 1.e-10
        if plot_type == 'loglog':
            x_lim = np.maximum(eps,x_lim)
            y_lim = np.maximum(eps,y_lim)
            x_pos = 10**(np.log10(x_lim[0]) + scale_x * (np.log10(x_lim[1]) - np.log10(x_lim[0])))
            y_pos = 10**(np.log10(y_lim[0]) + scale_y * (np.log10(y_lim[1]) - np.log10(y_lim[0])))
        elif plot_type == 'semilogx':
            x_lim = np.maximum(eps,x_lim)
            x_pos = 10**(np.log10(x_lim[0]) + scale_x * (np.log10(x_lim[1]) - np.log10(x_lim[0])))
            y_pos = y_lim[0] + scale_y * (y_lim[1] - y_lim[0])
        elif plot_type == 'semilogy':
            x_pos = x_lim[0] + scale_x * (x_lim[1] - x_lim[0])
            y_lim = np.maximum(eps,y_lim)
            y_pos = 10**(np.log10(y_lim[0]) + scale_y * (np.log10(y_lim[1]) - np.log10(y_lim[0])))
        else:
            x_pos = x_lim[0] + scale_x * (x_lim[1] - x_lim[0])
            y_pos = y_lim[0] + scale_y * (y_lim[1] - y_lim[0])

        ax.text(x_pos, y_pos, version_str, fontsize=font_size, fontname=font_name, verticalalignment='bottom', horizontalalignment='right')



def get_plot_style(style="fds"):
    """
    Returns a dictionary of plot style parameters based on the specified style.

    Parameters:
    - style (str): The style to use ('fds', 'paper', etc.). Default is 'fds'.

    Returns:
    - dict: A dictionary containing plot style parameters.
    """
    import logging
    # Suppress just the 'findfont' warnings from matplotlib's font manager
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

    if style == "fds":
        return {
            # Font properties
            "Font_Name": "Times",
            "Font_Interpreter": "TeX",
            "Key_Font_Size": 12,
            "Title_Font_Size": 16,
            "Label_Font_Size": 16,
            "Scat_Title_Font_Size": 14,
            "Scat_Label_Font_Size": 14,
            "Marker_Size": 4,
            "D1_Marker_Size": 4,
            "D2_Marker_Size": 4,

            # Line properties
            "Line_Width": 1.0,
            "D1_Line_Width": 1.0,
            "D2_Line_Width": 1.0,

            # Plot properties
            "Plot_Units": "inches",
            "Plot_Width": 5.0,
            "Plot_Height": 3.4,
            "Plot_X": 1.15,
            "Plot_Y": 0.75,
            "Scat_Plot_Width": 4.75,
            "Scat_Plot_Height": 4.75,
            "Scat_Plot_X": 1.00,
            "Scat_Plot_Y": 0.75,
            "Subtitle_Text_Offset": 0.05,
            "VerStr_Scale_X": 0.60,
            "VerStr_Scale_Y": 1.05,

            # Paper properties
            "Paper_Units": "inches",
            "Paper_Width": 6.5,
            "Paper_Height": 4.5,
            "Scat_Paper_Height": 6.0,
            "Scat_Paper_Width": 6.0,

            # Print properties
            "Figure_Visibility": "off",
            "Image_File_Type": "-dpdf",
        }

    elif style == "paper":
        return {
            # Font properties
            "Font_Name": "Helvetica",
            "Font_Interpreter": "TeX",
            "Key_Font_Size": 16,
            "Title_Font_Size": 20,
            "Label_Font_Size": 20,
            "Scat_Title_Font_Size": 14,
            "Scat_Label_Font_Size": 14,
            "Marker_Size": 10,
            "D1_Marker_Size": 10,
            "D2_Marker_Size": 10,

            # Line properties
            "Line_Width": 1.0,
            "D1_Line_Width": 1.0,
            "D2_Line_Width": 1.0,

            # Plot properties
            "Plot_Units": "normalized",
            "Plot_X": 0.1500,
            "Plot_Y": 0.1500,
            "Plot_Width": 0.7750,
            "Plot_Height": 0.8150 * 0.95,  # Adjusted for exponential y-axis tick labels
            "Scat_Plot_Width": 4.75,
            "Scat_Plot_Height": 4.75,
            "Scat_Plot_X": 0.75,
            "Scat_Plot_Y": 0.75,
            "Subtitle_Text_Offset": 0.05,
            "VerStr_Scale_X": 0.60,
            "VerStr_Scale_Y": 1.05,

            # Paper properties
            "Paper_Units": "inches",
            "Paper_Width": 8.0,
            "Paper_Height": 6.0,
            "Scat_Paper_Height": 6.0,
            "Scat_Paper_Width": 6.0,

            # Print properties
            "Figure_Visibility": "on",
            "Image_File_Type": "-dpdf",
        }

    else:
        raise ValueError(f"Unknown style '{style}'. Please choose 'fds' or 'paper'.")


def matlab_legend_to_matplotlib(position):
    """
    Convert a MATLAB legend position string to the corresponding matplotlib location string.

    Parameters:
        position (str): MATLAB-style legend position (e.g., 'northeast', 'southwestoutside')

    Returns:
        str: Matplotlib-compatible legend location (e.g., 'upper right')
    """
    mapping = {
        'north': 'upper center',
        'south': 'lower center',
        'east': 'center right',
        'west': 'center left',
        'northeast': 'upper right',
        'southeast': 'lower right',
        'southwest': 'lower left',
        'northwest': 'upper left',
        'northeastoutside': 'center left',    # rough equivalent
        'northwestoutside': 'center right',
        'southeastoutside': 'center left',
        'southwestoutside': 'center right',
        'best': 'best'
    }

    if not isinstance(position, str):
        return 'best'

    return mapping.get(position.strip().lower(), 'best')


def define_plot_parameters(D, irow, lightweight=False):
    import numpy as np

    class plot_parameters:
        def __init__(self):
            pass

        def __repr__(self):
            return str(self.__dict__)

    # --- FAST PATH ----------------------------------------------------------
    if lightweight:
        col_idx = {col: i for i, col in enumerate(D.columns)}
        row = D.iloc[irow].values

        def get(col, default=None):
            idx = col_idx.get(col)
            return row[idx] if idx is not None else default

        d = plot_parameters()

        # Core identifiers
        d.switch_id       = get('switch_id')
        d.Dataname        = get('Dataname')
        d.VerStr_Filename = get('VerStr_Filename')
        d.Plot_Filename   = get('Plot_Filename')
        d.Plot_Title      = get('Plot_Title')
        d.Quantity        = get('Quantity')
        d.Metric          = get('Metric')
        d.Error_Tolerance = get('Error_Tolerance')

        # File and column info
        d.d1_Filename       = get('d1_Filename')
        d.d1_Col_Name_Row   = get('d1_Col_Name_Row', 1)
        d.d1_Data_Row       = get('d1_Data_Row', 2)
        d.d1_Ind_Col_Name   = get('d1_Ind_Col_Name')
        d.d1_Dep_Col_Name   = get('d1_Dep_Col_Name')
        d.d1_Key            = get('d1_Key', '')
        d.d1_Style          = get('d1_Style', '')
        d.d1_Comp_Start     = get('d1_Comp_Start', np.nan)
        d.d1_Comp_End       = get('d1_Comp_End', np.nan)
        d.d1_Dep_Comp_Start = get('d1_Dep_Comp_Start', np.nan)
        d.d1_Dep_Comp_End   = get('d1_Dep_Comp_End', np.nan)
        d.d1_Initial_Value  = get('d1_Initial_Value', 0.0)

        d.d2_Filename       = get('d2_Filename')
        d.d2_Col_Name_Row   = get('d2_Col_Name_Row', 1)
        d.d2_Data_Row       = get('d2_Data_Row', 2)
        d.d2_Ind_Col_Name   = get('d2_Ind_Col_Name')
        d.d2_Dep_Col_Name   = get('d2_Dep_Col_Name')
        d.d2_Key            = get('d2_Key', '')
        d.d2_Style          = get('d2_Style', '')
        d.d2_Comp_Start     = get('d2_Comp_Start', np.nan)
        d.d2_Comp_End       = get('d2_Comp_End', np.nan)
        d.d2_Dep_Comp_Start = get('d2_Dep_Comp_Start', np.nan)
        d.d2_Dep_Comp_End   = get('d2_Dep_Comp_End', np.nan)
        d.d2_Initial_Value  = get('d2_Initial_Value', 0.0)

        # Plot formatting
        d.Ind_Title        = get('Ind_Title', '')
        d.Dep_Title        = get('Dep_Title', '')
        d.Min_Ind          = get('Min_Ind')
        d.Max_Ind          = get('Max_Ind')
        d.Min_Dep          = get('Min_Dep')
        d.Max_Dep          = get('Max_Dep')
        d.Scale_Ind        = get('Scale_Ind', 1.0)
        d.Scale_Dep        = get('Scale_Dep', 1.0)
        d.Flip_Axis        = get('Flip_Axis', '')
        d.Plot_Type        = get('Plot_Type', 'linear')
        d.Key_Position     = get('Key_Position', 'best')
        d.Title_Position   = get('Title_Position', '')
        d.Legend_XYWidthHeight = get('Legend_XYWidthHeight', '')
        d.Paper_Width_Factor   = get('Paper_Width_Factor', 1.0)

        # Grouping / style info
        d.Group_Key_Label  = get('Group_Key_Label')
        d.Group_Style      = get('Group_Style')
        d.Fill_Color       = get('Fill_Color')
        d.Font_Interpreter = get('Font_Interpreter')

        # --- sanitization for human-facing strings ---
        d.Plot_Title      = sanitize(safe_strip(d.Plot_Title))
        d.Ind_Title       = sanitize(safe_strip(d.Ind_Title))
        d.Dep_Title       = sanitize(safe_strip(d.Dep_Title))
        d.Quantity        = sanitize(safe_strip(d.Quantity))
        d.Metric          = sanitize(safe_strip(d.Metric))
        d.Group_Key_Label = sanitize(safe_strip(d.Group_Key_Label))
        d.d1_Key          = sanitize(safe_strip(d.d1_Key))
        d.d2_Key          = sanitize(safe_strip(d.d2_Key))

        return d

    # --- FULL PATH ----------------------------------------------------------
    class plot_parameters_full(plot_parameters):
        def __init__(self):
            self.switch_id            = D.values[irow,D.columns.get_loc('switch_id')]
            self.Dataname             = D.values[irow,D.columns.get_loc('Dataname')]
            self.VerStr_Filename      = D.values[irow,D.columns.get_loc('VerStr_Filename')]
            self.d1_Filename          = D.values[irow,D.columns.get_loc('d1_Filename')]
            self.d1_Col_Name_Row      = D.values[irow,D.columns.get_loc('d1_Col_Name_Row')]
            self.d1_Data_Row          = D.values[irow,D.columns.get_loc('d1_Data_Row')]
            self.d1_Ind_Col_Name      = D.values[irow,D.columns.get_loc('d1_Ind_Col_Name')]
            self.d1_Dep_Col_Name      = D.values[irow,D.columns.get_loc('d1_Dep_Col_Name')]
            self.d1_Key               = D.values[irow,D.columns.get_loc('d1_Key')]
            self.d1_Style             = D.values[irow,D.columns.get_loc('d1_Style')]
            self.d1_Start             = D.values[irow,D.columns.get_loc('d1_Start')]
            self.d1_End               = D.values[irow,D.columns.get_loc('d1_End')]
            self.d1_Tick              = D.values[irow,D.columns.get_loc('d1_Tick')]
            self.d1_Comp_Start        = D.values[irow,D.columns.get_loc('d1_Comp_Start')]
            self.d1_Comp_End          = D.values[irow,D.columns.get_loc('d1_Comp_End')]
            self.d1_Dep_Comp_Start    = D.values[irow,D.columns.get_loc('d1_Dep_Comp_Start')]
            self.d1_Dep_Comp_End      = D.values[irow,D.columns.get_loc('d1_Dep_Comp_End')]
            self.d1_Initial_Value     = D.values[irow,D.columns.get_loc('d1_Initial_Value')]
            self.d2_Filename          = D.values[irow,D.columns.get_loc('d2_Filename')]
            self.d2_Col_Name_Row      = D.values[irow,D.columns.get_loc('d2_Col_Name_Row')]
            self.d2_Data_Row          = D.values[irow,D.columns.get_loc('d2_Data_Row')]
            self.d2_Ind_Col_Name      = D.values[irow,D.columns.get_loc('d2_Ind_Col_Name')]
            self.d2_Dep_Col_Name      = D.values[irow,D.columns.get_loc('d2_Dep_Col_Name')]
            self.d2_Key               = D.values[irow,D.columns.get_loc('d2_Key')]
            self.d2_Style             = D.values[irow,D.columns.get_loc('d2_Style')]
            self.d2_Start             = D.values[irow,D.columns.get_loc('d2_Start')]
            self.d2_End               = D.values[irow,D.columns.get_loc('d2_End')]
            self.d2_Tick              = D.values[irow,D.columns.get_loc('d2_Tick')]
            self.d2_Comp_Start        = D.values[irow,D.columns.get_loc('d2_Comp_Start')]
            self.d2_Comp_End          = D.values[irow,D.columns.get_loc('d2_Comp_End')]
            self.d2_Dep_Comp_Start    = D.values[irow,D.columns.get_loc('d2_Dep_Comp_Start')]
            self.d2_Dep_Comp_End      = D.values[irow,D.columns.get_loc('d2_Dep_Comp_End')]
            self.d2_Initial_Value     = D.values[irow,D.columns.get_loc('d2_Initial_Value')]
            self.Plot_Title           = D.values[irow,D.columns.get_loc('Plot_Title')]
            self.Ind_Title            = D.values[irow,D.columns.get_loc('Ind_Title')]
            self.Dep_Title            = D.values[irow,D.columns.get_loc('Dep_Title')]
            self.Min_Ind              = D.values[irow,D.columns.get_loc('Min_Ind')]
            self.Max_Ind              = D.values[irow,D.columns.get_loc('Max_Ind')]
            self.Scale_Ind            = D.values[irow,D.columns.get_loc('Scale_Ind')]
            self.Min_Dep              = D.values[irow,D.columns.get_loc('Min_Dep')]
            self.Max_Dep              = D.values[irow,D.columns.get_loc('Max_Dep')]
            self.Scale_Dep            = D.values[irow,D.columns.get_loc('Scale_Dep')]
            self.Flip_Axis            = D.values[irow,D.columns.get_loc('Flip_Axis')]
            self.Title_Position       = D.values[irow,D.columns.get_loc('Title_Position')]
            self.Key_Position         = D.values[irow,D.columns.get_loc('Key_Position')]
            self.Legend_XYWidthHeight = D.values[irow,D.columns.get_loc('Legend_XYWidthHeight')]
            self.Paper_Width_Factor   = D.values[irow,D.columns.get_loc('Paper_Width_Factor')]
            self.Plot_Type            = D.values[irow,D.columns.get_loc('Plot_Type')]
            self.Plot_Filename        = D.values[irow,D.columns.get_loc('Plot_Filename')]
            self.Quantity             = D.values[irow,D.columns.get_loc('Quantity')]
            self.Metric               = D.values[irow,D.columns.get_loc('Metric')]
            self.Error_Tolerance      = D.values[irow,D.columns.get_loc('Error_Tolerance')]
            self.Group_Key_Label      = D.values[irow,D.columns.get_loc('Group_Key_Label')]
            self.Group_Style          = D.values[irow,D.columns.get_loc('Group_Style')]
            self.Fill_Color           = D.values[irow,D.columns.get_loc('Fill_Color')]
            self.Font_Interpreter     = D.values[irow,D.columns.get_loc('Font_Interpreter')]

    d = plot_parameters_full()

    # --- sanitization block (unchanged) ---
    d.Plot_Title      = sanitize(safe_strip(d.Plot_Title))
    d.Ind_Title       = sanitize(safe_strip(d.Ind_Title))
    d.Dep_Title       = sanitize(safe_strip(d.Dep_Title))
    d.Quantity        = sanitize(safe_strip(d.Quantity))
    d.Metric          = sanitize(safe_strip(d.Metric))
    d.Group_Key_Label = sanitize(safe_strip(d.Group_Key_Label))
    d.d1_Key          = sanitize(safe_strip(d.d1_Key))
    d.d2_Key          = sanitize(safe_strip(d.d2_Key))

    return d


def define_qrow_variables(Q, j):
    """
    Define scatterplot parameters from the Scatterplot_Inputs.csv row j.

    Mirrors the MATLAB 'define_qrow_variables.m' behavior and returns
    a simple Python object with all scatterplot configuration fields.

    Parameters
    ----------
    Q : pandas.DataFrame
        The scatterplot input file loaded by pandas.read_csv().
    j : int
        Row index (0-based).

    Returns
    -------
    q : object
        Object with attributes corresponding to scatterplot input fields.
    """

    class qrow:
        def __init__(self):
            self.Scatter_Plot_Title  = Q.loc[j, "Scatter_Plot_Title"]
            self.Ind_Title           = Q.loc[j, "Ind_Title"]
            self.Dep_Title           = Q.loc[j, "Dep_Title"]
            self.Plot_Min            = Q.loc[j, "Plot_Min"]
            self.Plot_Max            = Q.loc[j, "Plot_Max"]
            self.Title_Position      = Q.loc[j, "Title_Position"]
            self.Key_Position        = Q.loc[j, "Key_Position"]
            self.Paper_Width_Factor  = Q.loc[j, "Paper_Width_Factor"]
            self.Sigma_E             = Q.loc[j, "Sigma_E"]
            self.Weight_Data         = Q.loc[j, "Weight_Data"]
            self.Plot_Type           = Q.loc[j, "Plot_Type"]
            self.Plot_Filename       = Q.loc[j, "Plot_Filename"]

        def __repr__(self):
            return str(self.__dict__)

    q = qrow()

    # Sanitize only human-readable fields
    q.Scatter_Plot_Title = sanitize(safe_strip(q.Scatter_Plot_Title))
    q.Ind_Title = sanitize(safe_strip(q.Ind_Title))
    q.Dep_Title = sanitize(safe_strip(q.Dep_Title))
    q.Plot_Filename = safe_strip(q.Plot_Filename)
    q.Key_Position = safe_strip(q.Key_Position)
    q.Plot_Type = safe_strip(q.Plot_Type)

    # Parse numeric fields
    def to_float(val):
        try:
            return float(val)
        except Exception:
            return np.nan

    q.Plot_Min = to_float(q.Plot_Min)
    q.Plot_Max = to_float(q.Plot_Max)
    q.Paper_Width_Factor = to_float(q.Paper_Width_Factor)
    q.Sigma_E = to_float(q.Sigma_E)

    # Parse Title_Position as [x, y] floats
    if isinstance(q.Title_Position, str):
        try:
            vals = [float(x) for x in q.Title_Position.split()]
            if len(vals) == 2:
                q.Title_Position = vals
            else:
                q.Title_Position = [0.03, 0.95]
        except Exception:
            q.Title_Position = [0.03, 0.95]
    else:
        q.Title_Position = [0.03, 0.95]

    # Normalize Weight_Data (yes/no → bool)
    if isinstance(q.Weight_Data, str):
        q.Weight_Data = q.Weight_Data.strip().lower() == "yes"
    else:
        q.Weight_Data = bool(q.Weight_Data)

    return q


def sanitize(text: str) -> str:
    """Escape LaTeX specials outside math mode ($...$)."""
    if not isinstance(text, str):
        return text

    specials = {
        "&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#",
        "$": r"\$", "{": r"\{", "}": r"\}", "^": r"\^{}", "~": r"\~{}",
    }

    # Split into math and text segments
    import re
    parts = re.split(r"(\$.*?\$)", text)
    sanitized = []
    for part in parts:
        if part.startswith("$") and part.endswith("$"):
            sanitized.append(part)  # math untouched
        else:
            s = part
            for k, v in specials.items():
                s = s.replace(k, v)
            sanitized.append(s)
    return "".join(sanitized)


def safe_strip(val):
    """Strip whitespace safely from strings; return empty string otherwise."""
    return val.strip() if isinstance(val, str) else ""


def scatplot(saved_data, drange, **kwargs):
    """
    Generate scatter plots and compute validation/verification statistics.
    Faithful translation of MATLAB scatplot.m behavior:
      - validation_statistics.tex
      - validation_histograms.tex
      - validation_scatterplot_output.csv
    """
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import fdsplotlib

    Manuals_Dir = kwargs.get("Manuals_Dir", "")
    Scatterplot_Inputs_File = kwargs.get("Scatterplot_Inputs_File", "")
    Stats_Output = kwargs.get("Stats_Output", "Validation")
    Scatterplot_Dir = kwargs.get("Scatterplot_Dir", "")
    verbose = kwargs.get("verbose", True)

    plot_style = get_plot_style("fds")
    scat_figure_size = (plot_style["Scat_Paper_Width"],plot_style["Scat_Paper_Height"])
    scat_plot_size = (plot_style["Scat_Plot_Width"],plot_style["Scat_Plot_Height"])
    scat_plot_origin = (plot_style["Scat_Plot_X"],plot_style["Scat_Plot_Y"])

    if not os.path.exists(Scatterplot_Inputs_File):
        raise FileNotFoundError(f"[scatplot] Missing input file: {Scatterplot_Inputs_File}")
    os.makedirs(Scatterplot_Dir, exist_ok=True)

    # Output filenames (same logic as MATLAB)
    if Stats_Output.lower() == "validation":
        Statistics_Tex_Output = os.path.join(Scatterplot_Dir, "validation_statistics.tex")
        Output_File = os.path.join(Scatterplot_Dir, "validation_scatterplot_output.csv")
        Histogram_Tex_Output = os.path.join(Scatterplot_Dir, "validation_histograms.tex")
    elif Stats_Output.lower() == "verification":
        Statistics_Tex_Output = os.path.join(Scatterplot_Dir, "verification_statistics.tex")
        Output_File = os.path.join(Scatterplot_Dir, "verification_scatterplot_output.csv")
        Histogram_Tex_Output = os.path.join(Scatterplot_Dir, "verification_histograms.tex")
    else:
        Statistics_Tex_Output = os.path.join(Scatterplot_Dir, f"Scatterplot_Tables_{Stats_Output}.tex")
        Output_File = os.path.join(Scatterplot_Dir, f"Scatterplot_Stats_{Stats_Output}.csv")
        Histogram_Tex_Output = os.path.join(Scatterplot_Dir, f"Scatterplot_Histograms_{Stats_Output}.tex")

    # --- Unpack saved_data (dataplot output) ---
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
    ) = saved_data

    Q = pd.read_csv(Scatterplot_Inputs_File)
    if verbose:
        print(f"[scatplot] Loaded {len(Q)} scatterplot definitions")

    Output_Histograms = []

    if Stats_Output.lower() == "verification":
        output_stats = [[
            "Dataplot Line Number",
            "Verification Group",
            "Case Name",
            "Type of Metric",
            "Expected Quantity",
            "Expected Value",
            "Predicted Quantity",
            "Predicted Value",
            "Dependent Variable",
            "Type of Error",
            "Error",
            "Error Tolerance",
            "Within Specified Error Tolerance",
            "Plot Filename",
        ]]
    else:  # validation
        output_stats = [[
            "Quantity",
            "Number of Datasets",
            "Number of Points",
            "Sigma_Experiment",
            "Sigma_Model",
            "Bias",
        ]]

    for _, row in Q.iterrows():
        plt.close('all')
        plt.figure().clear()

        Scatter_Plot_Title = row["Scatter_Plot_Title"]
        Plot_Filename = row["Plot_Filename"]
        Plot_Min = float(row["Plot_Min"])
        Plot_Max = float(row["Plot_Max"])
        Plot_Type = str(row["Plot_Type"]).strip().lower()

        # --- MATLAB parity: Sigma_E only required for Validation ---
        if Stats_Output.lower() == "validation":
            Sigma_E_input = float(row["Sigma_E"]) if "Sigma_E" in row and not pd.isna(row["Sigma_E"]) else 0.0
        else:
            Sigma_E_input = 0.0

        #if verbose:
        #    print(f"[scatplot] Processing {Scatter_Plot_Title}")

        # Match dataplot entries
        match_idx = [i for i, q in enumerate(Save_Quantity)
                     if Scatter_Plot_Title.strip().lower() in str(q).lower()]
        if not match_idx:
            print(f"[scatplot] No dataplot entries for {Scatter_Plot_Title}")
            continue

        # --- Split logic: Verification vs Validation ---
        # --- VERIFICATION branch ---
        if Stats_Output.lower() == "verification":

            # Loop through each dataplot entry that matches this scatterplot
            for idx in match_idx:

                # --- Extract measured/predicted numeric values ---
                mvals = np.array(Save_Measured_Metric[idx], dtype=float).flatten()
                pvals = np.array(Save_Predicted_Metric[idx], dtype=float).flatten()

                # Keep only finite pairs
                mask = np.isfinite(mvals) & np.isfinite(pvals)
                mvals = mvals[mask]
                pvals = pvals[mask]

                # --- Extract quantity labels (MATLAB behavior: pipe = multiple quantities) ---
                def _split_pipe_list(x):
                    """
                    MATLAB dataplot uses '|' to separate multiple quantity labels.
                    Convert "a|b|c" → ["a", "b", "c"] exactly.
                    """
                    if x is None:
                        return [""]

                    # Try to extract a single string out of an object/array
                    try:
                        s = str(np.ravel(np.array(x, dtype=object))[0])
                    except Exception:
                        s = str(x)

                    # Properly split on pipes and strip whitespace
                    return [item.strip() for item in s.split("|")]

                meas_labels = _split_pipe_list(Save_Measured_Quantity[idx])
                pred_labels = _split_pipe_list(Save_Predicted_Quantity[idx])

                # Assign label by index (MATLAB style)
                def _label_at(k, labels):
                    if len(labels) == 0:
                        return ""
                    if len(labels) == 1:
                        return labels[0]
                    if k < len(labels):
                        return labels[k]
                    return labels[-1]  # MATLAB tolerance fallback

                # Case metadata
                case_name = Save_Dataname[idx]
                group     = Save_Group_Key_Label[idx]
                metric    = Save_Metric_Type[idx]
                depvar    = Save_Dep_Title[idx]
                err_tol   = float(Save_Error_Tolerance[idx] or 0.0)
                err_type  = str(Save_Quantity[idx] or "")
                plot_file = Save_Plot_Filename[idx]

                # --- Generate one row per point ---
                for k in range(min(len(mvals), len(pvals))):

                    m = float(mvals[k])
                    p = float(pvals[k])
                    mq = _label_at(k, meas_labels)
                    pq = _label_at(k, pred_labels)

                    # Compute error (Relative or Absolute)
                    if err_type.lower().startswith("relative"):
                        err = abs((p - m) / m) if m != 0 else np.nan
                    else:
                        err = abs(p - m)

                    within = "Yes" if (np.isfinite(err) and err <= err_tol) else "Out of Tolerance"

                    # Order MUST match MATLAB exactly
                    output_stats.append([
                        idx + 1,           # Dataplot line number
                        group,             # Verification group
                        case_name,         # Case name (string)
                        metric,            # Metric type (e.g. "end")
                        mq,                # Expected quantity label (MATLAB: one per row)
                        m,                 # Expected value
                        pq,                # Predicted quantity label
                        p,                 # Predicted value
                        depvar,            # Dependent variable title
                        err_type,          # Type of Error (Relative or Absolute)
                        f"{err:1.2e}",     # Error
                        f"{err_tol:1.2e}", # Error Tolerance
                        within,            # Within Spec?
                        plot_file,         # Plot filename
                    ])

            continue  # Move to next scatterplot definition

        Measured_Values = np.concatenate(  [np.ravel(np.array(Save_Measured_Metric[i], dtype=float)) for i in match_idx]  )
        Predicted_Values = np.concatenate(  [np.ravel(np.array(Save_Predicted_Metric[i], dtype=float)) for i in match_idx]  )

        # --- Ensure equal-length measured and predicted arrays before masking ---
        m_len = len(Measured_Values)
        p_len = len(Predicted_Values)
        if m_len != p_len:
            print(f"[scatplot] Skipping '{Scatter_Plot_Title}' "
                  f"due to length mismatch: measured={m_len}, predicted={p_len}")
            print(f"   Example Measured sample: {Measured_Values[:5]}")
            print(f"   Example Predicted sample: {Predicted_Values[:5]}")
            continue  # skip to next scatterplot, as MATLAB does

        mask = np.isfinite(Measured_Values) & np.isfinite(Predicted_Values)
        Measured_Values = Measured_Values[mask]
        Predicted_Values = Predicted_Values[mask]


        # --- Call histogram BEFORE filtering to match MATLAB (includes zeros/Infs) ---
        try:
            hist_file, pval = statistics_histogram(
                Measured_Values, Predicted_Values,
                Plot_Filename, Manuals_Dir, Scatter_Plot_Title
            )
            if hist_file:
                Output_Histograms.append(hist_file)
        except Exception as e:
            print(f"[scatplot] Histogram error for {Scatter_Plot_Title}: {e}")

        # --- Now filter for plotting & statistics ---
        in_range = (
            (Measured_Values >= Plot_Min) & (Measured_Values <= Plot_Max) &
            (Predicted_Values >= Plot_Min) & (Predicted_Values <= Plot_Max)
        )
        positive = (Measured_Values > 0) & (Predicted_Values > 0)
        mask = in_range & positive
        Measured_Values  = Measured_Values[mask]
        Predicted_Values = Predicted_Values[mask]

        if len(Measured_Values) == 0:
            print(f"[scatplot] Skipping {Scatter_Plot_Title} (no valid data)")
            continue

        # --- Compute statistics ---
        weight = np.ones_like(Measured_Values)
        log_E_bar = np.sum(np.log(Measured_Values) * weight) / np.sum(weight)
        log_M_bar = np.sum(np.log(Predicted_Values) * weight) / np.sum(weight)
        u2 = np.sum((((np.log(Predicted_Values) - np.log(Measured_Values))
                      - (log_M_bar - log_E_bar)) ** 2) * weight) / (np.sum(weight) - 1)
        u = np.sqrt(u2)
        Sigma_E = min(u / np.sqrt(2), Sigma_E_input / 100.0)
        Sigma_M = np.sqrt(max(0.0, u ** 2 - Sigma_E ** 2))
        delta = np.exp(log_M_bar - log_E_bar + 0.5 * Sigma_M ** 2 - 0.5 * Sigma_E ** 2)

        # --- Scatter Plot ---
        fig = fdsplotlib.plot_to_fig(x_data=[-1], y_data=[-1],
                                figure_size=scat_figure_size,
                                plot_size=scat_plot_size,
                                plot_origin=scat_plot_origin,
                                plot_type=Plot_Type,
                                x_min=Plot_Min, x_max=Plot_Max, y_min=Plot_Min, y_max=Plot_Max,
                                x_label=row["Ind_Title"],
                                y_label=row["Dep_Title"],
                                legend_location='outside',
                                legend_expand=row["Paper_Width_Factor"])

        fdsplotlib.plot_to_fig(x_data=[Plot_Min, Plot_Max], y_data=[Plot_Min, Plot_Max], line_style="k-", figure_handle=fig)
        if Sigma_E > 0:
            fdsplotlib.plot_to_fig(x_data=[Plot_Min, Plot_Max], y_data=np.array([Plot_Min, Plot_Max]) * (1 + 2 * Sigma_E), line_style="k--", figure_handle=fig)
            fdsplotlib.plot_to_fig(x_data=[Plot_Min, Plot_Max], y_data=np.array([Plot_Min, Plot_Max]) / (1 + 2 * Sigma_E), line_style="k--", figure_handle=fig)
            fdsplotlib.plot_to_fig(x_data=[Plot_Min, Plot_Max], y_data=np.array([Plot_Min, Plot_Max]) * delta, line_style="r-", figure_handle=fig)
            fdsplotlib.plot_to_fig(x_data=[Plot_Min, Plot_Max], y_data=np.array([Plot_Min, Plot_Max]) * delta * (1 + 2 * Sigma_M), line_style="r--", figure_handle=fig)
            fdsplotlib.plot_to_fig(x_data=[Plot_Min, Plot_Max], y_data=np.array([Plot_Min, Plot_Max]) * delta / (1 + 2 * Sigma_M), line_style="r--", figure_handle=fig)

        # --- Plot each dataset with its own marker and color ---
        seen_labels = set()

        for idx in match_idx:
            style = str(Save_Group_Style[idx]).strip() if Save_Group_Style[idx] else "ko"
            fill = str(Save_Fill_Color[idx]).strip() if Save_Fill_Color[idx] else "none"
            label = str(Save_Group_Key_Label[idx]).strip() if Save_Group_Key_Label[idx] else ""

            # Flatten valid points for this dataset
            mvals = np.array(Save_Measured_Metric[idx], dtype=float).flatten()
            pvals = np.array(Save_Predicted_Metric[idx], dtype=float).flatten()

            if len(mvals) != len(pvals):
                print(f"[DEBUG] Mismatch for {Scatter_Plot_Title} @ idx={idx}: "
                      f"Measured={len(mvals)}, Predicted={len(pvals)}, "
                      f"Group={Save_Group_Key_Label[idx]}")
                print(f"   Measured metric sample: {mvals[:5]}")
                print(f"   Predicted metric sample: {pvals[:5]}")

            mask = (
                (mvals >= Plot_Min) & (mvals <= Plot_Max) &
                (pvals >= Plot_Min) & (pvals <= Plot_Max) &
                (mvals > 0) & (pvals > 0)
            )
            mvals = mvals[mask]
            pvals = pvals[mask]

            if len(mvals) == 0:
                continue

            # Only assign a legend label once per experiment
            data_label = label if label and label not in seen_labels else None
            if label:
                seen_labels.add(label)

            fdsplotlib.plot_to_fig(
                x_data=mvals,
                y_data=pvals,
                marker_style=style,
                marker_fill_color=fill,
                figure_handle=fig,
                data_label=data_label,
            )

        pdf_path = os.path.join(Manuals_Dir, Plot_Filename + ".pdf")
        os.makedirs(os.path.dirname(pdf_path), exist_ok=True)
        fig.savefig(pdf_path)
        plt.close(fig)
        plt.figure().clear()

        # --- Collect statistics for CSV/TeX ---
        group_labels = []
        for i in match_idx:
            try:
                gl = Save_Group_Key_Label[i]
                if isinstance(gl, (list, tuple)):
                    gl = gl[0] if len(gl) > 0 else ""
                gl = str(gl).strip() if gl is not None else ""
                group_labels.append(gl)
            except Exception:
                pass
        unique_datasets = len(set(g for g in group_labels if g))

        output_stats.append([
            Scatter_Plot_Title,
            unique_datasets,
            len(Measured_Values),
            f"{Sigma_E:.2f}",
            f"{Sigma_M:.2f}",
            f"{delta:.2f}",
        ])

    # --- Write statistics outputs ---
    statistics_output(
        Stats_Output=Stats_Output,
        output_stats=output_stats,
        Output_File=Output_File,
        Statistics_Tex_Output=Statistics_Tex_Output,
        Histogram_Tex_Output=Histogram_Tex_Output,
        Output_Histograms=Output_Histograms,
    )

    print("[scatplot] Completed successfully.")
    return


def statistics_output(
    Stats_Output,
    output_stats,
    Output_File,
    Statistics_Tex_Output=None,
    Histogram_Tex_Output=None,
    Output_Histograms=None,
):
    """
    Python translation of MATLAB statistics_output.m

    - For 'Verification': writes CSV + verification_statistics.tex
    - For 'Validation' : writes CSV + validation_statistics.tex + histograms
    """

    import os
    import pandas as pd
    import numpy as np

    if Stats_Output is None or str(Stats_Output).lower() == "none":
        print("[statistics_output] Skipping (Stats_Output=None)")
        return

    # Ensure output directory exists
    out_dir = os.path.dirname(Output_File)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # VERIFICATION BRANCH  (Stats_Output == 'Verification')
    # ------------------------------------------------------------------
    if str(Stats_Output).lower() == "verification":

        # --- 1) Write CSV exactly from output_stats (header + rows) ---
        header = output_stats[0]
        rows   = output_stats[1:]
        df_csv = pd.DataFrame(rows, columns=header)
        df_csv.to_csv(Output_File, index=False)
        print(f"[statistics_output] Wrote verification CSV: {Output_File}")

        # --- 2) Write LaTeX table (MATLAB-style verification_statistics.tex) ---
        # --- Write LaTeX verification table (MATLAB-faithful) ---
        if Statistics_Tex_Output:

            def _escape(s):
                return (str(s)
                        .replace('\\', '\\textbackslash{}')
                        .replace('_', '\\_')
                        .replace('%', '\\%')
                        .replace('&', '\\&')
                        .replace('#', '\\#')
                        .replace('{', '\\{')
                        .replace('}', '\\}'))

            def _safe_float(x):
                try:
                    return float(x)
                except:
                    return None

            header = output_stats[0]
            rows = output_stats[1:]

            # ! DO NOT SORT — MATLAB preserves dataplot order
            # rows = sorted(rows, key=lambda r: str(r[2]).lower())

            with open(Statistics_Tex_Output, "w") as fid:

                fid.write("\\scriptsize\n")
                fid.write("\\begin{longtable}{|p{2.5in}|l|p{1in}|l|p{1in}|l|l|l|l|l|}\n")
                fid.write("\\hline\n")
                fid.write("Case Name & Section & Expected & Expected & Predicted & Predicted & "
                          "Type of & Error & Error & Within \\\\\n")
                fid.write(" & & Quantity & Value & Quantity & Value & Error &  & Tolerance & Tol. "
                          "\\\\ \\hline \\hline\n")
                fid.write("\\endfirsthead\n\\hline\n")
                fid.write("Case Name & Section & Expected & Expected & Predicted & Predicted & "
                          "Type of & Error & Error & Within \\\\\n")
                fid.write(" & & Quantity & Value & Quantity & Value & Error &  & Tolerance & Tol. "
                          "\\\\ \\hline \\hline\n")
                fid.write("\\endhead\n\\hline\n\\endfoot\n\\hline\n\\endlastfoot\n")

                for r in rows:

                    case = str(r[2])
                    if str(r[13])[:14]=='FDS_User_Guide':
                        section = f"\\ref{{{'UG-'+case}}}"
                    else:
                        section = f"\\ref{{{case}}}"

                    # One row per datapoint; no splitting, no combining
                    exp_q  = _escape(r[4])
                    pred_q = _escape(r[6])

                    exp_val_f = _safe_float(r[5])
                    pred_val_f = _safe_float(r[7])
                    err_val_f = _safe_float(r[10])
                    tol_f = _safe_float(r[11])

                    exp_val  = f"{exp_val_f:1.2e}" if exp_val_f is not None else _escape(r[5])
                    pred_val = f"{pred_val_f:1.2e}" if pred_val_f is not None else _escape(r[7])
                    err_val  = f"{err_val_f:1.2e}" if err_val_f is not None else _escape(r[10])
                    tol_val  = f"{tol_f:1.2e}"     if tol_f is not None else _escape(r[11])

                    err_type = str(r[9]).replace(" Error", "")
                    within   = _escape(r[12])

                    fid.write(
                        f"{_escape(case)} & {section} & "
                        f"{exp_q} & {exp_val} & "
                        f"{pred_q} & {pred_val} & "
                        f"{err_type} & {err_val} & "
                        f"{tol_val} & {within} \\\\\n"
                    )

                fid.write("\\end{longtable}\n\\normalsize\n")

            print(f"[statistics_output] Wrote LaTeX Verification table: {Statistics_Tex_Output}")

        return  # Done with verification branch

    # ------------------------------------------------------------------
    # VALIDATION BRANCH  (unchanged in spirit from your version)
    # ------------------------------------------------------------------
    # Build DataFrame from output_stats
    df = pd.DataFrame(output_stats[1:], columns=output_stats[0])

    # Ensure correct types for the two numeric count columns
    if "Number of Datasets" in df.columns:
        df["Number of Datasets"] = pd.to_numeric(
            df["Number of Datasets"], errors="coerce"
        ).fillna(0).astype(int)

    if "Number of Points" in df.columns:
        df["Number of Points"] = pd.to_numeric(
            df["Number of Points"], errors="coerce"
        ).fillna(0).astype(int)

    # Format last three numeric-looking columns as strings w/ 2 decimals
    for col in ["Sigma_Experiment", "Sigma_Model", "Bias"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").map(
                lambda x: f"{x:0.2f}" if np.isfinite(x) else ""
            )

    # Write CSV
    df.to_csv(Output_File, index=False)
    print(f"[statistics_output] Wrote CSV: {Output_File}")

    # -------- LaTeX Validation Table --------
    if str(Stats_Output).lower() == "validation" and Statistics_Tex_Output:
        with open(Statistics_Tex_Output, "w") as fid:
            fid.write("\\begin{longtable}[c]{|l|c|c|c|c|c|c|}\n")
            fid.write(
                "\\caption[Summary statistics]{Summary statistics for all quantities of interest}\n"
            )
            fid.write("\\label{summary_stats}\n")
            fid.write("\\\\ \\hline\n")
            fid.write(
                "Quantity & Section   & Datasets  & Points    & "
                "$\\widetilde{\\sigma}_{\\rm E}$ & "
                "$\\widetilde{\\sigma}_{\\rm M}$ & Bias "
                "\\\\ \\hline \\hline\n"
            )
            fid.write("\\endfirsthead\n\\hline\n")
            fid.write(
                "Quantity & Section   & Datasets  & Points    & "
                "$\\widetilde{\\sigma}_{\\rm E}$ & "
                "$\\widetilde{\\sigma}_{\\rm M}$ & Bias "
                "\\\\ \\hline \\hline\n"
            )
            fid.write("\\endhead\n")

            for _, r in df.iterrows():
                try:
                    sigma_e = float(r["Sigma_Experiment"])
                    if sigma_e < 0:
                        continue
                    quantity = str(r["Quantity"])
                    section = f"\\ref{{{quantity}}}"
                    fid.write(
                        f"{quantity} & {section} & "
                        f"{int(r['Number of Datasets'])} & "
                        f"{int(r['Number of Points'])} & "
                        f"{float(r['Sigma_Experiment']):0.2f} & "
                        f"{float(r['Sigma_Model']):0.2f} & "
                        f"{float(r['Bias']):0.2f} "
                        "\\\\ \\hline\n"
                    )
                except Exception as e:
                    print(f"[statistics_output] Skipped row due to error: {e}")

            fid.write("\\end{longtable}\n")

        print(f"[statistics_output] Wrote LaTeX Validation table: {Statistics_Tex_Output}")

    # -------- Histogram LaTeX --------
    if str(Stats_Output).lower() == "validation" and Output_Histograms:
        with open(Histogram_Tex_Output, "w") as fid:
            n = len(Output_Histograms)
            pages = int(np.ceil(n / 8.0))
            for i in range(pages):
                fid.write("\\begin{figure}[p]\n")
                fid.write("\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}r}\n")
                for j in range(i * 8, min((i + 1) * 8, n)):
                    end = "&" if (j % 2) == 0 else "\\\\"
                    fid.write(
                        f"\\includegraphics[height=2.2in]"
                        f"{{SCRIPT_FIGURES/Scatterplots/{Output_Histograms[j]}}} {end}\n"
                    )
                fid.write("\\end{tabular*}\n")
                fid.write(f"\\label{{Histogram_{i + 1}}}\n")
                fid.write("\\end{figure}\n\n")
        print(f"[statistics_output] Wrote LaTeX histograms: {Histogram_Tex_Output}")


def histogram_output(Histogram_Tex_Output, Output_Histograms):
    """
    Replicates MATLAB validation_histograms.tex layout exactly.
    Produces 8 histograms per page, 2 columns per row.
    """
    import numpy as np

    if not Output_Histograms:
        print("[histogram_output] No histograms to write.")
        return

    with open(Histogram_Tex_Output, "w") as fid:
        n = len(Output_Histograms)
        pages = int(np.ceil(n / 8))

        for i in range(pages):
            fid.write("\\begin{figure}[p]\n")
            fid.write("\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}r}\n")

            for j in range(i * 8, min((i + 1) * 8, n)):
                end = "&" if j % 2 == 0 else "\\\\"
                fid.write(f"\\includegraphics[height=2.2in]"
                          f"{{SCRIPT_FIGURES/Scatterplots/{Output_Histograms[j]}}} {end}\n")

            fid.write("\\end{tabular*}\n")
            fid.write(f"\\label{{Histogram_{i + 1}}}\n")
            fid.write("\\end{figure}\n\n")

    print(f"[histogram_output] Wrote LaTeX histogram file: {Histogram_Tex_Output}")


def spiegel_test(x):
    """Exact translation of MATLAB spiegel_test.m with Inf/NaN handling identical to MATLAB."""
    import numpy as np
    from math import sqrt, erf

    x = np.asarray(x, dtype=float)

    # MATLAB arithmetic allows Inf and NaN to propagate but doesn't remove them
    xm = np.nanmean(x)
    xs = np.nanstd(x, ddof=0)
    xz = (x - xm) / xs
    xz2 = xz ** 2

    # MATLAB behavior: Inf * 0 -> NaN, which is ignored by sum()
    with np.errstate(invalid='ignore', divide='ignore'):
        term = xz2 * np.log(xz2)

    N = np.nansum(term)
    n = np.sum(np.isfinite(xz2))  # count finite entries only
    ts = (N - 0.73 * n) / (0.8969 * sqrt(n))
    pval = 1 - abs(erf(ts / sqrt(2)))  # 2-sided
    return pval


def statistics_histogram(Measured_Values, Predicted_Values,
                         Plot_Filename, Manuals_Dir, Scatter_Plot_Title,
                         Figure_Visibility='off', Paper_Width=5.0,
                         Paper_Height=3.0, Paper_Units='inches',
                         Image_File_Type='-dpdf'):
    """Faithful translation of statistics_histogram.m (Overholt 2013)."""
    import numpy as np
    import matplotlib.pyplot as plt
    import os

    # --- Compute ln(M/E) exactly as MATLAB
    with np.errstate(divide='ignore', invalid='ignore'):
        ln_M_E = np.log(Predicted_Values) - np.log(Measured_Values)

    # Normality test before filtering — this is the key
    if len(ln_M_E) >= 4:
        pval = spiegel_test(ln_M_E)
    else:
        print(f"[statistics_histogram] Not enough data for {Scatter_Plot_Title}")
        return None, None

    # MATLAB's hist() ignores NaN/Inf implicitly
    valid = np.isfinite(ln_M_E)
    ln_M_E = ln_M_E[valid]

    n, xout = np.histogram(ln_M_E, bins=10)
    xcenters = 0.5 * (xout[:-1] + xout[1:])

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.bar(xcenters, n, width=(xout[1]-xout[0]), linewidth=1,
           edgecolor='k', color=(0.7, 0.7, 0.7))

    dx = xout[1] - xout[0]
    x_lim = [xout[0] - dx, xout[-1] + dx]
    ix = np.arange(x_lim[0], x_lim[1], 1e-3)
    mu = np.mean(ln_M_E)
    sd = np.std(ln_M_E, ddof=0)
    iy = (1 / (sd * np.sqrt(2 * np.pi))) * np.exp(-(ix - mu) ** 2 / (2 * sd ** 2))
    ax.plot(ix, iy * np.trapz(n, xcenters), 'k', linewidth=2)

    ax.set_xlim(x_lim)
    y0, y1 = ax.get_ylim()
    ax.set_ylim([y0, y1 * 1.25])
    ax.set_xlabel("Interval Number")
    ax.set_ylabel("Number of Data Points")
    ax.set_xticks(xcenters)
    ax.set_xticklabels([str(i) for i in range(1, len(xcenters) + 1)])
    ax.text(0.03, 0.90, Scatter_Plot_Title, transform=ax.transAxes)
    ax.text(0.03, 0.82, "Normality Test", transform=ax.transAxes)
    ax.text(0.03, 0.74, f"p-value = {pval:.2f}", transform=ax.transAxes)

    outpath = os.path.join(Manuals_Dir, f"{Plot_Filename}_Histogram.pdf")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    plt.figure().clear()

    print(f"[statistics_histogram] {Scatter_Plot_Title}: p = {pval:.2f}")
    return f"{os.path.basename(Plot_Filename)}_Histogram", pval


def set_ticks_like_matlab(fig):
    ax = fig.axes[0]
    ax.tick_params(axis="both", direction="in", top=True, right=True, width=0.5)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.5)


