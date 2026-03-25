#!/usr/bin/env python3
"""
One-dimensional mesh stretching utility.

The script computes a geometric cell-width sequence after the first m cells,
plots the resulting boundary locations, and writes FDS `&TRNX`, `&TRNY`, or
`&TRNZ` input lines.

Usage notes
-----------
Run with defaults:
    python3 mesh_quad_stretch.py

Write output to a specific directory and print the generated ramp lines:
    python3 mesh_quad_stretch.py --basedir ./out --print-trn

Select a different transform direction and mesh definition:
    python3 mesh_quad_stretch.py --trn TRNX --xs 0.0 --lx 2.0 --m 6 --dxs 0.02 --nx 40

Key outputs
-----------
- `TRN.dat` by default, containing the interior `&TRN*` lines
- `mesh_quad_stretch.png` by default, showing the stretched grid and mapping

Use `--show` to display the figure interactively instead of only saving it.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import tempfile
import time
from pathlib import Path

cache_root = Path(os.environ.get("TMPDIR", tempfile.gettempdir())) / "fds_python_cache"
cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))

import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate 1-D stretched mesh TRN lines and plots."
    )
    parser.add_argument(
        "--basedir",
        default=".",
        help="Directory where the TRN file and figure are written.",
    )
    parser.add_argument(
        "--trn",
        default="TRNZ",
        choices=("TRNX", "TRNY", "TRNZ"),
        help="FDS transform record to emit.",
    )
    parser.add_argument(
        "--trn-id",
        default="MY TRANSFORM",
        help="Transform ID written on each TRN line.",
    )
    parser.add_argument("--xs", type=float, default=1.0, help="Mesh start coordinate.")
    parser.add_argument("--lx", type=float, default=1.0, help="Mesh length.")
    parser.add_argument(
        "--m",
        type=int,
        default=4,
        help="Number of low-side cells with uniform width DXS.",
    )
    parser.add_argument(
        "--dxs",
        type=float,
        default=0.025,
        help="Width of the first m cells.",
    )
    parser.add_argument("--nx", type=int, default=20, help="Total number of cells.")
    parser.add_argument(
        "--etol-rel",
        type=float,
        default=1.0e-10,
        help="Relative convergence tolerance for Newton iteration.",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=1000,
        help="Maximum number of Newton iterations.",
    )
    parser.add_argument(
        "--output",
        default="TRN.dat",
        help="Output file name for TRN lines.",
    )
    parser.add_argument(
        "--figure",
        default="mesh_quad_stretch.png",
        help="Output figure file name.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure after saving it.",
    )
    parser.add_argument(
        "--print-trn",
        action="store_true",
        help="Also print the generated TRN lines to stdout.",
    )
    return parser.parse_args()


def normalize_trn_id(trn_id: str) -> str:
    stripped = trn_id.strip()
    if stripped.startswith("'") and stripped.endswith("'"):
        return stripped
    return f"'{stripped}'"


def validate_inputs(xs: float, lx: float, m: int, dxs: float, nx: int) -> None:
    if lx <= 0.0:
        raise ValueError("Lx must be positive.")
    if dxs <= 0.0:
        raise ValueError("DXS must be positive.")
    if nx < 2:
        raise ValueError("Nx must be at least 2.")
    if m < 1:
        raise ValueError("m must be at least 1.")
    if m >= nx:
        raise ValueError("m must be smaller than Nx so a stretched region exists.")
    if not math.isfinite(xs):
        raise ValueError("xs must be finite.")


def solve_stretch_factor(
    lx: float, m: int, dxs: float, nx: int, etol_rel: float, max_iter: int
) -> tuple[float, int, float, float]:
    print("Solving for stretching parameter..")

    cpu_start = time.process_time()
    an = 1.0
    progress_stride = max(1, math.ceil(max_iter / 1000))

    for n in range(1, max_iter + 1):
        suma = 0.0
        sumap = 0.0
        for i in range(m + 1, nx + 1):
            exponent = i - m
            suma += an**exponent
            sumap += exponent * an ** (exponent - 1)

        f_an = (m - lx / dxs) + suma
        fp_an = sumap
        if fp_an == 0.0:
            raise RuntimeError("Newton iteration hit a zero derivative.")

        an1 = an - f_an / fp_an
        if not math.isfinite(an1):
            raise RuntimeError("Newton iteration produced a non-finite stretch factor.")

        rel_err = abs(an1 - an) / max(abs(an), np.finfo(float).eps)
        if rel_err < etol_rel:
            elapsed = time.process_time() - cpu_start
            print(
                f"Iter {n:04d}, Convergence found, stretch factor a={an1}, "
                f"relative error={rel_err}."
            )
            print(f"Time taken :{elapsed} sec.")
            return an1, n, rel_err, elapsed

        if n % progress_stride == 0:
            print(f"Iter {n:04d}, Relative error={rel_err:18.12f}.")

        an = an1

    raise RuntimeError("Newton iteration did not converge within max_iter.")


def compute_boundaries(
    xs: float, lx: float, m: int, dxs: float, nx: int, a: float
) -> tuple[np.ndarray, np.ndarray]:
    widths = np.empty(nx, dtype=float)
    widths[:m] = dxs

    dx = dxs
    for i in range(m, nx):
        dx *= a
        widths[i] = dx

    xpl = xs + np.concatenate(([0.0], np.cumsum(widths)))
    xpl[-1] = xs + lx
    xipl = xs + np.linspace(0.0, lx, nx + 1)

    if np.any(np.diff(xpl) <= 0.0):
        raise RuntimeError("Generated physical coordinates are not strictly increasing.")

    return xipl, xpl


def build_trn_lines(trn: str, trn_id: str, xipl: np.ndarray, xpl: np.ndarray) -> list[str]:
    lines = []
    for i in range(1, len(xipl) - 1):
        lines.append(
            f"&{trn} ID={trn_id}, CC={xipl[i]:20.15f}, PC={xpl[i]:20.15f} /"
        )
    return lines


def make_figure(xs: float, lx: float, xipl: np.ndarray, xpl: np.ndarray, figure_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 6))

    ax0 = axes[0]
    ax0.plot([0.1 * lx, 0.1 * lx], [xs, xs + lx], "k")
    ax0.plot([0.3 * lx, 0.3 * lx], [xs, xs + lx], "k")
    ax0.plot([0.1 * lx, 0.3 * lx], [xs, xs], "k")
    for boundary in xpl[1:]:
        ax0.plot([0.1 * lx, 0.3 * lx], [boundary, boundary], "k")
    ax0.set_ylabel("Stretching direction", fontsize=16)
    ax0.set_xticks([])
    ax0.tick_params(labelsize=14)
    ax0.set_xlim(0.0, 0.4 * lx)
    ax0.set_ylim(xs - 0.1 * lx, xs + 1.1 * lx)
    ax0.set_box_aspect(None)
    ax0.grid(False)

    ax1 = axes[1]
    ax1.plot(xipl, xpl, "+k", linewidth=2)
    ax1.set_xlabel(r"$\xi$", fontsize=16)
    ax1.set_ylabel("Physical coordinate", fontsize=16)
    ax1.tick_params(labelsize=14)
    ax1.axis("equal")
    ax1.grid(True)

    for ax in axes:
        ax.set_frame_on(True)

    fig.tight_layout()
    fig.savefig(figure_path, dpi=150)
    plt.close(fig)


def write_trn_file(output_path: Path, lines: list[str]) -> None:
    output_path.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> int:
    args = parse_args()
    validate_inputs(args.xs, args.lx, args.m, args.dxs, args.nx)

    if args.dxs >= args.lx / args.nx:
        print(
            "Warning: DXS >= Lx/Nx, so the resulting mesh may not increase in size."
        )

    output_dir = Path(args.basedir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    stretch_factor, _, _, _ = solve_stretch_factor(
        args.lx, args.m, args.dxs, args.nx, args.etol_rel, args.max_iter
    )
    xipl, xpl = compute_boundaries(
        args.xs, args.lx, args.m, args.dxs, args.nx, stretch_factor
    )

    trn_id = normalize_trn_id(args.trn_id)
    trn_lines = build_trn_lines(args.trn, trn_id, xipl, xpl)

    figure_path = output_dir / args.figure
    make_figure(args.xs, args.lx, xipl, xpl, figure_path)

    print("")
    print("Writing stretching input file...")
    output_path = output_dir / args.output
    print(output_path)
    write_trn_file(output_path, trn_lines)

    if args.print_trn:
        print("")
        for line in trn_lines:
            print(line)

    print(f"Saved figure: {figure_path}")
    print("Done.")

    if args.show:
        img = plt.imread(figure_path)
        plt.figure(figsize=(10, 6))
        plt.imshow(img)
        plt.axis("off")
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
