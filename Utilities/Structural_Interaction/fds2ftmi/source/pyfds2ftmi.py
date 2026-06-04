# -*- coding: utf-8 -*-
"""
Created on Fri Oct 17 15:15:02 2025

@author: jcsilva
"""

#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
from typing import Dict, Any, Tuple
import pyfdstools as fds
import numpy as np
import pandas as pd
import math

# -------------------------------------------------
# Parsing helpers
# -------------------------------------------------

def parse_average_token(token: str) -> Dict[str, Any]:
    """
    Accepts:
      - 'time:t0:t1' -> time window with t1>t0
      - 'avg:none' or 'average:none' -> no averaging
      - 'avg:mean:N' -> last N steps (N>=1)
    Also accepts bare 'time:t0:t1' token.
    """
    s = token.strip().lower()

    if s.startswith("time:"):
        parts = s.split(":")
        if len(parts) == 3:
            try:
                t0 = float(parts[1]); t1 = float(parts[2])
                if t1 <= t0:
                    raise ValueError
                return {"mode": "time_window", "t0": t0, "t1": t1}
            except Exception:
                raise argparse.ArgumentTypeError("time:t0:t1 requires numeric t0,t1 with t1>t0")
    
    if s.startswith("avg:") or s.startswith("average:"):
        _, _, rest = s.partition(":")
        rest = rest.strip()
        
        if rest in ("none", "no", "off", "0"):
            return {"mode": "none"}
        
        if rest.startswith("mean:"):
            try:
                n = int(rest.split(":", 1)[1])
                if n < 1:
                    raise ValueError
                return {"mode": "mean_last_n", "n": n}
            except Exception:
                raise argparse.ArgumentTypeError("avg:mean:N requires N>=1")
    
    raise argparse.ArgumentTypeError(
        "Average token must be 'time:t0:t1', 'avg:none', or 'avg:mean:N'"
    )

def parse_var_token(token: str) -> Dict[str, Any]:
    """
    Accepts:
      - 'var:1' or 'var:ast'         -> AST only, requires h
      - 'var:1:h:<value>'            -> AST only with prescribed h
      - 'var:2' or 'var:ast+h'       -> AST and h extracted
    Returns {'mode': 'AST_only'|'AST_and_h', 'h': float|None}
    """
    if not token.lower().startswith("var:"):
        raise argparse.ArgumentTypeError("Variable token must start with 'var:'")
    s = token[4:].strip().lower()

    if s in ("1", "ast"):
        return {"mode": "AST_only", "h": None}
    if s in ("2", "ast+h", "ast_h", "both"):
        return {"mode": "AST_and_h", "h": None}
    if s.startswith("1:h:"):
        _, _, hv = s.partition("1:h:")
        try:
            h = float(hv)
            if h <= 0:
                raise ValueError
        except Exception:
            raise argparse.ArgumentTypeError("For 'var:1:h:<value>', <value> must be positive float (W/m^2/K)")
        return {"mode": "AST_only", "h": h}
    raise argparse.ArgumentTypeError("var must be '1', '2', 'ast', 'ast+h', or '1:h:<value>'")

def parse_cell_size_token(token: str) -> float:
    """
    Accepts:
      - 'cell:<value>' -> positive float cell size
    """
    if not token.lower().startswith("cell:"):
        raise argparse.ArgumentTypeError("Cell size token must start with 'cell:'")
    s = token.split(":", 1)[1].strip()
    try:
        val = float(s)
        if val <= 0:
            raise ValueError
        return val
    except Exception:
        raise argparse.ArgumentTypeError("cell:<value> requires a positive float")

def parse_key_value_token(prefix: str, token: str) -> str:
    if not token.lower().startswith(prefix + ":"):
        raise argparse.ArgumentTypeError(f"Token must start with '{prefix}:'")
    return token.split(":", 1)[1]

def parse_bool_token(val: str) -> bool:
    s = val.strip().lower()
    if s in ("1", "true", "yes", "on"): return True
    if s in ("0", "false", "no", "off"): return False
    raise argparse.ArgumentTypeError("Boolean must be one of: true/false, yes/no, on/off, 1/0")

def parse_bbox_token(val: str) -> Dict[str, float]:
    # Accept "xmin,xmax,ymin,ymax,zmin,zmax"
    parts = [p.strip() for p in val.split(",")]
    if len(parts) != 6:
        raise argparse.ArgumentTypeError("bbox must be 6 comma-separated floats: xmin,xmax,ymin,ymax,zmin,zmax")
    try:
        xmin, xmax, ymin, ymax, zmin, zmax = map(float, parts)
    except Exception:
        raise argparse.ArgumentTypeError("bbox values must be floats")
    if not (xmin <= xmax and ymin <= ymax and zmin <= zmax):
        raise argparse.ArgumentTypeError("bbox must satisfy xmin<=xmax, ymin<=ymax, zmin<=zmax")
    return {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax, "zmin": zmin, "zmax": zmax}

def parse_choice(val: str, choices: Tuple[str, ...], name: str) -> str:
    s = val.strip().lower()
    if s not in choices:
        raise argparse.ArgumentTypeError(f"{name} must be one of: {', '.join(choices)}")
    return s

# -------------------------------------------------
# CLI that accepts flexible flag-like tokens
# -------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fds2ftmi_py",
        description="Process FDS results and export for FTMI/FE codes using in-memory access."
    )
    p.add_argument(
        "tokens",
        nargs="*",
        help=(
            "Flag-style tokens in any order. Required: chid:, var:, code:, file:, cell:. "
            "Optional: time: or avg:. "
            "Example: chid:caseA var:1:h:35 code:ftmi file:export.csv cell:0.05 time:120:360"
        ),
    )
    # Optional flags (override tokens if given)
    p.add_argument("--chid", type=str, help="FDS CHID (required).")
    p.add_argument("--var", type=str, help="Variables mode like '1', '2', or '1:h:35' (required).")
    p.add_argument("--time", type=str, help="Time window 't0:t1' (optional).")
    p.add_argument("--avg", type=str, help="Average spec: 'none' | <seconds> | 'mean:N' (optional).")
    p.add_argument("--code", type=str, choices=["ftmi", "abaqus", "ansys", "opensees", "custom"], help="Target code (required).")
    p.add_argument("--file", type=str, help="Output/input file path (required).")
    p.add_argument("--cell", type=float, help="Cell size, positive float (required, e.g., 0.05).")
    p.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity.")
    p.add_argument("--root", type=Path, default=Path("."), help="Root directory for FDS outputs.")
    
    return p

def compile_from_tokens(args) -> Dict[str, Any]:
    # Defaults
    out: Dict[str, Any] = {
        "chid": None,
        "variables": None,   # dict with mode,h
        "average": {"mode": "none"},  # default: no averaging (optional)
        "code": None,
        "file": None,
        "cell_size": None,
        "root": args.root.resolve(),
        "verbose": int(args.verbose),
    }

    # Tokens (optional support)
    for tok in args.tokens or []:
        low = tok.lower()
        if low.startswith("chid:"):
            out["chid"] = parse_key_value_token("chid", tok)
        elif low.startswith("var:"):
            out["variables"] = parse_var_token(tok)
        elif low.startswith("time:") or low.startswith("avg:") or low.startswith("average:"):
            out["average"] = parse_average_token(tok)
        elif low.startswith("code:"):
            out["code"] = parse_key_value_token("code", tok).lower()
        elif low.startswith("file:"):
            out["file"] = parse_key_value_token("file", tok)
        elif low.startswith("cell:"):
            out["cell_size"] = parse_cell_size_token(tok)
        else:
            raise argparse.ArgumentTypeError(
                f"Unrecognized token '{tok}'. Expected chid:, var:, time:, avg:, code:, file:, cell: "
            )

    # Overrides via flags
    if args.chid is not None:
        out["chid"] = args.chid
    if args.var is not None:
        out["variables"] = parse_var_token(f"var:{args.var}")
    if args.time is not None:
        out["average"] = parse_average_token(f"time:{args.time}")
    if args.avg is not None:
        out["average"] = parse_average_token(f"avg:{args.avg}")
    if args.code is not None:
        out["code"] = args.code.lower()
    if args.file is not None:
        out["file"] = args.file
    if args.cell is not None:
        if args.cell <= 0:
            raise argparse.ArgumentTypeError("--cell must be a positive float")
        out["cell_size"] = args.cell

    # Required keys (unchanged)
    missing = [k for k in ("chid", "variables", "code", "file", "cell_size") if out[k] in (None, "")]
    if missing:
        raise argparse.ArgumentTypeError(
            "Missing required tokens/flags: " + ", ".join(missing) + ". "
            "Provide: chid:caseA var:1:h:35 code:ftmi file:export.csv cell:0.05 "
            "Optionally add time:t0:t1 or avg:none."
        )

    if out["variables"]["mode"] == "AST_only" and out["variables"]["h"] is None:
        raise argparse.ArgumentTypeError("AST-only mode requires prescribing h via 'var:1:h:<value>'.")

    # Normalize file paths
    file_path = Path(out["file"])
    out["file"] = (file_path if file_path.is_absolute() else (out["root"] / file_path)).resolve()

    return out

# -------------------------------------------------
# Public API
# -------------------------------------------------

def get_config(argv=None) -> Dict[str, Any]:
    """
    Parse CLI args/tokens and return a config dict for further processing or tests.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        cfg = compile_from_tokens(args)
    except argparse.ArgumentTypeError as e:
        parser.error(str(e))
    return cfg

# -------------------------------------------------
# Main entry point
# -------------------------------------------------

def main(argv=None):
    """
    CLI entry: returns (exit_code, config_dict).
    """
    try:
        cfg = get_config(argv)
    except SystemExit as e:
        return e.code, None

    if cfg["verbose"]:
        print(f"CHID        : {cfg['chid']}")
        print(f"Variables   : {cfg['variables']}")
        print(f"Averaging   : {cfg['average']}")
        print(f"Code        : {cfg['code']}")
        print(f"File        : {cfg['file']}")
        print(f"Cell size   : {cfg['cell_size']}")
        print(f"Root        : {cfg['root']}")

    # Return cfg so downstream code can use it immediately
    return 0, cfg


###########################

def compute_face_centroids(vertices, faces):
    """
    Compute centroid coordinates for each face
    
    Parameters:
    vertices: np.array shape (n_vertices, 3)
    faces: np.array shape (n_faces, 3) - triangular face indices
    
    Returns:
    face_centroids: np.array shape (n_faces, 3)
    """
    face_centroids = np.zeros((faces.shape[0], 3))
    faces=faces-1
    for i, face in enumerate(faces):
        # Get vertices of this face
        face_vertices = vertices[face]
        # Compute centroid (mean of vertex coordinates)
        centroid = face_vertices.mean(axis=0)
        face_centroids[i] = centroid
    return face_centroids

def infer_gcf_from_be(be_file_path):
    """
    Infer gcf filename from .be filename using FDS naming convention:
    
    Pattern: CHID_mesh_variable.be -> CHID_mesh.gcf
    Examples:
      h_profile_geom_1_1.be -> h_profile_geom_1.gcf
      h_profile_geom_1_2.be -> h_profile_geom_1.gcf
      h_profile_geom_3_5.be -> h_profile_geom_3.gcf
    
    Parameters:
    be_file_path: str or Path - path to .be file
    
    Returns:
    str - inferred path to .gcf file
    """
    be_path = Path(be_file_path)
    stem = be_path.stem  # filename without extension
    
    # Split by '_' from right to find last underscore numbering (variable number)
    parts = stem.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        base_stem = parts[0]  # remove last variable number
    else:
        base_stem = stem  # no variable number found, use as-is
    
    gcf_name = base_stem + '.gcf'
    gcf_path = be_path.parent / gcf_name
    return str(gcf_path)


def read_nodes_and_elements(nodes_file='nodes.dat', elements_file='elements.dat'):
    """
    Read nodes and elements from input files.
    
    Parameters:
    -----------
    nodes_file : str
        Path to nodes data file
    elements_file : str
        Path to elements data file
        
    Returns:
    --------
    dict : Dictionary containing all processed data
    """
    
    # Read nodes file
    nodes_path = Path(nodes_file)
    
    # Read header: number of nodes and highest node number
    with open(nodes_path, 'r') as f:
        first_line = f.readline().split()
        nnode = int(first_line[0])
        highnode = int(first_line[1])
    
    # Read node data: node_id, x, y, z
    data = np.loadtxt(nodes_path, skiprows=1)
    
    if config['code']=='ftmi':
        nos = data[:, :5]  # All 4 columns
    else:
        nos = data[:, :4]  # All 3 columns
    
    # Read elements file with multiple groups
    element_groups = read_element_groups(elements_file)
    
    # Combine all elements into a single array
    all_elements = []
    element_metadata = []
    
    for group in element_groups:
        all_elements.append(group['elements'])
        element_metadata.append({
            'shell': group['shell'],
            'layer': group['layer'],
            'el': group['el'],
            'numel': group['numel'],
            'start_index': len(all_elements) - 1
        })
    
    # Stack all elements together
    elements = np.vstack(all_elements)
    total_numel = len(elements)
    
    # Initialize arrays
    nomaster = np.zeros((total_numel, 4))
    n = np.zeros((total_numel, 4))
    
    # Normal orientation (1 for outside, -1 for inside)
    n_or = 1
    
    # Process each element to create master nodes and normals
    nomaster, n = create_nomaster_nodes(
        elements, nos, nnode, highnode, n_or, total_numel
    )
    
    return {
        'nnode': nnode,
        'highnode': highnode,
        'nos': nos,
        'numel': total_numel,
        'element_groups': element_metadata,
        'elements': elements,
        'nomaster': nomaster,
        'n': n,
        'n_or': n_or
    }


def read_element_groups(elements_file):
    """
    Read element groups from elements file.
    Each group has a header line (numel, shell, layer, el) followed by element data.
    Groups are separated by "END" keyword.
    """
    
    element_groups = []
    elements_path = Path(elements_file)
    
    with open(elements_path, 'r') as f:
        while True:
            # Read next line
            line = f.readline()
            
            if not line:  # End of file
                break
            
            line = line.strip()
            
            if not line or line.upper() == 'END':  # Empty line or END marker
                continue
            
            # Parse header: numel, shell, layer, el
            header_values = line.split()
            
            # Validate this is a header (should have exactly 4 values)
            if len(header_values) != 4:
                continue
            
            numel = int(header_values[0])
            shell = int(header_values[1])
            layer = int(header_values[2])
            el = int(header_values[3])
            
            # Read exactly numel element lines
            elements_data = []
            for i in range(numel):
                element_line = f.readline().strip()
                
                if not element_line or element_line.upper() == 'END':
                    break
                
                # Parse element data (9 columns)
                element_values = [float(x) for x in element_line.split()]
                
                # Validate this is an element line (should have 9 values)
                if len(element_values) == 9:
                    elements_data.append(element_values)
            
            # Only add group if we read elements
            if elements_data:
                elements_array = np.array(elements_data)
                
                element_groups.append({
                    'numel': len(elements_data),  # Use actual count
                    'shell': shell,
                    'layer': layer,
                    'el': el,
                    'elements': elements_array
                })
    
    return element_groups


def create_nomaster_nodes(elements, nos, nnode, highnode, n_or, numel):
    """
    Create master nodes at element centroids and calculate normal vectors.
    
    Parameters:
    -----------
    elements : ndarray
        Element connectivity array (numel x 9)
    nos : ndarray
        Node coordinates array (nnode x 4)
    nnode : int
        Number of nodes
    highnode : int
        Highest node number
    n_or : int
        Normal orientation (+1 or -1)
    numel : int
        Number of elements
        
    Returns:
    --------
    nomaster : ndarray
        Master node data (numel x 4): [node_id, x, y, z]
    n : ndarray
        Normal vectors (numel x 4): [node_id, nx, ny, nz]
    """
    
    nomaster = np.zeros((numel, 4))
    n = np.zeros((numel, 4))
    
    # Create a dictionary for fast node lookup: {node_id: [x, y, z]}
    node_dict = {int(nos[i, 0]): nos[i, 1:4] for i in range(nnode)}
    
    for i in range(numel):
        highnode += 1
        
        # Get node IDs for this element
        no1 = int(elements[i, 1]) 
        no2 = int(elements[i, 2]) 
        no3 = int(elements[i, 3]) 
        no4 = int(elements[i, 4]) 
        
        # Get coordinates for nodes
        coord1 = node_dict[no1]
        coord2 = node_dict[no2]
        coord3 = node_dict[no3]
        
        # Store coordinates for normal calculation
        a = coord1
        b = coord2
        c = coord3
        
        # Calculate centroid
        if no4 == no3:
            # Triangle element (3 nodes)
            centroid = (coord1 + coord2 + coord3) / 3.0
        else:
            # Quadrilateral element (4 nodes)
            coord4 = node_dict[no4]
            centroid = (coord1 + coord2 + coord3 + coord4) / 4.0
        
        # Store master node data
        nomaster[i, 0] = highnode
        nomaster[i, 1:4] = centroid
        
        # Calculate normal vector using cross product
        if config['code']=='ftmi':
            if nos[i,5]==1:
                n[i, 0] = highnode
                n[i, 1] = 1
                n[i, 2] = 0
                n[i, 3] = 0
            elif nos[i,5]==-1:
                n[i, 0] = highnode
                n[i, 1] = -1
                n[i, 2] = 0
                n[i, 3] = 0                
            elif nos[i,5]==2:
                n[i, 0] = highnode
                n[i, 1] = 0
                n[i, 2] = 1
                n[i, 3] = 0
            elif nos[i,5]==-2:
                n[i, 0] = highnode
                n[i, 1] = 0
                n[i, 2] = -1
                n[i, 3] = 0
            elif nos[i,5]==3:
                n[i, 0] = highnode
                n[i, 1] = 0
                n[i, 2] = 0
                n[i, 3] = 1
            elif nos[i,5]==-3:
                n[i, 0] = highnode
                n[i, 1] = 0
                n[i, 2] = 0
                n[i, 3] = -1                
        else:
            # Normal = (B-A) × (C-A)
            vec_ab = b - a
            vec_ac = c - a
            
            normal = np.cross(vec_ab, vec_ac)
            
            # Apply orientation
            normal = n_or * normal
                                    
        # Store normal vector
        n[i, 0] = highnode
        n[i, 1:4] = normal
    
    return nomaster, n

def average_results_by_orientation(data_list):
    """
    Average results from multiple dictionaries considering their orientations.
    
    For surfaces with the same orientation: simple average
    For surfaces with different orientations: vector composition method (FTMI approach)
    
    Parameters:
    -----------
    data_list : list of dict
        Each dictionary contains:
        - 'time': numpy array of time values
        - 'result': numpy array of result values (e.g., TAST, heat flux)
        - 'orientation': integer representing orientation
          (1=+x, -1=-x, 2=+y, -2=-y, 3=+z, -3=-z)
    
    Returns:
    --------
    time_avg : numpy array
        Time array (same as input)
    result_avg : numpy array
        Averaged result values
    
    Reference:
    ----------
    Silva et al. (2016), Fire Safety Journal 83:66-78
    Equations 10-13 for vector composition method
    """
    
    if not data_list:
        raise ValueError("Empty data list provided")
    
    # Get time array (should be same for all)
    time_avg = data_list[0]['times']
    
    # Check if all dictionaries have the same orientation
    orientations = [d['orientation'] for d in data_list]
    unique_orientations = set(orientations)
    
    if len(unique_orientations) == 1:
        # Case 1: All same orientation - simple average
        result_sum = np.zeros_like(data_list[0]['time_series'], dtype=float)
        
        for data in data_list:
            result_sum += data['time_series']
        
        result_avg = result_sum / len(data_list)
        
    else:
        # Case 2: Different orientations - vector composition method
        result_avg = vector_composition_average(data_list, time_avg)
    
    return time_avg, result_avg


def vector_composition_average(data_list, time_array):
    """
    Average results using vector composition for different orientations.
    
    Based on FTMI methodology (Silva et al., 2016):
    - Create vectors for each orientation contribution
    - Project onto the target normal direction
    
    Parameters:
    -----------
    data_list : list of dict
        List of dictionaries with 'time', 'result', 'orientation'
    time_array : numpy array
        Common time array
    
    Returns:
    --------
    result_avg : numpy array
        Averaged result using vector composition
    """
    
    # Initialize contribution vectors for each direction
    n_timesteps = len(time_array)
    contributions_x = np.zeros(n_timesteps, dtype=float)
    contributions_y = np.zeros(n_timesteps, dtype=float)
    contributions_z = np.zeros(n_timesteps, dtype=float)
    
    # Orientation convention:
    # 1 = +x direction, -1 = -x direction
    # 2 = +y direction, -2 = -y direction
    # 3 = +z direction, -3 = -z direction
    
    orientation_vectors = {
        1: np.array([1.0, 0.0, 0.0]),    # +x
        -1: np.array([-1.0, 0.0, 0.0]),  # -x
        2: np.array([0.0, 1.0, 0.0]),    # +y
        -2: np.array([0.0, -1.0, 0.0]),  # -y
        3: np.array([0.0, 0.0, 1.0]),    # +z
        -3: np.array([0.0, 0.0, -1.0])   # -z
    }
    
    # Accumulate contributions from each orientation
    for data in data_list:
        orientation = data['orientation']
        result = data['time_series']
        
        # Get unit vector for this orientation
        unit_vec = orientation_vectors[orientation]
        
        # Add contribution (Eq. 10 in Silva et al., 2016)
        contributions_x += result * unit_vec[0]
        contributions_y += result * unit_vec[1]
        contributions_z += result * unit_vec[2]
    
    # Create result vector
    v_result = np.array([contributions_x, contributions_y, contributions_z])
    
    # Calculate target normal direction (average of all normal directions)
    normal_sum = np.zeros(3, dtype=float)
    for data in data_list:
        normal_sum += orientation_vectors[data['orientation']]
    
    # Handle case where normals cancel out (opposite orientations)
    norm_magnitude = np.linalg.norm(normal_sum)
    if norm_magnitude < 1e-10:
        # Opposite orientations - use magnitude of result vector
        result_avg = np.linalg.norm(v_result, axis=0)
    else:
        normal_avg = normal_sum / norm_magnitude
        
        # Project result vector onto average normal direction (Eq. 12-13 in Silva et al., 2016)
        # result_avg = |v_result · normal_avg|
        result_avg = np.abs(np.dot(normal_avg, v_result))
    
    return result_avg

def average_multiple_results(data_list):
    """
    Simple average of results from multiple dictionaries.
    Used when multiple FDS results are found within the same FEM cell_size (for GEOMs).
    
    Parameters:
    -----------
    data_list : list of dict
        Each dictionary contains:
        - 'time': numpy array of time values
        - 'result': numpy array of result values (e.g., TAST, heat flux)
    
    Returns:
    --------
    time_avg : numpy array
        Time array (same as input)
    result_avg : numpy array
        Averaged result values (simple arithmetic mean)
    """
    
    if not data_list:
        raise ValueError("Empty data list provided")
    
    if len(data_list) == 1:
        # Only one result - return as is
        return data_list[0]['times'], data_list[0]['time_series']
    
    # Get time array (should be same for all)
    time_avg = data_list[0]['times']
    
    # Simple arithmetic average
    result_sum = np.zeros_like(data_list[0]['time_series'], dtype=float)
    
    for data in data_list:
        result_sum += data['time_series']
    
    result_avg = result_sum / len(data_list)
    
    return time_avg, result_avg

def apply_time_window_average(time_array, values_array, start_time, end_time, use_fast=True):
    """
    Calculate the average of a time series within a time window
    and prescribe this average for the remaining duration.
    
    Parameters:
    -----------
    time_array : numpy.ndarray
        Array containing time values
    values_array : numpy.ndarray
        Array containing time series values
    start_time : float
        Start of the averaging time window
    end_time : float
        End of the averaging time window
    use_fast : bool
        If True, uses searchsorted (faster for large arrays)
        If False, uses boolean masking
    
    Returns:
    --------
    dict : Dictionary with modified time and values arrays
           Original time array is preserved, values are set to 
           the windowed average from end_time onwards
    """
   
    if use_fast:
        # Fast method using searchsorted
        start_idx = np.searchsorted(time_array, start_time, side='left')
        end_idx = np.searchsorted(time_array, end_time, side='right')
    else:
        # Alternative method using boolean masking
        mask = (time_array >= start_time) & (time_array <= end_time)
        indices = np.where(mask)[0]
        if len(indices) == 0:
            raise ValueError("No time values found within the specified window")
        start_idx = indices[0]
        end_idx = indices[-1] + 1
    
    # Validate the window
    if start_idx >= end_idx:
        raise ValueError(f"No valid data points in time window [{start_time}, {end_time}]")
    
    # Extract values within the window and calculate average
    windowed_values = values_array[start_idx:end_idx]
    window_average = windowed_values.mean()
    
    # Find where to start prescribing the average (from end_time onwards)
    prescription_idx = end_idx
    
    # Set all values from prescription point onwards to the window average
    values_array[prescription_idx:] = window_average
    
    # Return modified data with original time but updated values
    return  window_average

def read_nodes_ftmi(nodes_file='nodes.dat'):
    """
    Read nodes from FTMI format input file.
    
    Format:
    - Line 1: number_of_nodes
    - Lines 2+: ID, x, y, z, orientation (CSV)
    
    Parameters:
    -----------
    nodes_file : str
        Path to nodes data file
        
    Returns:
    --------
    dict : Dictionary containing:
        - 'nnode': int - number of nodes
        - 'node_ids': list of str - position IDs
        - 'coordinates': np.array - (nnode, 3) with x, y, z
        - 'orientations': np.array - (nnode,) with orientation values
    """
    nodes_path = Path(nodes_file)
    
    # Read header: number of nodes
    with open(nodes_path, 'r') as f:
        first_line = f.readline().strip()
        nnode = int(first_line)
    
    # Read CSV data: ID, x, y, z, orientation
    data = pd.read_csv(nodes_path, skiprows=1, names=['ID', 'x', 'y', 'z', 'orientation'])
    
    return {
        'nnode': nnode,
        'node_ids': data['ID'].tolist(),
        'coordinates': data[['x', 'y', 'z']].to_numpy(),
        'orientations': data['orientation'].to_numpy()
    }


def write_ftmi_csv(node_data, config):
    """
    Write node data to CSV file in DEVC-like format.
    Groups nodes by time vectors to handle different time discretizations.
    
    Parameters:
    -----------
    node_data : dict
        Dictionary with node_id as key and dict of time, AST, h as values
    config : dict
        Configuration dictionary with file path and variable mode
    """
    if not node_data:
        print("Warning: No node data to write")
        return
    
    # Group nodes by identical time vectors
    time_groups = {}
    
    for node_id, data in node_data.items():
        if 'time' not in data or data['time'] is None:
            continue
            
        # Use hash of time array as key (round to avoid floating point issues)
        time_tuple = tuple(np.round(data['time'], 6))
        
        if time_tuple not in time_groups:
            time_groups[time_tuple] = {
                'time': data['time'],
                'nodes': []
            }
        
        time_groups[time_tuple]['nodes'].append({
            'id': node_id,
            'AST': data.get('AST'),
            'h': data.get('h')
        })
    
    if not time_groups:
        print("Warning: No valid time data found")
        return
    
    print(f"\nOrganized {len(node_data)} nodes into {len(time_groups)} time group(s)")
    
    # Build column structure
    columns = []
    data_arrays = []
    
    for group_idx, (time_tuple, group) in enumerate(sorted(time_groups.items(), 
                                                            key=lambda x: len(x[1]['nodes']), 
                                                            reverse=True)):
        time_array = group['time']
        
        # Add time column (with suffix if multiple groups)
        if len(time_groups) == 1:
            time_col = 'Time'
        else:
            time_col = f'Time_{group_idx+1}'
        
        columns.append(time_col)
        data_arrays.append(time_array)
        
        print(f"  Group {group_idx+1}: {len(group['nodes'])} nodes, {len(time_array)} time steps")
        
        # Add node data columns for this time group
        for node in sorted(group['nodes'], key=lambda x: x['id']):
            if node['AST'] is not None:
                columns.append(f"{node['id']}_AST")
                data_arrays.append(node['AST'])
            
            if node['h'] is not None:
                columns.append(f"{node['id']}_h")
                data_arrays.append(node['h'])
    
    # Find maximum length for padding
    max_len = max(len(arr) for arr in data_arrays)
    
    # Pad arrays with NaN for unequal lengths
    data_dict = {}
    for col, arr in zip(columns, data_arrays):
        if len(arr) < max_len:
            padded = np.full(max_len, np.nan)
            padded[:len(arr)] = arr
            data_dict[col] = padded
        else:
            data_dict[col] = arr
    
    # Create DataFrame
    df = pd.DataFrame(data_dict)
    
    # Write to CSV
    df.to_csv(config['file'], index=False, float_format='%.6e')
    
    print(f"\nCSV file written to: {config['file']}")
    print(f"  Columns: {len(columns)}")
    print(f"  Rows: {len(df)}")

def build_cell_centers(patch):
    """
    Calculate cell center coordinates for a BNDF patch.
    FDS stores data at cell centers, which are offset from grid nodes.
    
    Parameters:
    -----------
    patch : BNDF patch object
        Must have x, y, z grid coordinates already built
    
    Returns:
    --------
    None (modifies patch object in-place, adding x_centers, y_centers, z_centers)
    """
    if not hasattr(patch, 'x'):
        patch.buildSpace()
    
    # Calculate cell centers by averaging adjacent grid nodes
    # If grid is (n+1, m+1), centers are (n, m)
    patch.x_centers = (patch.x[:-1, :-1] + patch.x[1:, 1:]) / 2.0
    patch.y_centers = (patch.y[:-1, :-1] + patch.y[1:, 1:]) / 2.0
    patch.z_centers = (patch.z[:-1, :-1] + patch.z[1:, 1:]) / 2.0
    
    
    
#######################################################################    

if __name__ == "__main__":
    code, config = main()
    if isinstance(config, dict):
        print({"status": "ok", **config})
    else:
        sys.exit(code)

    case_path = Path(Path.cwd() / f"{config['chid']}.fds")
    if case_path.is_dir():
        # Find .smv file in directory
        smv_files = list(case_path.glob('*.smv'))
        if not smv_files:
            raise ValueError(f"No .smv file found in {case_path}")
        smv_file = smv_files[0]
        chid = smv_file.stem
        result_dir = case_path
    else:
        # Assume it's a case name, look for files
        if case_path.suffix == '.fds':
            chid = case_path.stem
            result_dir = case_path.parent
        else:
            chid = str(case_path)
            result_dir = Path('.')
        smv_file = result_dir / f"{chid}.smv"
        if not smv_file.exists():
            raise ValueError(f"SMV file not found: {smv_file}")

    # Parse SMV file to get mesh and boundary info
    smv_data = fds.smokeviewParser.parseSMVFile(str(smv_file))

    # Get list of boundary files for the quantity
    bf_files = fds.utilities.getFileList(str(result_dir), chid, 'bf')
    
    # Get list of .be files (geometry boundary data files)
    be_files = fds.utilities.getFileList(str(result_dir), chid, 'be')
    
    results_bf=[]
    results_be=[]
    
    if config['variables']['mode']=='AST_and_h':
        quantities=['ADIABATIC SURFACE TEMPERATURE','HEAT TRANSFER COEFFICIENT']
    elif config['variables']['mode']=="AST_only": 
        quantities=['ADIABATIC SURFACE TEMPERATURE']
    else:
        print("ERROR: Variables not defined correctly")
        sys.exit()                
        
    for quantity in quantities:
        
        for bf_file in bf_files:
        
            try:
                # Check if this file contains our quantity
                file_quantity, short_name, units, npatch = fds.extractBoundaryData.readBoundaryHeader(bf_file)
                if file_quantity != quantity:
                    continue
    
                # Determine mesh number from filename
                filename = Path(bf_file).name
                try:
                    mesh_num = int(filename.split('_')[-2].split('.')[0]) -1  # Convert to 0-based
                except:
                    mesh_num = 0  # Default to first mesh if can't parse
    
                print(f"Processing {filename} for quantity '{file_quantity}' (mesh {mesh_num})")
    
                # Import boundary data from this file
                times, patches, units = fds.extractBoundaryData.importBoundaryFile(
                    bf_file, str(smv_file), gridNum=mesh_num, smvData=smv_data)
                
                if patches is not None:
                    print(f"  Precomputing cell centers for {len(patches)} patches...")
                    for patch in patches:
                        try:
                            build_cell_centers(patch)
                        except Exception as e:
                            print(f"  Warning: Could not build cell centers for a patch: {e}")
            
                if patches is None:
                    continue
    
                results_bf.append({
                        'source_type': 'BNDF',
                        'bndf_file': bf_file,
                        'mesh_id': mesh_num,
                        'times': times,
                        'quantity': file_quantity,
                        'patches': patches,
                    })
    
            except Exception as e:
                print(f"Warning: Error processing {bf_file}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
   
    # Get geometry boundary file quantities
    for be_file in be_files:
        try:
            file_quantities = fds.extractGeomData.getBndeQuantities(str(smv_file))
            print(f"Found {len(file_quantities)} GEOM boundary file entries in SMV")
            filename = Path(be_file).name
            print(f"Processing: {filename}")
            
            # Check if this file is in the SMV file and has our quantity
            if filename not in file_quantities:
                print(f"  Skipping {filename} - not found in SMV file")
                continue
                
            file_quantity = file_quantities[filename]['quantity']
            if file_quantity not in quantities:
                print(f"  Skipping {filename} - quantity '{file_quantity}' != '{quantity}'")
                continue
            
            print(f"  Processing {filename} with quantity '{file_quantity}'")
            
            # Infer corresponding geometry file (.gcf) using naming convention
            gcf_file = infer_gcf_from_be(be_file)
            print(f"  Inferred .gcf file: {Path(gcf_file).name}")
            
            if not Path(gcf_file).exists():
                print(f"  Error: Inferred geometry file not found: {gcf_file}")
                continue
            
            # Read geometry data
            vertices, faces, geom_header = fds.extractGeomData.readGcfFile(gcf_file)
            print(f"  Loaded geometry: {len(vertices)} vertices, {len(faces)} faces")
            nfaces_fds = len(faces)
            nvertices_fds = len(vertices)
            
            # Read boundary data  
            times, values, data_header = fds.extractGeomData.readBeFile(be_file)
            print(f"  Loaded boundary data: {len(times)} time steps, {values.shape[0]} data points")
            
            # Verify that values correspond to faces
            if values.shape[0] != len(faces):
                print(f"  Warning: Data size {values.shape[0]} doesn't match face count {len(faces)}")
                print(f"    This might indicate vertex-based data or a different data structure")
                # Continue processing anyway - might still work
            
            # Compute face centroids
            face_centroids = compute_face_centroids(vertices, faces)
            print(f"  Computed {len(face_centroids)} face centroids")
            
            results_be.append({
                'source_type': 'GEOM',
                'geom_file': gcf_file,
                'data_file': be_file,
                'mesh_id': file_quantities[filename]['mesh'],
                'times': times,
                'quantity': file_quantity,
                'n_vertices': len(vertices),
                'n_faces': len(faces),
                'all_face_centroids': face_centroids,  # Include all centroids for reference
                'all_vertices': vertices,              # Include all vertices for reference
                'all_faces': faces,                    # Include all faces for reference
                'values': values
            })
                
        except Exception as e:
            print(f"  Error processing GEOM file {be_file}: {e}")
            import traceback
            traceback.print_exc()
            continue
        
        if not file_quantities:
            print("No GEOM boundary quantities found in SMV file")
            continue #return []
        
    
#########################################################    
    match config['code']:
        case 'ansys':
            nodes_elements = read_nodes_and_elements()
            
            print(f"Number of nodes: {nodes_elements['nnode']}")
            print(f"Total number of elements: {nodes_elements['numel']}")
            print(f"\nElement groups:")
            for i, group in enumerate(nodes_elements['element_groups']):
                print(f"  Group {i+1}: numel={group['numel']}, shell={group['shell']}, "
                      f"layer={group['layer']}, el={group['el']}")
                
            element_start_index=0
            nomaster_start_index=0
            normal_start_index=0
            with config['file'].open('w') as file70, open(config['file'].parent / f"{config['file'].stem}_loads.dat", 'w') as file71:
                for group in nodes_elements['element_groups']:
                    el=group['el']
                    shell=group['shell']
                    layer=group['layer']
                    numel=group['numel']
                    nomaster = nodes_elements['nomaster'][0+nomaster_start_index:numel+nomaster_start_index]
                    nomaster_start_index+=numel
                    normal = nodes_elements['n'][0+normal_start_index:numel+normal_start_index]
                    normal_start_index+=numel
                    file70.write("/PREP7\n")
                    for i in range(numel):
                        if shell==1:
                            n_or=1
                            vec_mult=config['cell_size']/2/(math.sqrt((normal[i,1]**2)+(normal[i,2]**2)+(normal[i,3]**2)))
                            nomaster[i,1]+=n_or*normal[i,1]*vec_mult
                            nomaster[i,2]+=n_or*normal[i,2]*vec_mult
                            nomaster[i,3]+=n_or*normal[i,3]*vec_mult
                        
                        # Format: N, integer(8 width), float(7.3), float(7.3), float(7.3), trailing commas
                        file70.write(f"N,{int(nomaster[i, 0]):8d},"
                                     f"{nomaster[i, 1]:7.3f},"
                                     f"{nomaster[i, 2]:7.3f},"
                                     f"{nomaster[i, 3]:7.3f},,,,,\n")
                    file70.write("!*\n")
                    file70.write(f"ET,{el:8d},SURF152\n")
                    file70.write("!*\n")
                    file70.write(f"KEYOPT,{el:8d},1,0\n")
                    file70.write(f"KEYOPT,{el:8d},2,0\n")
                    file70.write(f"KEYOPT,{el:8d},3,0\n")
                    file70.write(f"KEYOPT,{el:8d},4,0\n")
                    file70.write(f"KEYOPT,{el:8d},5,1\n")
                    file70.write(f"KEYOPT,{el:8d},6,0\n")
                    file70.write(f"KEYOPT,{el:8d},7,0\n")
                    file70.write(f"KEYOPT,{el:8d},8,4\n")
                    file70.write(f"KEYOPT,{el:8d},9,1\n")
                    
                    # IF IT IS SHELL ELEMENTS,
                    # 1 FOR TOP LAYER AND 2 FOR BOTTOM LAYER
                    if shell == 1:
                        if layer == 1:
                            file70.write(f"KEYOPT,{el:8d},11,1\n")
                        if layer == 2:
                            file70.write(f"KEYOPT,{el:8d},11,2\n")
                    
                    file70.write("!*\n")
                    file70.write("!*\n")
                    file70.write(f"R,{el:8d},1,5.67E-08, , , ,\n")
                    file70.write("RMORE, , , ,\n")
                    file70.write("RMORE, , ,\n")
                    file70.write("!*\n")
                        
                    #************* CREATE SURF152 ELEMENTS ********
                    print('LOOP_SURF152')
                    
                    file70.write("!*\n")
                    file70.write(f"TYPE,{el:8d}\n")
                    file70.write(f"MAT,{1:8d}\n")
                    file70.write(f"REAL,{el:8d}\n")
                    file70.write(f"ESYS,{0:8d}\n")
                    file70.write("SECNUM,\n")
                    file70.write("TSHAP,LINE\n")
                    file70.write("!*\n")
                    
                    # LOOP_SURF152
                    elementos = nodes_elements['elements'][0+element_start_index:numel+element_start_index]
                    element_start_index+=numel
                    
                    
                    for i in range(numel):
                        if elementos[i, 5] == 0:  # Column 6 in FORTRAN (0-indexed in Python)
                            file70.write(f"nsel,S,node,,{int(elementos[i, 1]):8d}\n")  # Column 2
                            file70.write(f"nsel,A,node,,{int(elementos[i, 2]):8d}\n")  # Column 3
                            file70.write(f"nsel,A,node,,{int(elementos[i, 3]):8d}\n")  # Column 4
                            file70.write(f"nsel,A,node,,{int(elementos[i, 4]):8d}\n")  # Column 5
                            file70.write(f"ESURF,{int(nomaster[i, 0]):8d}\n")
                            file70.write("!*\n")
                        else:
                            file70.write(f"nsel,S,node,,{int(elementos[i, 1]):8d}\n")  # Column 2
                            file70.write(f"nsel,A,node,,{int(elementos[i, 2]):8d}\n")  # Column 3
                            file70.write(f"nsel,A,node,,{int(elementos[i, 3]):8d}\n")  # Column 4
                            file70.write(f"nsel,A,node,,{int(elementos[i, 4]):8d}\n")  # Column 5
                            file70.write(f"nsel,A,node,,{int(elementos[i, 5]):8d}\n")  # Column 6
                            file70.write(f"nsel,A,node,,{int(elementos[i, 6]):8d}\n")  # Column 7
                            file70.write(f"nsel,A,node,,{int(elementos[i, 7]):8d}\n")  # Column 8
                            file70.write(f"nsel,A,node,,{int(elementos[i, 8]):8d}\n")  # Column 9
                            file70.write(f"ESURF,{int(nomaster[i, 0]):8d}\n")
                            file70.write("!*\n")
                    
                    file70.write("ALLSEL,ALL\n")
                    file70.write("!*\n")
                    #*******************************
                    #******PARSING RESULTS***********
                    for i in range(numel):
                        orientation=[]
                        # Initialize variables for the current element
                        curve_ast = []
                        curve_h = []
                        curve_geom_ast = []
                        curve_geom_h = []
                        
                        node_id=int(nomaster[i,0])
                        x=nomaster[i,1]
                        y=nomaster[i,2]
                        z=nomaster[i,3]
                        n_id=int(normal[i,0])
                        nx=normal[i,1]
                        ny=normal[i,2]
                        nz=normal[i,3]
                        if nx > 0:
                            orientation.append(1)
                        if nx < 0:
                            orientation.append(-1)
                        if ny > 0:
                            orientation.append(2)
                        if ny < 0:
                            orientation.append(-2)
                        if nz > 0:
                            orientation.append(3)
                        if nz < 0:
                            orientation.append(-3)
                                               
                        if results_bf:
                            for item in results_bf:
                                # filter quantity
                                for quantity in quantities:
                                    if item['quantity'] != quantity:
                                        continue
                                    for patch in item['patches']:                                    
                                        # Simple filter: skip if orientation doesn't match
                                        for o in orientation:
                                            if o is not None and patch.orientation != o:
                                                continue
                                            # Code to get data starts here
                                            # Build spatial coordinates for the patch (easier way to check node location)
                                            patch.buildSpace()
            
                                            lims = patch.lims  # [xmin, xmax, ymin, ymax, zmin, zmax]
            
                                            print(f"  Patch {i}: orientation {patch.orientation}, bounds {lims}")
            
                                            # Identify the planar axis (which pair in lims is equal since face is planar)
                                            planar_axis = None
                                            plane_coord = None
            
                                            if abs(lims[0] - lims[1]) < 1e-9:  # x is constant (x-normal face)
                                                planar_axis = 0  # x-axis
                                                plane_coord = x
                                                plane_value = lims[0]
            
                                            elif abs(lims[2] - lims[3]) < 1e-9:  # y is constant (y-normal face)
                                                planar_axis = 1  # y-axis
                                                plane_coord = y
                                                plane_value = lims[2]
            
                                            elif abs(lims[4] - lims[5]) < 1e-9:  # z is constant (z-normal face)
                                                planar_axis = 2  # z-axis
                                                plane_coord = z
                                                plane_value = lims[4]
            
                                            else:
                                                print(f"    Warning: Could not identify planar axis for patch {i}")
                                                continue
            
                                            # Tolerance for distance to plane
                                            tolerance = config['cell_size']
            
                                            # Check if distance between point and plane is less than tolerance
                                            plane_distance = abs(plane_coord - plane_value)
                                            if plane_distance > tolerance:
                                                print(f"    Point too far from plane: distance={plane_distance:.6f}, tolerance={tolerance:.6f}")
                                                continue
            
                                            # Check if point is within limits of other two directions
                                            in_bounds = True
            
                                            if planar_axis == 0:  # x-normal face, check y,z bounds
                                                if not (lims[2] <= y <= lims[3] and lims[4] <= z <= lims[5]):
                                                    in_bounds = False
                                            elif planar_axis == 1:  # y-normal face, check x,z bounds
                                                if not (lims[0] <= x <= lims[1] and lims[4] <= z <= lims[5]):
                                                    in_bounds = False
                                            elif planar_axis == 2:  # z-normal face, check x,y bounds
                                                if not (lims[0] <= x <= lims[1] and lims[2] <= y <= lims[3]):
                                                    in_bounds = False
            
                                            if not in_bounds:
                                                print(f"    Point outside patch bounds in non-planar directions")
                                                continue
            
                                            # Find closest CELL CENTER in the patch
                                            if planar_axis == 0:
                                                distances = np.sqrt((patch.y_centers - y) ** 2 + (patch.z_centers - z) ** 2)
                                            elif planar_axis == 1:
                                                distances = np.sqrt((patch.x_centers - x) ** 2 + (patch.z_centers - z) ** 2)
                                            elif planar_axis == 2:
                                                distances = np.sqrt((patch.x_centers - x) ** 2 + (patch.y_centers - y) ** 2)
                                            
                                            closest_idx = np.unravel_index(np.argmin(distances), distances.shape)
                                            closest_i, closest_j = closest_idx
                                            
                                            # Extract time series at this location (NO -1 needed now!)
                                            time_series = patch.data[closest_j, closest_i, :]

            
                                            print(f"    Grid indices: {closest_idx}")
                                            print(f"    Data range: {time_series.min():.2f} to {time_series.max():.2f}")
            
                                            if quantity=='ADIABATIC SURFACE TEMPERATURE': curve_ast.append({
                                                'times': times,
                                                'time_series': time_series,
                                                'orientation': o,
                                                # 'array_indices': closest_idx,
                                                # 'patch_bounds': lims,
                                                # 'units': units,
                                                # 'planar_axis': planar_axis,
                                                # 'plane_distance': plane_distance,
                                                # 'tolerance_used': tolerance
                                            })
                                            if quantity=='HEAT TRANSFER COEFFICIENT': curve_h.append({
                                                'times': times,
                                                'time_series': time_series,
                                                'orientation': o,
                                                # 'array_indices': closest_idx,
                                                # 'patch_bounds': lims,
                                                # 'units': units,
                                                # 'planar_axis': planar_axis,
                                                # 'plane_distance': plane_distance,
                                                # 'tolerance_used': tolerance
                                            })
                                
                        #***************************************
                        if results_be:
                            curve_geom_ast=[]
                            curve_geom_h=[]
                            for item in results_be:
                                # filter quantity
                                for quantity in quantities:
                                    if item['quantity'] != quantity:
                                        continue
                                    #if quantity=='ADIABATIC SURFACE TEMPERATURE': curve_geom_ast=[]
                                    #if quantity=='HEAT TRANSFER COEFFICIENT': curve_geom_h=[]
                                    values=item['values']
                                    query_point = np.array([x, y, z])
                                    # Calculate distances from query point to all face centroids
                                    distances = np.sqrt(np.sum((item['all_face_centroids'] - query_point)**2, axis=1))
                                    closest_face_idx = np.argmin(distances)
                                    closest_distance = distances[closest_face_idx]
                                    if closest_distance > config['cell_size']:
                                        continue
                                    closest_face_centroid = item['all_face_centroids'][closest_face_idx]
                                    # Get the time series for the closest face
                                    if closest_face_idx < values.shape[0]:
                                        time_series = values[closest_face_idx, :]
                                        print(f"  Closest face: index {closest_face_idx}, distance {closest_distance:.6f}")
                                        print(f"  Face centroid coords: ({closest_face_centroid[0]:.6f}, {closest_face_centroid[1]:.6f}, {closest_face_centroid[2]:.6f})")
                                        print(f"  Data range: {time_series.min():.2f} to {time_series.max():.2f}")
                                        # Get the vertices of the closest face for additional info
                                        closest_face_vertices = item['all_vertices'][item['all_faces'][closest_face_idx]-1]
                                        
                                        if quantity=='ADIABATIC SURFACE TEMPERATURE': 
                                            curve_geom_ast.append({
                                                # 'source_type': 'GEOM',
                                                # 'geom_file': gcf_file,
                                                # 'data_file': be_file,
                                                # 'mesh_id': file_quantities[filename]['mesh'],
                                                # 'face_index': closest_face_idx,
                                                'times': times,
                                                'time_series': time_series,
                                                # 'closest_face_centroid': closest_face_centroid,
                                                # 'closest_face_vertices': closest_face_vertices,
                                                # 'distance_to_query': closest_distance,
                                                # 'query_point': query_point,
                                                # 'quantity': file_quantity,
                                                # 'n_vertices': len(vertices),
                                                # 'n_faces': len(faces),
                                                # 'all_face_centroids': face_centroids,  # Include all centroids for reference
                                                # 'all_vertices': vertices,              # Include all vertices for reference
                                                # 'all_faces': faces                     # Include all faces for reference
                                            })
                                        if quantity=='HEAT TRANSFER COEFFICIENT': 
                                            curve_geom_h.append({
                                                # 'source_type': 'GEOM',
                                                # 'geom_file': gcf_file,
                                                # 'data_file': be_file,
                                                # 'mesh_id': file_quantities[filename]['mesh'],
                                                # 'face_index': closest_face_idx,
                                                'times': times,
                                                'time_series': time_series,
                                                # 'closest_face_centroid': closest_face_centroid,
                                                # 'closest_face_vertices': closest_face_vertices,
                                                # 'distance_to_query': closest_distance,
                                                # 'query_point': query_point,
                                                # 'quantity': file_quantity,
                                                # 'n_vertices': len(vertices),
                                                # 'n_faces': len(faces),
                                                # 'all_face_centroids': face_centroids,  # Include all centroids for reference
                                                # 'all_vertices': vertices,              # Include all vertices for reference
                                                # 'all_faces': faces                     # Include all faces for reference
                                            })
                        # Process AST Data
                        if curve_geom_ast and curve_ast:
                            if len(curve_ast) > 1: 
                                time_avg, tast_avg = average_results_by_orientation(curve_ast)
                            else:
                                time_avg = curve_ast[0]['times']
                                tast_avg = curve_ast[0]['time_series']
                            
                            time_avg = np.concatenate([[0.0], time_avg])
                            tast_avg = np.concatenate([[20.0], tast_avg])
                            
                            time_geom_avg, tast_geom_avg = average_multiple_results(curve_geom_ast)
                            curve = [{'times': time_avg, 'time_series': (tast_avg + tast_geom_avg) / 2}]
                            
                        elif curve_geom_ast:
                            time_geom_avg, tast_geom_avg = average_multiple_results(curve_geom_ast)
                            curve = [{'times': time_geom_avg, 'time_series': tast_geom_avg}]
                            
                        elif curve_ast:
                            if len(curve_ast) > 1: 
                                time_avg, tast_avg = average_results_by_orientation(curve_ast)
                            else:
                                time_avg = curve_ast[0]['times']
                                tast_avg = curve_ast[0]['time_series']
                                
                            time_avg = np.concatenate([[0.0], time_avg])
                            tast_avg = np.concatenate([[20.0], tast_avg])
                            curve = [{'times': time_avg, 'time_series': tast_avg, 'orientation': 0}]
                            
                        else:
                            print(f'Warning: No AST data found for element {i+1}')
                            # Provide a safe ambient fallback to prevent crashing downstream
                            curve = [{'times': np.array([0.0, 18000.0]), 'time_series': np.array([20.0, 20.0])}]

                        # Process Heat Transfer Coefficient (h) Data
                        if curve_geom_h and curve_h:
                            if len(curve_h) > 1: 
                                time_avg, h_avg = average_results_by_orientation(curve_h)
                            else:
                                time_avg = curve_h[0]['times']
                                h_avg = curve_h[0]['time_series']
                                
                            time_avg = np.concatenate([[0.0], time_avg])
                            h_avg = np.concatenate([[0.0], h_avg])
                            
                            time_geom_avg, h_geom_avg = average_multiple_results(curve_geom_h)
                            curve2 = [{'times': time_avg, 'time_series': (h_avg + h_geom_avg) / 2}]
                            
                        elif curve_geom_h:
                            time_geom_avg, h_geom_avg = average_multiple_results(curve_geom_h)
                            curve2 = [{'times': time_geom_avg, 'time_series': h_geom_avg}]
                            
                        elif curve_h:
                            if len(curve_h) > 1: 
                                time_avg, h_avg = average_results_by_orientation(curve_h)
                            else:
                                time_avg = curve_h[0]['times']
                                h_avg = curve_h[0]['time_series']
                                
                            time_avg = np.concatenate([[0.0], time_avg])
                            h_avg = np.concatenate([[0.0], h_avg])
                            curve2 = [{'times': time_avg, 'time_series': h_avg, 'orientation': 0}]
                            
                        else:
                            print(f'Warning: No h data found for element {i+1}')
                            # Provide a safe baseline fallback to prevent crashing downstream
                            curve2 = [{'times': np.array([0.0, 18000.0]), 'time_series': np.array([0.0, 0.0])}]
                        
                        #*********** CREATING TABLES (.DAT) **********
                        # Convert node ID to string with formatting (equivalent to g8.0 format)
                        intfile2 = f"{int(nomaster[i, 0]):g}"
                        print(intfile2)
                        
                        if config['average']['mode']=='none':
                            linhas = len(curve[0]['times'])
                        elif config['average']['mode']=='mean_last_n': 
                            linhas = len(curve[0]['times']) + 2
                        # elif config['average']['mode']=='window':
                        #     linhas = len(curve[0]['times']) + 1
                        elif config['average']['mode']=='time_window':
                            linhas =  2
                        
                        for quantity in quantities:
                            # Write temperature table definition
                            if quantity=='ADIABATIC SURFACE TEMPERATURE':
                                file70.write(f"*DIM,A{intfile2},TABLE,{linhas:8d},1,1,TIME,TEMP,\n")
                                file70.write("!*\n")
                                if config['average']['mode']=='none':
                                    for j in range(len(curve[0]['time_series'])):
                                        file70.write(f"*set,A{intfile2}({j+1:8d},0),{curve[0]['times'][j]:12.5e}\n")
                                        file70.write(f"*set,A{intfile2}({j+1:8d},1),{curve[0]['time_series'][j]:12.5e}\n")
                                elif config['average']['mode']=='time_window':
                                    time_array = curve[0]['times']
                                    series_array = curve[0]['time_series']
                                    start_time=config['average']['t0']
                                    end_time=config['average']['t1']
                                    # Calculate averaging parameters
                                    average_window = apply_time_window_average(curve[0]['times'], curve[0]['time_series'], start_time, end_time)
                                    file70.write(f"*set,A{intfile2}(1,0),{0.0:12.5e}\n")
                                    file70.write(f"*set,A{intfile2}(1,1),{average_window:12.5e}\n")
                                    file70.write(f"*set,A{intfile2}(2,0),{18000.0:12.5e}\n")
                                    file70.write(f"*set,A{intfile2}(2,1),{average_window:12.5e}\n")
                                elif config['average']['mode']=='mean_last_n':
                                    for j in range(len(curve[0]['time_series'])):
                                        file70.write(f"*set,A{intfile2}({j+1},0),{curve[0]['times'][j]:12.5e}\n")
                                        file70.write(f"*set,A{intfile2}({j+1},1),{curve[0]['time_series'][j]:12.5e}\n")
                                    avg_last_n = np.mean(curve[0]['time_series'][-int(config['average']['n']):])
                                    file70.write(f"*set,A{intfile2}({j+2},0),{curve[0]['times'][-1]+round(curve[0]['times'][1]-curve[0]['times'][0],1):12.5e}\n")
                                    file70.write(f"*set,A{intfile2}({j+2},1),{avg_last_n:12.5e}\n")
                                    file70.write(f"*set,A{intfile2}({j+3},0),{18000.0:12.5e}\n")
                                    file70.write(f"*set,A{intfile2}({j+3},1),{avg_last_n:12.5e}\n")
                            # Write heat transfer coefficient table if VARIABLE==2
                            if quantity=='HEAT TRANSFER COEFFICIENT':
                                file70.write(f"*DIM,H{intfile2},TABLE,{linhas:8d},1,1,TIME,,\n")
                                file70.write("!*\n")
                                if config['average']['mode']=='none':
                                    for j in range(len(curve2[0]['time_series'])):
                                        file70.write(f"*set,H{intfile2}({j+1:8d},0),{curve2[0]['times'][j]:12.5e}\n")
                                        file70.write(f"*set,H{intfile2}({j+1:8d},1),{curve2[0]['time_series'][j]:12.5e}\n")
                                elif config['average']['mode']=='time_window':
                                    time_array = curve2[0]['times']
                                    series_array = curve2[0]['time_series']
                                    start_time=config['average']['t0']
                                    end_time=config['average']['t1']
                                    # Calculate averaging parameters
                                    average_window = apply_time_window_average(curve2[0]['times'], curve2[0]['time_series'], start_time, end_time)
                                    file70.write(f"*set,H{intfile2}(1,0),{0.0:12.5e}\n")
                                    file70.write(f"*set,H{intfile2}(1,1),{average_window:12.5e}\n")
                                    file70.write(f"*set,H{intfile2}(2,0),{18000.0:12.5e}\n")
                                    file70.write(f"*set,H{intfile2}(2,1),{average_window:12.5e}\n")
                                elif config['average']['mode']=='mean_last_n':
                                    for j in range(len(curve2[0]['time_series'])):
                                        file70.write(f"*set,H{intfile2}({j+1},0),{curve2[0]['times'][j]:12.5e}\n")
                                        file70.write(f"*set,H{intfile2}({j+1},1),{curve2[0]['time_series'][j]:12.5e}\n")
                                    avg_last_n = np.mean(curve2[0]['time_series'][-int(config['average']['n']):])
                                    file70.write(f"*set,H{intfile2}({j+2},0),{curve2[0]['times'][-1]+round(curve2[0]['times'][1]-curve2[0]['times'][0],1):12.5e}\n")
                                    file70.write(f"*set,H{intfile2}({j+2},1),{avg_last_n:12.5e}\n")
                                    file70.write(f"*set,H{intfile2}({j+3},0),{18000.0:12.5e}\n")
                                    file70.write(f"*set,H{intfile2}({j+3},1),{avg_last_n:12.5e}\n")
                    

                # *********** APPLYING TABLES AS NODAL LOADS **********

                file71.write("/PREP7\n")
                print('LOOP_CARGAS')
                for i in range(nodes_elements['numel']):
                    intfile2 = int(nodes_elements['nomaster'][i, 0])
                    print(intfile2)
                    file71.write("!*\n")
                    file71.write(f"D,{intfile2:8d}, , %A{intfile2:8d}% , , , ,TEMP, , , , ,\n")

                # *******************************
                # *********** APPLYING TABLES AS CONVECTIVE HEAT TRANSFER COEFFICIENT **********
                if config['variables']['mode']=='AST_and_h':
                    file71.write("/PREP7\n")
                    print('LOOP_HEAT_TRANSFER')    
                    for i in range(nodes_elements['numel']):
                        intfile2 = int(nodes_elements['nomaster'][i, 0])
                        print(intfile2)
                        intfile3 = int(nodes_elements['elements'][i, 0])
                        print(intfile3)
                        file71.write(f"SFE,{intfile3:8d},1,CONV,0,%H{intfile2:8d}%\n")

                if config['variables']['mode']=='AST_only':
                    file71.write("/PREP7\n")
                    print('LOOP_HEAT_TRANSFER')
                    hc=config['variables']['h']
                    for i in range(nodes_elements['numel']):
                        intfile2 = int(nodes_elements['nomaster'][i, 0])
                        print(intfile2)
                        intfile3 = int(nodes_elements['elements'][i, 0])
                        print(intfile3)
                        file71.write(f"SFE,{intfile3:8d},1,CONV,0,{hc:7.3f}\n")
                
                
                file70.write(f'/input,{config["file"].stem}_loads,dat')
####################################################
### NEW CASE - FTMI #############
        case 'ftmi':
            # Read FTMI-specific node file
            nodes = read_nodes_ftmi()
            
            print(f"Number of nodes: {nodes['nnode']}")
            print(f"Node IDs: {nodes['node_ids']}")
            
            # Initialize data storage for all nodes
            node_data = {}
            
            # Process each node
            for idx in range(nodes['nnode']):
                node_id = nodes['node_ids'][idx]
                x, y, z = nodes['coordinates'][idx]
                if not np.isnan(nodes['orientations'][idx]):
                    orientation = [int(nodes['orientations'][idx])]
                    bndf_file_pass=False
                else:
                    orientation = [] 
                    bndf_file_pass=True
                
                print(f"\n{'='*60}")
                print(f"Processing node: {node_id}")
                print(f"  Position: ({x:.3f}, {y:.3f}, {z:.3f})")
                if bndf_file_pass==False: print(f"  Orientation: {orientation[0]}")
                
                # Initialize result storage for this node
                curve_ast = []
                curve_h = []
                curve_geom_ast = []
                curve_geom_h = []
                
                # Search in BNDF files (planar boundary files)
                if results_bf and bndf_file_pass==False:
                    for item in results_bf:
                        # Filter by quantity
                        for quantity in quantities:
                            if item['quantity'] != quantity:
                                continue
                            
                            for patch in item['patches']:
                                # Check if orientation matches
                                if patch.orientation not in orientation:
                                    continue
                                
                                # Build spatial coordinates for the patch
                                patch.buildSpace()
                                lims = patch.lims  # [xmin, xmax, ymin, ymax, zmin, zmax]
                                
                                # Identify the planar axis
                                planar_axis = None
                                plane_coord = None
                                
                                if abs(lims[0] - lims[1]) < 1e-9:  # x is constant
                                    planar_axis = 0
                                    plane_coord = x
                                    plane_value = lims[0]
                                elif abs(lims[2] - lims[3]) < 1e-9:  # y is constant
                                    planar_axis = 1
                                    plane_coord = y
                                    plane_value = lims[2]
                                elif abs(lims[4] - lims[5]) < 1e-9:  # z is constant
                                    planar_axis = 2
                                    plane_coord = z
                                    plane_value = lims[4]
                                else:
                                    continue
                                
                                # Check distance to plane
                                tolerance = config['cell_size']
                                plane_distance = abs(plane_coord - plane_value)
                                if plane_distance > tolerance:
                                    continue
                                
                                # Check if point is within bounds
                                in_bounds = True
                                if planar_axis == 0:  # x-normal
                                    if not (lims[2] <= y <= lims[3] and lims[4] <= z <= lims[5]):
                                        in_bounds = False
                                elif planar_axis == 1:  # y-normal
                                    if not (lims[0] <= x <= lims[1] and lims[4] <= z <= lims[5]):
                                        in_bounds = False
                                elif planar_axis == 2:  # z-normal
                                    if not (lims[0] <= x <= lims[1] and lims[2] <= y <= lims[3]):
                                        in_bounds = False
                                
                                if not in_bounds:
                                    continue
                                
                                # Find closest CELL CENTER in the patch
                                if planar_axis == 0:
                                    distances = np.sqrt((patch.y_centers - y) ** 2 + (patch.z_centers - z) ** 2)
                                elif planar_axis == 1:
                                    distances = np.sqrt((patch.x_centers - x) ** 2 + (patch.z_centers - z) ** 2)
                                elif planar_axis == 2:
                                    distances = np.sqrt((patch.x_centers - x) ** 2 + (patch.y_centers - y) ** 2)
                                
                                closest_idx = np.unravel_index(np.argmin(distances), distances.shape)
                                closest_i, closest_j = closest_idx
                                
                                # Extract time series at this location (NO -1 needed now!)
                                time_series = patch.data[closest_j, closest_i, :]

                                
                                print(f"  Found BNDF match: {quantity}")
                                print(f"    Range: {time_series.min():.2f} to {time_series.max():.2f}")
                                
                                # Store results
                                if quantity == 'ADIABATIC SURFACE TEMPERATURE':
                                    curve_ast.append({
                                        'times': times,
                                        'time_series': time_series,
                                        'orientation': orientation[0]
                                    })
                                
                                if quantity == 'HEAT TRANSFER COEFFICIENT':
                                    curve_h.append({
                                        'times': times,
                                        'time_series': time_series,
                                        'orientation': orientation[0]
                                    })
                
                # Search in BE files (complex geometry)
                if results_be:
                    for item in results_be:
                        # Filter by quantity
                        for quantity in quantities:
                            if item['quantity'] != quantity:
                                continue
                            
                            query_point = np.array([x, y, z])
                            
                            # Find closest face
                            distances = np.sqrt(np.sum(
                                (item['all_face_centroids'] - query_point)**2, axis=1))
                            closest_face_idx = np.argmin(distances)
                            closest_distance = distances[closest_face_idx]
                            
                            if closest_distance > config['cell_size']:
                                continue
                            
                            # Get time series
                            values = item['values']
                            if closest_face_idx < values.shape[0]:
                                time_series = values[closest_face_idx, :]
                                
                                print(f"  Found BE match: {quantity}")
                                print(f"    Distance: {closest_distance:.6f}")
                                print(f"    Range: {time_series.min():.2f} to {time_series.max():.2f}")
                                
                                if quantity == 'ADIABATIC SURFACE TEMPERATURE':
                                    curve_geom_ast.append({
                                        'times': times,
                                        'time_series': time_series
                                    })
                                
                                if quantity == 'HEAT TRANSFER COEFFICIENT':
                                    curve_geom_h.append({
                                        'times': times,
                                        'time_series': time_series
                                    })
                
                # Process and average results for AST
                final_ast = None
                final_h = None
                final_times = None
                
                try:
                    # Combine GEOM and regular results for AST
                    if curve_geom_ast and curve_ast:
                        if len(curve_ast) > 1:
                            time_avg, tast_avg = average_results_by_orientation(curve_ast)
                        else:
                            time_avg = curve_ast[0]['times']
                            tast_avg = curve_ast[0]['time_series']
                        
                        time_geom_avg, tast_geom_avg = average_multiple_results(curve_geom_ast)
                        
                        final_times = time_avg
                        final_ast = (tast_avg + tast_geom_avg) / 2
                    
                    elif curve_geom_ast:
                        time_geom_avg, tast_geom_avg = average_multiple_results(curve_geom_ast)
                        final_times = time_geom_avg
                        final_ast = tast_geom_avg
                    
                    elif curve_ast:
                        if len(curve_ast) > 1:
                            time_avg, tast_avg = average_results_by_orientation(curve_ast)
                        else:
                            time_avg = curve_ast[0]['times']
                            tast_avg = curve_ast[0]['time_series']
                        
                        final_times = time_avg
                        final_ast = tast_avg
                    
                    # Add initial conditions
                    if final_ast is not None:
                        final_times = np.concatenate([[0.0], final_times])
                        final_ast = np.concatenate([[20.0], final_ast])
                    
                except Exception as e:
                    print(f"  Error processing AST: {e}")
                
                # Process results for h
                try:
                    if curve_geom_h and curve_h:
                        if len(curve_h) > 1:
                            time_avg, h_avg = average_results_by_orientation(curve_h)
                        else:
                            time_avg = curve_h[0]['times']
                            h_avg = curve_h[0]['time_series']
                        
                        time_geom_avg, h_geom_avg = average_multiple_results(curve_geom_h)
                        final_h = (h_avg + h_geom_avg) / 2
                    
                    elif curve_geom_h:
                        time_geom_avg, h_geom_avg = average_multiple_results(curve_geom_h)
                        final_h = h_geom_avg
                    
                    elif curve_h:
                        if len(curve_h) > 1:
                            time_avg, h_avg = average_results_by_orientation(curve_h)
                        else:
                            time_avg = curve_h[0]['times']
                            h_avg = curve_h[0]['time_series']
                        
                        final_h = h_avg
                    
                    # Add initial conditions
                    if final_h is not None:
                        final_h = np.concatenate([[0.0], final_h])
                
                except Exception as e:
                    print(f"  Error processing h: {e}")
                
                # Apply time window averaging if requested
                if config['average']['mode'] == 'time_window' and final_times is not None:
                    start_time = config['average']['t0']
                    end_time = config['average']['t1']
                    
                    if final_ast is not None:
                        ast_avg = apply_time_window_average(
                            final_times, final_ast, start_time, end_time)
                        final_times = np.array([0.0, 18000.0])
                        final_ast = np.array([ast_avg, ast_avg])
                    
                    if final_h is not None:
                        h_avg = apply_time_window_average(
                            final_times, final_h, start_time, end_time)
                        final_h = np.array([h_avg, h_avg])
                
                # Apply mean last N averaging if requested
                elif config['average']['mode'] == 'mean_last_n' and final_times is not None:
                    n = int(config['average']['n'])
                    
                    if final_ast is not None and len(final_ast) >= n:
                        ast_avg = np.mean(final_ast[-n:])
                        dt = final_times[1] - final_times[0] if len(final_times) > 1 else 1.0
                        final_times = np.concatenate([
                            final_times, [final_times[-1] + dt, 18000.0]])
                        final_ast = np.concatenate([final_ast, [ast_avg, ast_avg]])
                    
                    if final_h is not None and len(final_h) >= n:
                        h_avg = np.mean(final_h[-n:])
                        final_h = np.concatenate([final_h, [h_avg, h_avg]])
                
                # Store node data
                if final_times is not None:
                    node_data[node_id] = {
                        'time': final_times,
                        'AST': final_ast if config['variables']['mode'] in ['AST_only', 'AST_and_h'] else None,
                        'h': final_h if config['variables']['mode'] == 'AST_and_h' else None
                    }
                    print(f"  ✓ Data extracted for {node_id}")
                else:
                    print(f"  ✗ No data found for {node_id}")
            
            # Write results to CSV
            print(f"\n{'='*60}")
            print("Writing results to CSV...")
            write_ftmi_csv(node_data, config)
            print(f"{'='*60}")
