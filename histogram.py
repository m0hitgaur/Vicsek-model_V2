#!/usr/bin/env python3
"""
Histogram of ALL inter-particle distances r_ij
------------------------------------------------
Reads the CSV snapshots written by the C++ code and
produces one histogram that contains every r_ij (i≠j)
from every recorded time step.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import os
import glob
import argparse

# ---------- utilities -------------------------------------------------------
def minimum_image(dx, dy, Lx, Ly):
    """Apply minimum-image convention to separation vector."""
    if(dx<-Lx/2):dx=dx + Lx
    if(dx>Lx/2):dx=dx - Lx 
    if(dx<-Ly/2):dy=dy + Ly
    if(dx>Ly/2):dy=dy - Ly

    return dx, dy

def all_distances_one_frame(df, Lx, Ly):
    """
    Return every r_ij (i≠j) for a single snapshot dataframe.
    Shape of returned array: (N*(N-1)/2,)
    """
    pos = df[['x', 'y']].to_numpy()          # (N,2)
    N   = len(pos)
    r_list = []
    for i, j in combinations(range(N), 2):
        dx = pos[j, 0] - pos[i, 0]
        dy = pos[j, 1] - pos[i, 1]
        dx, dy = minimum_image(dx, dy, Lx, Ly)
        r = np.hypot(dx, dy)
        r_list.append(r)
    return np.array(r_list, dtype=np.float32)

# ---------- main -------------------------------------------------------------
def build_histogram(data_dir, trial, Lx, Ly, bins=100, r_max=None):
    """
    Accumulate r_ij over every recorded snapshot for one trial.
    Returns: counts, bin_edges (compatible with plt.stairs or plt.hist)
    """
    # --- locate every snapshot file -----------------------------------------
    pattern = os.path.join(data_dir, "config_data", f"trial_{trial}", "config_*.csv")
    files     = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No snapshots found for trial {trial} -> {pattern}")

    # --- loop over snapshots --------------------------------------------------
    all_r = []
    for fname in files:
        df = pd.read_csv(fname)
        all_r.append(all_distances_one_frame(df, Lx, Ly))
    all_r = np.concatenate(all_r)        # flatten into one long 1-D array

    # --- optional upper bound -----------------------------------------------
    if r_max is None:
        r_max = all_r.max()

    counts, edges = np.histogram(all_r, bins=bins, range=(0.0, r_max))
    return counts, edges, all_r          # return raw distances too if wanted

# ---------- CLI / standalone usage ------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inter-particle distance histogram")
    parser.add_argument("--data_dir", default="data", help="top-level data folder")
    parser.add_argument("--trial", type=int, default=0, help="trial index")
    parser.add_argument("--Lx", type=float, default=15.0, help="box length x")
    parser.add_argument("--Ly", type=float, default=15.0, help="box length y")
    parser.add_argument("--bins", type=int, default=150, help="number of histogram bins")
    parser.add_argument("--r_max", type=float, help="upper histogram limit (auto if omitted)")
    parser.add_argument("--png", default="rij_histogram.png", help="output figure file")
    args = parser.parse_args()

    counts, edges, raw = build_histogram(args.data_dir, args.trial,
                                         args.Lx, args.Ly,
                                         bins=args.bins, r_max=args.r_max)

    # ---------- plot ---------------------------------------------------------
    plt.figure(figsize=(6, 4))
    plt.stairs(counts, edges, fill=True, alpha=0.7, color="steelblue", edgecolor="black")
    plt.xlabel("r_ij")
    plt.ylabel("N(r_ij)")
    plt.title(f"Distribution of inter-particle distances  (trial {args.trial})")
    plt.tight_layout()
    plt.savefig(args.png, dpi=200)
    plt.show()