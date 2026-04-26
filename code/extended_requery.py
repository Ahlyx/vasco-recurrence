#!/usr/bin/env python3
"""
extended_requery.py
====================
Extended cluster requery targeting high-transient-count grid positions
identified in the DASCH fast scan. Adapted from Cann (2026) 
dasch_cluster_requery.py to target new positions beyond the original
271 and 319.

Reads transients_all.csv to find same-date multi-source events at
target positions, then queries DASCH for per-source coordinates.

Output: data/extended_cluster_members.csv
"""

import os
import sys
import numpy as np
import pandas as pd
import astropy.units as u
import astropy.coordinates as coord
from pathlib import Path
import daschlab

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
TRANSIENTS_FILE = os.path.join(DATA_DIR, "transients_all.csv")
OUTDIR = Path(DATA_DIR)

# Target positions — high same-date cluster counts from NewData.py
# Excluding pos 144, 168, 169, 143 (southern sky, Dec < -20, sparse DASCH coverage)
# Prioritizing near-equatorial positions with good DASCH coverage
TARGET_POSITIONS = [270, 11, 269]  # pos 144 etc added later if these work

MIN_CLUSTER_SIZE = 2  # minimum sources on same date to qualify

def build_cluster_list(df, target_positions, min_size):
    """Find same-date multi-source events at target positions."""
    clusters = []
    for pos_id in target_positions:
        subset = df[df['pos_idx'] == pos_id]
        grid_ra = subset['ra'].iloc[0]
        grid_dec = subset['dec'].iloc[0]
        
        by_date = subset.groupby('obs_date')
        for date, group in by_date:
            if len(group) >= min_size:
                members = list(zip(
                    group['src_index'].astype(int).tolist(),
                    group['mag'].tolist()
                ))
                clusters.append((pos_id, grid_ra, grid_dec, date, members))
                print(f"  Pos {pos_id} | {date} | {len(members)} sources")
    return clusters

def lookup_position(pos_idx, grid_ra, grid_dec, outdir):
    print(f"\n--- pos_idx {pos_idx}: RA={grid_ra}, Dec={grid_dec:+.5f} ---")
    root = outdir / f"session_{pos_idx}"
    root.mkdir(exist_ok=True)
    sess = daschlab.open_session(root=str(root), interactive=False)
    sess.select_target(ra_deg=float(grid_ra), dec_deg=float(grid_dec))
    sess.select_refcat("apass")
    refcat = sess.refcat()
    exps = sess.exposures()
    print(f"  Refcat: {len(refcat)} sources | Exposures: {len(exps)}")
    return sess, refcat, exps

def get_ra_dec(refcat):
    if "pos" in refcat.colnames:
        return refcat["pos"].ra.deg, refcat["pos"].dec.deg
    for ra_col in ("ra_deg", "ra", "RA"):
        if ra_col in refcat.colnames:
            for dec_col in ("dec_deg", "dec", "Dec"):
                if dec_col in refcat.colnames:
                    return np.array(refcat[ra_col]), np.array(refcat[dec_col])
    raise ValueError("Cannot find RA/Dec columns in refcat")

def find_plate_for_date(exps, obs_date):
    dates = np.array([str(d)[:10] for d in exps["obs_date"]])
    return np.where(dates == obs_date)[0]

def plate_name(exp_row):
    try:
        return f"{exp_row['series']}{int(exp_row['platenum']):05d}"
    except Exception:
        return "?"

def run_requery(clusters, outdir):
    all_members = []
    
    # Group by position to avoid reopening sessions
    by_pos = {}
    for c in clusters:
        by_pos.setdefault(c[0], []).append(c)
    
    for pos_idx, cluster_list in by_pos.items():
        grid_ra = cluster_list[0][1]
        grid_dec = cluster_list[0][2]
        
        try:
            sess, refcat, exps = lookup_position(pos_idx, grid_ra, grid_dec, outdir)
            ra_all, dec_all = get_ra_dec(refcat)
        except Exception as e:
            print(f"  ERROR opening session for pos {pos_idx}: {e}")
            continue
        
        for cluster in cluster_list:
            _, _, _, obs_date, members = cluster
            print(f"\n  Cluster pos={pos_idx} date={obs_date} ({len(members)} members)")
            
            plate_idx = find_plate_for_date(exps, obs_date)
            plate_names = [plate_name(exps[int(i)]) for i in plate_idx] if len(plate_idx) > 0 else []
            print(f"  Plates: {';'.join(plate_names) if plate_names else 'none found'}")
            
            member_coords = []
            for src_index, mag in members:
                if src_index >= len(refcat):
                    print(f"    src {src_index} out of bounds (refcat size {len(refcat)})")
                    continue
                ra_s = float(ra_all[int(src_index)])
                dec_s = float(dec_all[int(src_index)])
                print(f"    src {src_index}: RA={ra_s:.5f} Dec={dec_s:+.5f} mag={mag:.2f}")
                member_coords.append((src_index, ra_s, dec_s, mag))
                all_members.append({
                    "pos_idx": pos_idx,
                    "grid_ra": grid_ra,
                    "grid_dec": grid_dec,
                    "obs_date": obs_date,
                    "src_index": src_index,
                    "src_ra_deg": ra_s,
                    "src_dec_deg": dec_s,
                    "plate_mag": mag,
                    "plate_names": ";".join(plate_names),
                })
            
            # Pairwise separations
            if len(member_coords) >= 2:
                seps = []
                for a in range(len(member_coords)):
                    for b in range(a+1, len(member_coords)):
                        sa = coord.SkyCoord(ra=member_coords[a][1]*u.deg, dec=member_coords[a][2]*u.deg)
                        sb = coord.SkyCoord(ra=member_coords[b][1]*u.deg, dec=member_coords[b][2]*u.deg)
                        seps.append(sa.separation(sb).arcsec)
                print(f"    Separations: min={min(seps):.1f}\" max={max(seps):.1f}\" median={np.median(seps):.1f}\"")
    
    return all_members

if __name__ == "__main__":
    print("Loading transients...")
    df = pd.read_csv(TRANSIENTS_FILE)
    
    session_dir = Path(DATA_DIR) / "requery_sessions"
    session_dir.mkdir(exist_ok=True)
    
    print(f"\nBuilding cluster list for positions: {TARGET_POSITIONS}")
    clusters = build_cluster_list(df, TARGET_POSITIONS, MIN_CLUSTER_SIZE)
    print(f"\nTotal clusters to query: {len(clusters)}")
    
    if len(clusters) == 0:
        print("No clusters found. Check TARGET_POSITIONS and MIN_CLUSTER_SIZE.")
        sys.exit(1)
    
    print("\nStarting DASCH requery...")
    members = run_requery(clusters, session_dir)
    
    if members:
        out_path = os.path.join(DATA_DIR, "extended_cluster_members.csv")
        pd.DataFrame(members).to_csv(out_path, index=False)
        print(f"\nWrote {len(members)} source records to {out_path}")
    else:
        print("\nNo members retrieved.")