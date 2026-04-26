import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from itertools import combinations

cm = pd.read_csv(r"DASCH_Kp_Replication\dasch_northern_reproduction_v2\dasch_northern_reproduction\cluster_requery\cluster_members.csv")

GEO_KM = 42164

# For pos 271, build a recurrence map
pos271 = cm[cm['pos_idx'] == 271].copy().reset_index(drop=True)
coords = SkyCoord(ra=pos271['src_ra_deg'].values*u.deg,
                  dec=pos271['src_dec_deg'].values*u.deg)

# Bin sources into 600 arcsec spatial cells and count cross-date recurrences
print("RECURRENCE ANALYSIS - Pos 271")
print("Sources within 200 arcsec across different dates:\n")

recurrence_groups = []
used = set()

for i, j in combinations(range(len(pos271)), 2):
    if pos271['obs_date'].iloc[i] == pos271['obs_date'].iloc[j]:
        continue
    sep = coords[i].separation(coords[j]).arcsec
    if sep < 200:
        print(f"  CLOSE PAIR: {pos271['obs_date'].iloc[i]} src{pos271['src_index'].iloc[i]} "
              f"<-> {pos271['obs_date'].iloc[j]} src{pos271['src_index'].iloc[j]}")
        print(f"    Sep: {sep:.1f} arcsec = {GEO_KM * np.tan(np.radians(sep/3600)):.1f} km at GEO")
        print(f"    Mags: {pos271['plate_mag'].iloc[i]:.2f} vs {pos271['plate_mag'].iloc[j]:.2f}")
        dt = abs(pd.Timestamp(pos271['obs_date'].iloc[i]) - 
                 pd.Timestamp(pos271['obs_date'].iloc[j])).days
        print(f"    Time separation: {dt} days ({dt/365.25:.2f} years)\n")