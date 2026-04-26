"""
Same-date scale check: separations between sources detected on the same date at the same position.

Reads:  data/cluster_members.csv
Outputs: for each (pos, date) group with >=2 sources, pairwise separations in arcsec and km at GEO.
"""
import os
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
CLUSTER_FILE = os.path.join(DATA_DIR, "cluster_members.csv")

cm = pd.read_csv(CLUSTER_FILE)

GEO_KM = 42164

for (pos, date), group in cm.groupby(['pos_idx', 'obs_date']):
    if len(group) < 2:
        continue
    coords = SkyCoord(ra=group['src_ra_deg'].values * u.deg,
                      dec=group['src_dec_deg'].values * u.deg)
    print(f"\nPos {pos} | {date} | {len(group)} sources | plates: {group['plate_names'].iloc[0][:30]}")
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            sep_arcsec = coords[i].separation(coords[j]).arcsec
            sep_km = GEO_KM * np.tan(np.radians(sep_arcsec / 3600))
            print(f"  src {group['src_index'].iloc[i]} <-> src {group['src_index'].iloc[j]}: "
                  f"{sep_arcsec:.1f} arcsec = {sep_km:.1f} km at GEO")
