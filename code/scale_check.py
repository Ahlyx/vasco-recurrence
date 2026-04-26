"""
Cross-date proximity check: source tables and close cross-date separations for pos 271 and 319.

Reads:  data/cluster_members.csv
Outputs: full source table for pos 271 and 319; cross-date separations <600 arcsec for pos 271
         in arcsec and km at geosynchronous altitude.
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

pos271 = cm[cm['pos_idx'] == 271].copy()
print("All pos 271 sources across dates:")
print(pos271[['obs_date', 'src_index', 'src_ra_deg', 'src_dec_deg', 'plate_mag']].to_string())

print("\n\nAll pos 319 sources across dates:")
pos319 = cm[cm['pos_idx'] == 319].copy()
print(pos319[['obs_date', 'src_index', 'src_ra_deg', 'src_dec_deg', 'plate_mag']].to_string())

print("\n\nCross-date separations for pos 271:")
coords271 = SkyCoord(ra=pos271['src_ra_deg'].values * u.deg,
                     dec=pos271['src_dec_deg'].values * u.deg)
dates271 = pos271['obs_date'].values

for i in range(len(coords271)):
    for j in range(i + 1, len(coords271)):
        if dates271[i] == dates271[j]:
            continue
        sep_arcsec = coords271[i].separation(coords271[j]).arcsec
        sep_km = GEO_KM * np.tan(np.radians(sep_arcsec / 3600))
        if sep_arcsec < 600:
            print(f"  {dates271[i]} src{pos271['src_index'].iloc[i]} "
                  f"<-> {dates271[j]} src{pos271['src_index'].iloc[j]}: "
                  f"{sep_arcsec:.1f} arcsec = {sep_km:.1f} km")
