import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

cm = pd.read_csv(r"DASCH_Kp_Replication\dasch_northern_reproduction_v2\dasch_northern_reproduction\cluster_requery\cluster_members.csv")

# Check pos 271 across all dates — do sources reappear near same coordinates?
pos271 = cm[cm['pos_idx'] == 271].copy()
print("All pos 271 sources across dates:")
print(pos271[['obs_date','src_index','src_ra_deg','src_dec_deg','plate_mag']].to_string())

print("\n\nAll pos 319 sources across dates:")
pos319 = cm[cm['pos_idx'] == 319].copy()
print(pos319[['obs_date','src_index','src_ra_deg','src_dec_deg','plate_mag']].to_string())

# Cross-date proximity check for pos 271
print("\n\nCross-date separations for pos 271:")
coords271 = SkyCoord(ra=pos271['src_ra_deg'].values*u.deg,
                     dec=pos271['src_dec_deg'].values*u.deg)
dates271 = pos271['obs_date'].values
GEO_KM = 42164

for i in range(len(coords271)):
    for j in range(i+1, len(coords271)):
        if dates271[i] == dates271[j]:
            continue  # skip same-date, already done
        sep_arcsec = coords271[i].separation(coords271[j]).arcsec
        sep_km = GEO_KM * np.tan(np.radians(sep_arcsec/3600))
        if sep_arcsec < 600:  # only show close cross-date pairs
            print(f"  {dates271[i]} src{pos271['src_index'].iloc[i]} <-> {dates271[j]} src{pos271['src_index'].iloc[j]}: {sep_arcsec:.1f} arcsec = {sep_km:.1f} km")