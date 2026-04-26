import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from itertools import combinations

cm = pd.read_csv(r"DASCH_Kp_Replication\dasch_northern_reproduction_v2\dasch_northern_reproduction\cluster_requery\cluster_members.csv")

GEO_KM = 42164

for pos_id in [271, 319]:
    subset = cm[cm['pos_idx'] == pos_id].copy().reset_index(drop=True)
    coords = SkyCoord(ra=subset['src_ra_deg'].values*u.deg,
                      dec=subset['src_dec_deg'].values*u.deg)
    
    print(f"\n{'='*60}")
    print(f"POS {pos_id} — cross-date pairs by separation threshold")
    print(f"{'='*60}")
    
    for threshold in [200, 400, 600]:
        pairs = []
        for i, j in combinations(range(len(subset)), 2):
            if subset['obs_date'].iloc[i] == subset['obs_date'].iloc[j]:
                continue
            sep = coords[i].separation(coords[j]).arcsec
            if sep < threshold:
                dt = abs(pd.Timestamp(subset['obs_date'].iloc[i]) - 
                         pd.Timestamp(subset['obs_date'].iloc[j])).days
                dmag = abs(subset['plate_mag'].iloc[i] - subset['plate_mag'].iloc[j])
                pairs.append((sep, dt, dmag,
                             subset['obs_date'].iloc[i], subset['src_index'].iloc[i],
                             subset['obs_date'].iloc[j], subset['src_index'].iloc[j]))
        
        print(f"\n  < {threshold} arcsec ({GEO_KM * np.tan(np.radians(threshold/3600)):.0f} km): {len(pairs)} pairs")
        for sep, dt, dmag, d1, s1, d2, s2 in sorted(pairs):
            print(f"    {d1} s{s1} <-> {d2} s{s2}: {sep:.0f}\" {sep*GEO_KM/206265:.0f}km Δt={dt}d Δmag={dmag:.1f}")