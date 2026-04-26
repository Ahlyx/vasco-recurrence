import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from itertools import combinations

cm = pd.read_csv(r"DASCH_Kp_Replication\dasch_northern_reproduction_v2\dasch_northern_reproduction\cluster_requery\cluster_members.csv")

GEO_KM = 42164
THRESHOLD = 400
N_ITER = 100000

def count_cross_date_pairs(subset, threshold_arcsec):
    ra = subset['src_ra_deg'].values
    dec = subset['src_dec_deg'].values
    dates = subset['obs_date'].values
    count = 0
    for i in range(len(ra)):
        for j in range(i+1, len(ra)):
            if dates[i] == dates[j]:
                continue
            dra = (ra[i] - ra[j]) * np.cos(np.radians((dec[i]+dec[j])/2))
            ddec = dec[i] - dec[j]
            sep = np.sqrt(dra**2 + ddec**2) * 3600
            if sep < threshold_arcsec:
                count += 1
    return count

rng = np.random.default_rng(42)

for pos_id in [271, 319]:
    subset = cm[cm['pos_idx'] == pos_id].copy().reset_index(drop=True)
    ra_center = subset['src_ra_deg'].mean()
    dec_center = subset['src_dec_deg'].mean()
    observed = count_cross_date_pairs(subset, THRESHOLD)

    for field_radius_deg in [2.0, 1.0, 0.5]:
        null_counts = []
        for _ in range(N_ITER):
            fake = subset.copy()
            fake['src_ra_deg'] = ra_center + rng.uniform(-field_radius_deg, field_radius_deg, len(subset))
            fake['src_dec_deg'] = dec_center + rng.uniform(-field_radius_deg, field_radius_deg, len(subset))
            null_counts.append(count_cross_date_pairs(fake, THRESHOLD))
        null_counts = np.array(null_counts)
        p = (null_counts >= observed).mean()
        print(f"Pos {pos_id} | radius={field_radius_deg}deg | observed={observed} | "
              f"null={null_counts.mean():.3f}±{null_counts.std():.3f} | p={p:.6f} | N={N_ITER}")