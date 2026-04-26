import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from itertools import combinations

cm = pd.read_csv(r"DASCH_Kp_Replication\dasch_northern_reproduction_v2\dasch_northern_reproduction\cluster_requery\cluster_members.csv")

GEO_KM = 42164
THRESHOLD = 400  # arcsec

# For each pos, compute observed close pair count
# Then Monte Carlo: randomize source positions within the field and count pairs

def count_cross_date_pairs(subset, threshold_arcsec):
    coords = SkyCoord(ra=subset['src_ra_deg'].values*u.deg,
                      dec=subset['src_dec_deg'].values*u.deg)
    count = 0
    for i, j in combinations(range(len(subset)), 2):
        if subset['obs_date'].iloc[i] == subset['obs_date'].iloc[j]:
            continue
        if coords[i].separation(coords[j]).arcsec < threshold_arcsec:
            count += 1
    return count

N_ITER = 10000
rng = np.random.default_rng(42)

for pos_id in [271, 319]:
    subset = cm[cm['pos_idx'] == pos_id].copy().reset_index(drop=True)
    
    # Field bounds
    ra_min, ra_max = subset['src_ra_deg'].min(), subset['src_ra_deg'].max()
    dec_min, dec_max = subset['src_dec_deg'].min(), subset['src_dec_deg'].max()
    
    # Pad field to full DASCH grid cell (~2 deg radius)
    field_radius_deg = 2.0
    ra_center = subset['src_ra_deg'].mean()
    dec_center = subset['src_dec_deg'].mean()
    
    observed = count_cross_date_pairs(subset, THRESHOLD)
    
    # Monte Carlo: randomize positions within field, keep dates fixed
    null_counts = []
    for _ in range(N_ITER):
        fake = subset.copy()
        # Random positions within field radius
        fake['src_ra_deg'] = ra_center + rng.uniform(-field_radius_deg, field_radius_deg, len(subset))
        fake['src_dec_deg'] = dec_center + rng.uniform(-field_radius_deg, field_radius_deg, len(subset))
        null_counts.append(count_cross_date_pairs(fake, THRESHOLD))
    
    null_counts = np.array(null_counts)
    p_value = (null_counts >= observed).mean()
    
    print(f"Pos {pos_id}: observed={observed} pairs, "
          f"null mean={null_counts.mean():.2f} ± {null_counts.std():.2f}, "
          f"p={p_value:.4f}")