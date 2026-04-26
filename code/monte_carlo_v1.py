"""
Monte Carlo permutation test for spatial recurrence of VASCO transients.

Reads:  data/cluster_members.csv
Outputs: p-value for cross-date close pair count vs circular null model (10k iterations, 2.0 deg radius)
"""
import os
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from itertools import combinations

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
CLUSTER_FILE = os.path.join(DATA_DIR, "cluster_members.csv")

cm = pd.read_csv(CLUSTER_FILE)

GEO_KM = 42164
THRESHOLD = 400  # arcsec


def count_cross_date_pairs(subset, threshold_arcsec):
    coords = SkyCoord(ra=subset['src_ra_deg'].values * u.deg,
                      dec=subset['src_dec_deg'].values * u.deg)
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

    field_radius_deg = 2.0
    ra_center = subset['src_ra_deg'].mean()
    dec_center = subset['src_dec_deg'].mean()

    observed = count_cross_date_pairs(subset, THRESHOLD)

    null_counts = []
    for _ in range(N_ITER):
        fake = subset.copy()
        r = field_radius_deg * np.sqrt(rng.uniform(0, 1, len(subset)))
        theta = rng.uniform(0, 2 * np.pi, len(subset))
        fake['src_ra_deg'] = ra_center + r * np.cos(theta) / np.cos(np.radians(dec_center))
        fake['src_dec_deg'] = dec_center + r * np.sin(theta)
        null_counts.append(count_cross_date_pairs(fake, THRESHOLD))

    null_counts = np.array(null_counts)
    p_value = (null_counts >= observed).mean()

    print(f"Pos {pos_id}: observed={observed} pairs, "
          f"null mean={null_counts.mean():.2f} ± {null_counts.std():.2f}, "
          f"p={p_value:.4f}")
