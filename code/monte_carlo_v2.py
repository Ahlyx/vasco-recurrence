"""
Monte Carlo permutation test with multiple field radii for VASCO transient spatial recurrence.

Reads:  data/cluster_members.csv
Outputs: p-values for cross-date close pair counts and unique source counts across field radii
         [0.5, 1.0, 2.0 deg] at 100k iterations, for positions 271 and 319.
"""
import os
import pandas as pd
import numpy as np
from itertools import combinations

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
CLUSTER_FILE = os.path.join(DATA_DIR, "cluster_members.csv")

cm = pd.read_csv(CLUSTER_FILE)

GEO_KM = 42164
THRESHOLD = 400
N_ITER = 100000

# Field radii tested and their physical interpretation:
#   2.0 deg — conservative bound: DASCH grid cells are ~3-4 deg across, so 2.0 deg is roughly
#              half a cell and represents the maximum plausible randomization region.
#   1.0 deg — moderately conservative: half of the 2.0 deg bound, a tighter region within the
#              same grid cell.
#   0.5 deg — most conservative: quarter of the cell width. Clustering at this radius is the
#              strongest signal because the null hypothesis has the least area to produce pairs.


def count_cross_date_pairs(subset, threshold_arcsec):
    ra = subset['src_ra_deg'].values
    dec = subset['src_dec_deg'].values
    dates = subset['obs_date'].values
    count = 0
    for i in range(len(ra)):
        for j in range(i + 1, len(ra)):
            if dates[i] == dates[j]:
                continue
            dra = (ra[i] - ra[j]) * np.cos(np.radians((dec[i] + dec[j]) / 2))
            ddec = dec[i] - dec[j]
            sep = np.sqrt(dra ** 2 + ddec ** 2) * 3600
            if sep < threshold_arcsec:
                count += 1
    return count


def count_unique_sources_in_close_pairs(subset, threshold_arcsec):
    """Count unique row indices involved in at least one cross-date close pair.

    Uses unique source count rather than pair count to avoid overcounting from
    non-independent pairs (one source can belong to many pairs).
    """
    ra = subset['src_ra_deg'].values
    dec = subset['src_dec_deg'].values
    dates = subset['obs_date'].values
    involved = set()
    for i in range(len(ra)):
        for j in range(i + 1, len(ra)):
            if dates[i] == dates[j]:
                continue
            dra = (ra[i] - ra[j]) * np.cos(np.radians((dec[i] + dec[j]) / 2))
            ddec = dec[i] - dec[j]
            sep = np.sqrt(dra ** 2 + ddec ** 2) * 3600
            if sep < threshold_arcsec:
                involved.add(i)
                involved.add(j)
    return len(involved)


rng = np.random.default_rng(42)

for pos_id in [271, 319]:
    subset = cm[cm['pos_idx'] == pos_id].copy().reset_index(drop=True)
    ra_center = subset['src_ra_deg'].mean()
    dec_center = subset['src_dec_deg'].mean()
    observed_pairs = count_cross_date_pairs(subset, THRESHOLD)
    observed_unique = count_unique_sources_in_close_pairs(subset, THRESHOLD)

    for field_radius_deg in [2.0, 1.0, 0.5]:
        null_pairs = []
        null_unique = []
        for _ in range(N_ITER):
            fake = subset.copy()
            r = field_radius_deg * np.sqrt(rng.uniform(0, 1, len(subset)))
            theta = rng.uniform(0, 2 * np.pi, len(subset))
            fake['src_ra_deg'] = ra_center + r * np.cos(theta) / np.cos(np.radians(dec_center))
            fake['src_dec_deg'] = dec_center + r * np.sin(theta)
            null_pairs.append(count_cross_date_pairs(fake, THRESHOLD))
            null_unique.append(count_unique_sources_in_close_pairs(fake, THRESHOLD))

        null_pairs = np.array(null_pairs)
        null_unique = np.array(null_unique)
        p_pairs = (null_pairs >= observed_pairs).mean()
        p_unique = (null_unique >= observed_unique).mean()
        print(
            f"Pos {pos_id} | radius={field_radius_deg}deg | "
            f"pairs: obs={observed_pairs} null={null_pairs.mean():.3f}±{null_pairs.std():.3f} p={p_pairs:.6f} | "
            f"unique_srcs: obs={observed_unique} null={null_unique.mean():.3f}±{null_unique.std():.3f} p={p_unique:.6f} | "
            f"N={N_ITER}"
        )
