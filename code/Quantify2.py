"""
Cross-date pair proximity table with Mann-Whitney delta-mag test for VASCO cluster members.

Reads:  data/cluster_members.csv
Outputs: pair counts at separation thresholds 200/400/600 arcsec; Mann-Whitney U test
         comparing delta-mag of close (<400 arcsec) vs far (>=400 arcsec) cross-date pairs.
         One-sided test: H1 = close pairs have lower delta-mag (more consistent brightness).
"""
import os
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from itertools import combinations
from scipy.stats import mannwhitneyu

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
CLUSTER_FILE = os.path.join(DATA_DIR, "cluster_members.csv")

cm = pd.read_csv(CLUSTER_FILE)

GEO_KM = 42164
CLOSE_THRESHOLD = 400  # arcsec — boundary between close and far pairs

for pos_id in [271, 319]:
    subset = cm[cm['pos_idx'] == pos_id].copy().reset_index(drop=True)
    coords = SkyCoord(ra=subset['src_ra_deg'].values * u.deg,
                      dec=subset['src_dec_deg'].values * u.deg)

    print(f"\n{'='*60}")
    print(f"POS {pos_id} - cross-date pairs by separation threshold")
    print(f"{'='*60}")

    # Collect ALL cross-date pairs for the Mann-Whitney test
    all_pairs = []
    for i, j in combinations(range(len(subset)), 2):
        if subset['obs_date'].iloc[i] == subset['obs_date'].iloc[j]:
            continue
        sep = coords[i].separation(coords[j]).arcsec
        dt = abs(pd.Timestamp(subset['obs_date'].iloc[i]) -
                 pd.Timestamp(subset['obs_date'].iloc[j])).days
        dmag = abs(subset['plate_mag'].iloc[i] - subset['plate_mag'].iloc[j])
        all_pairs.append((sep, dt, dmag,
                          subset['obs_date'].iloc[i], subset['src_index'].iloc[i],
                          subset['obs_date'].iloc[j], subset['src_index'].iloc[j]))

    # Threshold-based pair tables
    for threshold in [200, 400, 600]:
        pairs = [(sep, dt, dmag, d1, s1, d2, s2)
                 for sep, dt, dmag, d1, s1, d2, s2 in all_pairs if sep < threshold]
        print(f"\n  < {threshold} arcsec ({GEO_KM * np.tan(np.radians(threshold/3600)):.0f} km): {len(pairs)} pairs")
        for sep, dt, dmag, d1, s1, d2, s2 in sorted(pairs):
            print(f"    {d1} s{s1} <-> {d2} s{s2}: {sep:.0f}\" {sep*GEO_KM/206265:.0f}km dt={dt}d dmag={dmag:.1f}")

    # Mann-Whitney U: delta-mag of close (<400 arcsec) vs far (>=400 arcsec) cross-date pairs
    print(f"\n  --- Mann-Whitney U: delta-mag close vs far (split at {CLOSE_THRESHOLD}\") ---")
    close_dmag = [dmag for sep, dt, dmag, *_ in all_pairs if sep < CLOSE_THRESHOLD]
    far_dmag   = [dmag for sep, dt, dmag, *_ in all_pairs if sep >= CLOSE_THRESHOLD]

    if len(close_dmag) >= 2 and len(far_dmag) >= 2:
        # One-sided: alternative='less' tests H1: close delta-mag < far delta-mag
        stat, p = mannwhitneyu(close_dmag, far_dmag, alternative='less')
        print(f"  close pairs (n={len(close_dmag)}): mean dmag={np.mean(close_dmag):.3f} +/- {np.std(close_dmag):.3f}")
        print(f"  far pairs   (n={len(far_dmag)}):   mean dmag={np.mean(far_dmag):.3f} +/- {np.std(far_dmag):.3f}")
        print(f"  U={stat:.1f}, p={p:.4f} (one-sided: H1 = close < far)")
    else:
        print(f"  Insufficient pairs for test: close n={len(close_dmag)}, far n={len(far_dmag)}")
