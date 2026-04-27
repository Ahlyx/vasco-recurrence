"""
Cross-date recurrence analysis and Monte Carlo for extended positions.

Reads:  data/extended_cluster_members.csv
Outputs: pair counts, top pairs, and Monte Carlo p-values at three field radii
         (0.5, 1.0, 2.0 deg, 100k iterations) for each position in the dataset.
"""
import pandas as pd
import numpy as np
from itertools import combinations
import os

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
df = pd.read_csv(os.path.join(DATA_DIR, "extended_cluster_members.csv"))

GEO_KM = 42164
THRESHOLD = 100  # arcsec

def cross_date_pairs(subset, threshold):
    ra = subset['src_ra_deg'].values
    dec = subset['src_dec_deg'].values
    dates = subset['obs_date'].values
    mags = subset['plate_mag'].values
    pairs = []
    for i in range(len(ra)):
        for j in range(i+1, len(ra)):
            if dates[i] == dates[j]:
                continue
            dra = (ra[i]-ra[j]) * np.cos(np.radians((dec[i]+dec[j])/2))
            ddec = dec[i]-dec[j]
            sep = np.sqrt(dra**2 + ddec**2) * 3600
            if sep < threshold:
                dt = abs(pd.Timestamp(dates[i]) - pd.Timestamp(dates[j])).days
                dmag = abs(mags[i]-mags[j])
                pairs.append((sep, dt, dmag, dates[i], subset['src_index'].iloc[i],
                             dates[j], subset['src_index'].iloc[j]))
    return pairs

def monte_carlo(subset, threshold, n_iter=100000):
    ra_center = subset['src_ra_deg'].mean()
    dec_center = subset['src_dec_deg'].mean()
    observed = len(cross_date_pairs(subset, threshold))
    rng = np.random.default_rng(42)
    null_counts = []
    for field_radius in [2.0, 1.0, 0.5]:
        null = []
        for _ in range(n_iter):
            fake = subset.copy()
            r = field_radius * np.sqrt(rng.uniform(0, 1, len(subset)))
            theta = rng.uniform(0, 2*np.pi, len(subset))
            fake['src_ra_deg'] = ra_center + r*np.cos(theta)/np.cos(np.radians(dec_center))
            fake['src_dec_deg'] = dec_center + r*np.sin(theta)
            null.append(len(cross_date_pairs(fake, threshold)))
        null = np.array(null)
        p = (null >= observed).mean()
        print(f"  r={field_radius}deg | observed={observed} | null={null.mean():.3f}±{null.std():.3f} | p={p:.6f}")

print(f"Extended dataset: {len(df)} sources across {df['pos_idx'].nunique()} positions\n")

for pos_id in sorted(df['pos_idx'].unique()):
    subset = df[df['pos_idx'] == pos_id].copy().reset_index(drop=True)
    pairs = cross_date_pairs(subset, THRESHOLD)
    print(f"{'='*60}")
    print(f"POS {pos_id} | RA={subset['grid_ra'].iloc[0]:.1f} Dec={subset['grid_dec'].iloc[0]:.2f}")
    print(f"  {len(subset)} sources | {subset['obs_date'].nunique()} dates | {len(pairs)} cross-date pairs <{THRESHOLD}\"")
    
    if pairs:
        print(f"  Top pairs by separation:")
        for sep, dt, dmag, d1, s1, d2, s2 in sorted(pairs)[:5]:
            print(f"    {d1} s{s1} <-> {d2} s{s2}: {sep:.0f}\" {sep*GEO_KM/206265:.0f}km Δt={dt}d Δmag={dmag:.1f}")
        
        print(f"\n  Monte Carlo (N=100k):")
        monte_carlo(subset, THRESHOLD)
    print()