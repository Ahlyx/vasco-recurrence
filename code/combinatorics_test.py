"""
Tests whether close pair counts are driven by source density at extended positions.

Reads:  data/extended_cluster_members.csv
Outputs: per-position source involvement counts and top sources by pair count
         at 400 arcsec threshold for pos 11 and 270.
"""
import pandas as pd
import numpy as np
import os
from itertools import combinations
from collections import Counter

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
df = pd.read_csv(os.path.join(DATA_DIR, "extended_cluster_members.csv"))

THRESHOLD = 400

for pos_id in [11, 270]:
    subset = df[df['pos_idx'] == pos_id].copy().reset_index(drop=True)
    ra = subset['src_ra_deg'].values
    dec = subset['src_dec_deg'].values
    dates = subset['obs_date'].values
    
    source_pair_counts = Counter()
    pair_dates = []
    
    for i in range(len(ra)):
        for j in range(i+1, len(ra)):
            if dates[i] == dates[j]:
                continue
            dra = (ra[i]-ra[j]) * np.cos(np.radians((dec[i]+dec[j])/2))
            ddec = dec[i]-dec[j]
            sep = np.sqrt(dra**2 + ddec**2) * 3600
            if sep < THRESHOLD:
                source_pair_counts[i] += 1
                source_pair_counts[j] += 1
                pair_dates.append((dates[i], dates[j]))
    
    print(f"\nPos {pos_id}: {len(source_pair_counts)} unique sources involved in close pairs")
    print(f"Top 10 sources by pair count:")
    for idx, count in source_pair_counts.most_common(10):
        print(f"  src {subset['src_index'].iloc[idx]} ({dates[idx]}): {count} close pairs, "
              f"RA={ra[idx]:.4f} Dec={dec[idx]:.4f} mag={subset['plate_mag'].iloc[idx]:.2f}")
    
    # How many dates are involved in cross-date pairs?
    all_dates = set([d for pair in pair_dates for d in pair])
    print(f"Dates involved in close pairs: {sorted(all_dates)}")