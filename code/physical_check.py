import pandas as pd
import numpy as np
import os
from itertools import combinations
from collections import Counter

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
cm = pd.read_csv(os.path.join(DATA_DIR, "cluster_members.csv"))

THRESHOLD = 400
GEO_KM = 42164

for pos_id in [271, 319]:
    subset = cm[cm['pos_idx'] == pos_id].copy().reset_index(drop=True)
    ra = subset['src_ra_deg'].values
    dec = subset['src_dec_deg'].values
    dates = subset['obs_date'].values
    
    source_pair_counts = Counter()
    
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
    
    n_sources = len(subset)
    n_involved = len(source_pair_counts)
    print(f"\nPos {pos_id}: {n_involved}/{n_sources} sources involved in close pairs")
    print(f"Top sources by pair count:")
    for idx, count in source_pair_counts.most_common(5):
        print(f"  {dates[idx]} src{subset['src_index'].iloc[idx]}: "
              f"{count} pairs | RA={ra[idx]:.4f} Dec={dec[idx]:.4f} "
              f"mag={subset['plate_mag'].iloc[idx]:.2f}")