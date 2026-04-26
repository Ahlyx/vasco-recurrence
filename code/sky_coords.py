"""
RA/Dec summary for target grid positions with transient counts.

Reads:  data/transients_all.csv
Outputs: for each target pos_idx, prints RA, Dec, total transient count, and unique date count.
"""
import os
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
TRANSIENTS_FILE = os.path.join(DATA_DIR, "transients_all.csv")

df = pd.read_csv(TRANSIENTS_FILE)

target_positions = [144, 270, 168, 11, 169, 269, 143]

for pos in target_positions:
    subset = df[df['pos_idx'] == pos][['ra', 'dec']].iloc[0]
    dates = df[df['pos_idx'] == pos]['obs_date'].nunique()
    total = len(df[df['pos_idx'] == pos])
    print(f"Pos {pos}: RA={subset['ra']:.1f} Dec={subset['dec']:.1f} | {total} transients across {dates} dates")
