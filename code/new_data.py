"""
Find top candidate positions for requery by counting same-date multi-source events.

Reads:  data/transients_all.csv
Outputs: ranked list of (pos_idx, obs_date, count) for positions with >1 transient on the same date.
"""
import os
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
TRANSIENTS_FILE = os.path.join(DATA_DIR, "transients_all.csv")

df = pd.read_csv(TRANSIENTS_FILE)

same_date = df.groupby(['pos_idx', 'obs_date']).size()
multi = same_date[same_date > 1].reset_index()
multi.columns = ['pos_idx', 'obs_date', 'count']
print(multi.sort_values('count', ascending=False).head(20))
