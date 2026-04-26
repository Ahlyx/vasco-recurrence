"""
Inspect detections on the highest-count date in the VASCO transients catalog.

Reads:  data/transients_all.csv
Outputs: sorted listing of all detections on 1950-05-23 (ra, dec, obs_jd, mag).
"""
import os
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
TRANSIENTS_FILE = os.path.join(DATA_DIR, "transients_all.csv")

df = pd.read_csv(TRANSIENTS_FILE)

top_date = "1950-05-23"
day = df[df['obs_date'] == top_date][['ra', 'dec', 'obs_jd', 'mag']].sort_values('obs_jd')
print(f"Detections on {top_date}: {len(day)}")
print(day.to_string())
