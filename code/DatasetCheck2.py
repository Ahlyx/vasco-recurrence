import pandas as pd
import numpy as np

df = pd.read_csv("transients_all.csv")

# Pull the highest count date
top_date = "1950-05-23"
day = df[df['obs_date'] == top_date][['ra', 'dec', 'obs_jd', 'mag']].sort_values('obs_jd')
print(f"Detections on {top_date}: {len(day)}")
print(day.to_string())