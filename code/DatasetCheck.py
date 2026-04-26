# explore.py
import pandas as pd
import numpy as np

df = pd.read_csv("transients_all.csv")

print("=== SHAPE ===")
print(df.shape)

print("\n=== COLUMNS ===")
print(df.columns.tolist())

print("\n=== FIRST 5 ROWS ===")
print(df.head())

print("\n=== NULLS ===")
print(df.isnull().sum())

print("\n=== DATE RANGE ===")
print(f"obs_jd min: {df['obs_jd'].min()}")
print(f"obs_jd max: {df['obs_jd'].max()}")

print("\n=== RA/DEC RANGE ===")
print(f"RA:  {df['ra'].min():.2f} to {df['ra'].max():.2f}")
print(f"Dec: {df['dec'].min():.2f} to {df['dec'].max():.2f}")

print("\n=== KP BIN DISTRIBUTION ===")
print(df['kp_bin'].value_counts().sort_index())

print("\n=== NUCLEAR WINDOW ===")
print(df['nuclear_window'].value_counts())

# add to explore.py or run interactively
from astropy.time import Time

# Check actual time precision on a few
sample_jd = [2432922.341765, 2433098.0, 2433151.0]
for jd in df['obs_jd'].head(10):
    t = Time(jd, format='jd')
    print(f"JD {jd:.6f} -> {t.iso}")

date_counts = df.groupby('obs_date').size().sort_values(ascending=False)
print(date_counts.head(20))