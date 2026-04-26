"""
Basic dataset exploration for the VASCO transients catalog.

Reads:  data/transients_all.csv
Outputs: shape, column names, date range, RA/Dec range, kp_bin distribution,
         nuclear_window counts, JD-to-ISO samples, top dates by detection count.
"""
import os
import pandas as pd
import numpy as np
from astropy.time import Time

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
TRANSIENTS_FILE = os.path.join(DATA_DIR, "transients_all.csv")

df = pd.read_csv(TRANSIENTS_FILE)

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

print("\n=== JD SAMPLES ===")
for jd in df['obs_jd'].head(10):
    t = Time(jd, format='jd')
    print(f"JD {jd:.6f} -> {t.iso}")

print("\n=== TOP DATES BY DETECTION COUNT ===")
date_counts = df.groupby('obs_date').size().sort_values(ascending=False)
print(date_counts.head(20))
