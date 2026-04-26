"""
Print all columns and rows of cluster members for visual inspection.

Reads:  data/cluster_members.csv
Outputs: shape, column list, full table.
"""
import os
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
CLUSTER_FILE = os.path.join(DATA_DIR, "cluster_members.csv")

cm = pd.read_csv(CLUSTER_FILE)
print(cm.shape)
print(cm.columns.tolist())
print(cm.to_string())
