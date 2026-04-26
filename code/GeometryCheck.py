import pandas as pd
cm = pd.read_csv(r"DASCH_Kp_Replication\dasch_northern_reproduction_v2\dasch_northern_reproduction\cluster_requery\cluster_members.csv")
print(cm.shape)
print(cm.columns.tolist())
print(cm.to_string())