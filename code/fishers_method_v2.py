"""
Fisher's combined probability test for pos 271 and 319 only.

Reads:  hardcoded p-values from monte_carlo_v2.py at 1-degree field radius
Outputs: combined chi-squared statistic, p-value, and sigma equivalent.
"""
import numpy as np
from scipy.stats import chi2, norm

# Only pos 271 and 319 — original pre-selected positions
p_values = [0.000010, 0.000520]
positions = [271, 319]

p_floor = 1 / 100000
p_values_safe = [max(p, p_floor) for p in p_values]

chi2_stat = -2 * sum(np.log(p) for p in p_values_safe)
df_val = 2 * len(p_values)
p_combined = 1 - chi2.cdf(chi2_stat, df_val)
sigma = norm.ppf(1 - p_combined) if p_combined > 0 else float('inf')
print(f"Pos 271+319 only: chi2={chi2_stat:.2f}, df={df_val}, p={p_combined:.2e}, sigma={sigma:.1f}")
