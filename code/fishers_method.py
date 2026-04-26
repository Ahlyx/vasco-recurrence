"""
Fisher's combined probability test across all five target positions.

Reads:  hardcoded p-values from monte_carlo_v2.py runs at 1-degree field radius
Outputs: combined chi-squared statistic, p-value, and sigma equivalent.
"""
import numpy as np
from scipy.stats import chi2, norm

# p-values at 1 degree radius (middle ground), from monte_carlo_v2.py
p_values = [0.000000, 0.000820, 0.000000, 0.000010, 0.000520]
positions = [11, 269, 270, 271, 319]

# Use conservative floor of 1/N_ITER for zeros
p_floor = 1 / 100000
p_values_safe = [max(p, p_floor) for p in p_values]

chi2_stat = -2 * sum(np.log(p) for p in p_values_safe)
df = 2 * len(p_values)
p_combined = 1 - chi2.cdf(chi2_stat, df)

sigma = norm.ppf(1 - p_combined) if p_combined > 0 else float('inf')
print(f"Fisher combined: chi2={chi2_stat:.2f}, df={df}, p={p_combined:.2e}, sigma={sigma:.1f}")
