"""
Proper motion drift calculation for a Gaia source over 75 years.

Reads:  hardcoded proper motion values (pm_ra=-7.029 mas/yr, pm_dec=3.175 mas/yr)
Outputs: RA drift, Dec drift, and total drift in arcsec over the plate epoch baseline.
"""
pm_ra = -7.029  # mas/yr
pm_dec = 3.175  # mas/yr
years = 75

drift_ra = pm_ra * years / 1000  # convert to arcsec
drift_dec = pm_dec * years / 1000

total_drift = (drift_ra**2 + drift_dec**2)**0.5

print(f"RA drift:    {drift_ra:.3f} arcsec")
print(f"Dec drift:   {drift_dec:.3f} arcsec")
print(f"Total drift: {total_drift:.3f} arcsec")
