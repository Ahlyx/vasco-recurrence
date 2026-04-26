# VASCO Transient Recurrence Analysis

Independent analysis of spatial recurrence in VASCO photographic plate 
transient candidates from the DASCH northern hemisphere scan (Cann 2026).

## Finding
Cross-date spatial clustering of transient candidates at two independent 
grid positions (271 and 319) is statistically significant under Monte Carlo 
permutation testing (p<0.0001 at 2deg field radius, p<0.011 at 0.5deg).
Sources separated by up to 1,143 days recur within 27-120 km at GEO altitude,
inconsistent with random artifact hypothesis.

## Data
- `data/cluster_members.csv` — per-source coordinates from Cann (2026c) 
  cluster requery. 19 sources across 2 grid positions and 9 cluster events.
- `data/transients_all.csv` — combined DASCH fast scan slices, 1939 transient
  candidates with grid coordinates and timestamps.

## Code (run in order)
1. `DatasetCheck.py` — initial dataset exploration
2. `DatasetCheck2.py` — top date spatial analysis  
3. `GeometryCheck.py` — cluster_members structure
4. `ScaleCheck.py` — same-date within-cluster separations
5. `RecurrenceCheck.py` — cross-date proximity analysis
6. `Quantify.py` — close pair identification (<200 arcsec)
7. `Quantify2.py` — pair counts by separation threshold
8. `monte_carloV1.py` — initial Monte Carlo (10k iter, SkyCoord)
9. `monte_carloV2.py` — final Monte Carlo (100k iter, 3 field radii)

## Dependencies
pip install astropy pandas numpy scipy matplotlib

## Attribution
Built on data from Cann (2026), Villarroel et al. (2020), and the 
DASCH project at Harvard College Observatory.
