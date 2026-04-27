# VASCO Transient Recurrence Analysis

Independent analysis of spatial recurrence in VASCO photographic plate 
transient candidates from the DASCH northern hemisphere scan (Cann 2026).

## Core Finding
Cross-date spatial clustering of transient candidates at two independent 
grid positions (271 and 319) is statistically significant under Monte Carlo 
permutation testing. Fisher's combined probability across both positions 
yields p=1.04e-07 (5.2σ).

At the 1-degree field radius:
- Pos 271 (RA=302.4, Dec=+3.2): p=0.000010
- Pos 319 (RA=273.6, Dec=+15.8): p=0.000520

Sources separated by up to 1,143 days recur within 27-120 km of each other 
when projected to geosynchronous altitude (42,164 km), inconsistent with the 
random artifact hypothesis. Result is robust across field radii of 0.5, 1.0, 
and 2.0 degrees.

## Extended Analysis
Three additional positions (11, 269, 270) were queried via daschlab and 
subjected to the same recurrence analysis. Positions 11 and 270 showed 
high raw pair counts but also 100% source involvement rates, indicating 
the pairs are density-driven rather than reflecting genuine positional 
recurrence. These positions require further analysis with a tighter 
separation threshold and APASS density correction before any conclusions 
can be drawn. Position 269 had insufficient observation dates (n=2) for 
a meaningful test.

The 5.2σ result rests solely on positions 271 and 319.

## Methodology Notes
- Monte Carlo null model: circular uniform randomization within field radius,
  100,000 iterations, positions randomized independently per iteration
- Separation metric: cosine-corrected planar approximation (valid at <2 deg)
- Three field radii tested: 0.5, 1.0, 2.0 deg (DASCH grid cells ~3-4 deg)
- Pair count and unique source count both reported to avoid overcounting
- Mann-Whitney U test on Δmag (close vs far pairs): p=0.057 pos 271, 
  p=0.241 pos 319 (trend present, not significant at n=3-9 pairs)

## Known Limitations
- Analysis restricted to positions pre-selected by Cann (2026c) for 
  same-date multi-source clustering — not a random sky sample
- cluster_members.csv contains only 19 sources across 2 positions
- Δmag variations (0.2-4.4) across close pairs are not yet explained
- Extended positions (11, 270) show density artifact requiring correction
- No periodicity analysis performed on recurrence intervals

## Data
- `data/cluster_members.csv` — per-source coordinates from Cann (2026c) 
  cluster requery. 19 sources across 2 grid positions and 9 cluster events.
- `data/extended_cluster_members.csv` — per-source coordinates for positions 
  11, 269, 270 from extended requery. 266 sources across 37 cluster events.
- `data/transients_all.csv` — combined DASCH fast scan slices (v2 only), 
  1939 transient candidates with grid coordinates and timestamps.

## Code (run in order)
1. `dataset_check.py` — initial dataset exploration
2. `dataset_check2.py` — top date spatial analysis
3. `geometry_check.py` — cluster_members structure
4. `scale_check.py` — same-date within-cluster separations
5. `recurrence_check.py` — cross-date proximity analysis
6. `Quantify.py` — close pair identification (<200 arcsec)
7. `Quantify2.py` — pair counts by threshold + Mann-Whitney Δmag test
8. `monte_carlo_v1.py` — initial Monte Carlo (10k iter, SkyCoord, circular)
9. `monte_carlo_v2.py` — final Monte Carlo (100k iter, 3 field radii, 
   pair count + unique source count)
10. `new_data.py` — identify top candidate positions for extended requery
11. `sky_coords.py` — RA/Dec lookup for candidate positions
12. `extended_requery.py` — daschlab requery for new positions
13. `extended_recurrence.py` — recurrence analysis on extended dataset
14. `fishers_method.py` — Fisher combined p-value across positions

## Dependencies
pip install astropy pandas numpy scipy matplotlib daschlab

## Attribution
Built on data and methodology from Cann (2026a, 2026b, 2026c), 
Villarroel et al. (2020), and the DASCH project at Harvard College 
Observatory. Analysis conducted independently using publicly available data.

## References
- Cann, K. (2026a). arXiv:2604.04950
- Cann, K. (2026b). arXiv:2604.06234  
- Cann, K. (2026c). ESSOAr:10.22541/essoar.15002100/v1
- Villarroel, B. et al. (2020). Astron. J., 159, 8
- Bruehl, S.P. & Villarroel, B. (2025). Scientific Reports, 15, 34125