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
At 400 arcsec threshold, positions 11 and 270 showed 100% source involvement
rates indicating density-driven pair inflation. At 100 arcsec threshold the
signal persists strongly (pos 11: 76 pairs, null mean 1.1, p=0.000000 at all
radii; pos 270: 83 pairs, null mean 1.5, p=0.000000 at all radii), but source
involvement remains high (82/112 at pos 11, 91/135 at pos 270), indicating the
bulk statistics are not yet clean. Individual noteworthy pairs exist — pos 270
s2173/s2200 at 29 arcsec across 1,037 days dmag=1.2 has no known object within
44 arcsec — but one pair (s9785/s9827, 33 arcsec) has a known BY Draconis
variable star at 12.93 arcsec offset. These positions are flagged as requiring
density-corrected analysis and individual SIMBAD verification before any bulk
claims can be made. The 5.2σ result rests solely on positions 271 and 319.

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
6. `quantify.py` — close pair identification by threshold (200/400/600 arcsec)
7. `quantify2.py` — pair counts by threshold + Mann-Whitney dmag test
8. `physical_check.py` — source involvement at pos 271 and 319
9. `monte_carlo_v1.py` — initial Monte Carlo (10k iter, SkyCoord, circular)
10. `monte_carlo_v2.py` — final Monte Carlo (100k iter, 3 field radii,
    pair count + unique source count)
11. `fishers_method_v2.py` — Fisher combined p-value for pos 271 and 319
12. `new_data.py` — identify top candidate positions for extended requery
13. `sky_coords.py` — RA/Dec lookup for candidate positions
14. `extended_requery.py` — daschlab requery for new positions
15. `extended_recurrence.py` — recurrence analysis and Monte Carlo on extended dataset
16. `combinatorics_test.py` — source density / involvement check at 400 arcsec
17. `source_involvement.py` — source involvement at 100 arcsec for pos 11 and 270
18. `gaia_check.py` — proper motion drift check for a candidate Gaia source
19. `fishers_method.py` — Fisher combined p-value across all five positions

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