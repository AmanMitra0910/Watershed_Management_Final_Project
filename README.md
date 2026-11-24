# Groundwater Recharge Potential and Recharge-Structure Prioritisation: Upper Bhima-type Sub-watershed (Synthetic Case Study)

AHP + weighted-overlay assessment of groundwater recharge potential, with rule-based siting and ranking of check dams, percolation tanks and farm ponds, and a consistency check against observation wells.

**IMPORTANT: all data are synthetic.** The study area is a hypothetical 100 km² (10 x 10 km, 100 m cells) Upper Bhima-type sub-watershed in Maharashtra. No real IMD, Bhuvan, CGWB or Sentinel-2 data were used. Results illustrate the workflow and must not be used for real siting decisions.

## Structure
```
code/gwrp_analysis.py               full pipeline (data generation -> AHP -> overlay -> siting -> figures)
data/synthetic_wells.csv            28 synthetic observation wells
results/Fig1_recharge_potential_map.png
results/Fig2_well_validation.png
results/ranked_recharge_sites.csv   15 ranked structures (type, location, score, priority class)
results/run_log.txt                 console output (weights, CR, zone areas, Spearman rho)
docs/writeup.md                     project write-up
requirements.txt  LICENSE  .gitignore
```

## Run
```
git clone <this-repo-url>
cd upper-bhima-groundwater-recharge
pip install -r requirements.txt
python code/gwrp_analysis.py
```
Outputs are written to `data/` and `results/`. Seeds are fixed, so results are reproducible.

## Method summary
- Layers (scored 1-5): geology/geomorphology, slope, soil, lineament density, rainfall, drainage density, LULC.
- AHP weights (principal eigenvector): 34.7, 23.2, 14.8, 9.0, 9.0, 5.6, 3.7 %; CR = 0.016.
- GWRPI = sum(w_i * r_i); five equal-interval zones.
- Priority score = 0.5 GWRPI + 0.3 runoff availability + 0.2 groundwater need (min-max normalised).
- Validation: Spearman correlation of GWRPI vs. well water-level rise (synthetic wells, so partly circular).

## Using real data
Replace the synthetic layer-generation block in `code/gwrp_analysis.py` with rasters resampled to a common grid (SRTM/Cartosat DEM, Sentinel-2 LULC, Bhuvan soil/geology/lineaments, IMD rainfall, CGWB/India-WRIS wells). The reclassification, AHP, overlay, siting and plotting steps apply unchanged.

## References
Sahu et al. (2022); Jadhav et al. (2022); Saaty (1980).
