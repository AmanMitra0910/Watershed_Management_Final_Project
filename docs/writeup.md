# Project write-up (synthetic Upper Bhima-type sub-watershed)

## 1. Problem statement
Groundwater is the main source of irrigation and drinking water in the Upper Bhima basin, but its availability varies spatially with rainfall, terrain, soil, Deccan basalt geology and land use. This project maps groundwater recharge potential in a hypothetical 100 km² Upper Bhima-type sub-watershed (monsoon season, June-September) and ranks sites for check dams, percolation tanks and farm ponds. Stakeholders: GSDA, CGWB, soil and water conservation department, gram panchayats, farmers. Following proposal feedback, siting rules and the well-based consistency check were made explicit.

## 2. Data
All data are synthetic (Python, fixed seed), 100 m cells: DEM (700-900 m), slope and D8 drainage, rainfall (650-1100 mm, west to east), LULC (5 classes), soil (4 texture classes), geology/geomorphology (5 units), lineament density (45 lineaments), and 28 observation wells (pre-monsoon depth, pre- to post-monsoon rise) generated from an independently perturbed weight set plus noise. One monsoon season.

## 3. Methodology
XLRM: X = AHP judgments, reclassification, rainfall variability, classification error; L = structure type, number (5 check dams, 4 percolation tanks, 6 farm ponds) and location; R = AHP weights, weighted overlay, siting rules, priority score; M = zone area shares, Spearman correlation with wells, priority rank.

GWRPI = Σ (w_i × r_i), i = 1..7, where r_i is the 1-5 score of layer i and w_i its AHP weight (Σw_i = 1; Saaty, 1980). CR = [(λmax − n)/(n − 1)]/RI, n = 7, RI = 1.32.

Weights: geology/geomorphology 34.7%, slope 23.2%, soil 14.8%, lineament 9.0%, rainfall 9.0%, drainage density 5.6%, LULC 3.7%; λmax = 7.13, CR = 0.016.

Siting rules: check dams on streams with 0.8-10 km² catchment, slope <10%, at least moderate zone; percolation tanks within 300 m of a stream, slope <5%, high/very high zone, soil and geomorphology score ≥3; farm ponds on agricultural land, slope 1-6%, soil score ≤3, small catchment, away from streams.
Priority score = 0.5 GWRPI + 0.3 runoff availability + 0.2 need (IDW of pre-monsoon depth). Ranks 1-5 Priority 1, 6-10 Priority 2, 11-15 Priority 3.

## 4. Results
Fig. 1 (results/Fig1_recharge_potential_map.png): zones cover very low 5%, low 24%, moderate 35%, high 30%, very high 6%. High-potential belts follow valley-fill and weathered basalt in the west and south-west; very low zones lie on steep ridges in the north-east and south-east. Check dams take ranks 1-4 and 6, percolation tanks 5 and 7-9, farm ponds 10-15.
Fig. 2 (results/Fig2_well_validation.png): Spearman ρ = 0.61 (p < 0.001, n = 28); mean water-level rise ≈3.1 m (low), 4.1 m (moderate), 4.6 m (high zones). Very low and very high zones have one well each.

## 5. Discussion
Results depend mainly on the AHP weights (geology + slope ≈58%), the zone classification and the siting thresholds; fragmented D8 streams in the smoothed synthetic DEM affect check-dam candidates. The well check is partly circular because the wells are synthetic. Not decision-ready: needs real layers, multi-year verified well data, Monte Carlo weight sensitivity, structure sizing, cost-benefit and land-ownership constraints. Farm ponds are mainly water-supply measures with limited recharge.

## 6. References
1. Sahu, S. et al. (2022). Geospatial analysis and AHP for recharge potential zones and structures, Ajani-Jhiri watershed, north Maharashtra.
2. Jadhav et al. (2022). Groundwater potential zones, Kamina sub-watershed, Bhima basin, Maharashtra.
3. Saaty, T. L. (1980). The Analytic Hierarchy Process. McGraw-Hill.
