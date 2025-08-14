## Repository Contents

| File / Folder             | Description |
|---------------------------|-------------|
| `Railway_track.cpg`, `.dbf`, `.prj`, `.qmd` | GIS shapefile components for UK railway track data. |
| `euro_pr.tif`              | Climate raster data (projected rainfall, from IPCC Interactive Atlas). |
| `rail_with_rain.csv`       | Merged dataset of railway segments with corresponding rainfall exposure values. |
| `trackpr.m`                | Main MATLAB script for processing railway and rainfall data, and running analyses. |

---

## Methods Implemented

The repository includes MATLAB code for:

1. **Climate Exposure Mapping**  
   Mapping projected rainfall data to UK railway track segments using GIS midpoint matching.

2. **Fragility Function Development**  
   Creating rainfall–failure probability curves (logistic model).

3. **Monte Carlo Simulation of Failure Frequency**  
   Simulating network failure events under varying rainfall intensities.

4. **Failure Scenario Simulation (Logistic)**  
   Modelling rainfall-driven disruption scenarios.

5. **Forecasting Material Needs under Adaptation Scenarios**  
   Estimating additional material demand (track, ballast, etc.) under climate adaptation strategies.
