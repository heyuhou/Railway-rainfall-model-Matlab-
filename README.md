# Railway Rainfall Model (MATLAB)

The code used in this research has been modified by **Li et al. (2024)** and **Zhivomirov (2025)**.  
Li, Q., Punzo, G., Robson, C., Arbabi, H. and Mayfield, M. (2024) 'A systematic approach to climate resilience assessment of infrastructure networks', IEEE Systems Journal, 18(1), pp. 24–35. Available at: https://doi.org/10.1109/JSYST.2023.3329765.
Zhivomirov, H. (2025) Monte Carlo Estimation Examples with Matlab. MATLAB Central File Exchange. Available at: https://www.mathworks.com/matlabcentral/fileexchange/55306-monte-carlo-estimation-examples-with-matlab (Accessed: 6 July 2025).
---

## Repository Contents

| File / Folder             | Description |
|---------------------------|-------------|
| `Railway_track.cpg`, `.dbf`, `.prj`, `.qmd` | GIS shapefile components for UK railway track data. RailEasyUK (2022) Railway GIS data. Available at: https://github.com/raileasyuk/railway-gis-data (Accessed: 24 June 2025). |
| `euro_pr.tif`              | Climate raster data (projected rainfall, from IPCC Interactive Atlas). IPCC (n.d.) Interactive Atlas. Available at: https://interactive-atlas.ipcc.ch/ (Accessed: 11 July 2025). |
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
   Estimating additional material demand under climate adaptation strategies.
