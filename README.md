# Railway Rainfall Model (MATLAB)

This repository contains the MATLAB scripts and data used in the dissertation project:

> **Rail-ising the Potential of UK Infrastructure Stocks: Quantifying the Material Stock and Assessing Climate Vulnerability of the UK Railway System**  
> MSc Cognitive and Computational Neuroscience, University of Sheffield  
> Author: Yuhou He  
> Supervisor: Dr. Hadi Arbabi

The project applies **Material Stock and Flow Analysis (MSFA)**, climate exposure mapping, and fragility function modelling to assess the resilience of the UK national railway network under extreme rainfall scenarios.

---

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

---

## Requirements

- **MATLAB** (tested with R2024a or later)
- Mapping Toolbox
- Statistics and Machine Learning Toolbox

---

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/heyuhou/Railway-rainfall-model-Matlab.git
