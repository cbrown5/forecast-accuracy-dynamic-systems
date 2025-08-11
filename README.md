# Assessing predictive accuracy of species abundance models in dynamic systems

Code to support the paper: [Christopher J. Brown, Christina Buelow, Rick D. Stuart-Smith, Neville Barrett, Graham Edgar, Elizabeth Oh2025, "Assessing predictive accuracy of species abundance models in dynamic systems"  Methods in Ecology and Evolution, https://doi.org/10.1111/2041-210X.70105](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.70105). Please cite that paper if using these methods. 

### Data availability 
Data included here are for reproducibility only. If you are wanting to publish original work with this data please see the sources below
Biological data available via https://portal.aodn.org.au/ (Edgar and Stuart-Smith 2014; Edgar and Barrett 2012). The ecological data used for this study are managed through, and were sourced from, Australia’s Integrated Marine Observing System (IMOS) – IMOS is enabled by the National Collaborative Research Infrastructure Strategy (NCRIS). Physical oceanographic data available via https://www.marine.csiro.au/data/trawler/ (Lynch et al. 2014)


## The Motivation
Can we predict change in species abundances when environmental changes are unprecedented in historical data?


## Data Sources  

Biological data from the Australian Temperate Reef Collaboration and Reef Life Survey are available at: https://portal.aodn.org.au/ (Edgar and Stuart-Smith 2014; Edgar and Barrett 2012)
Physical oceanographic data are available at: https://www.marine.csiro.au/data/trawler/ (Lynch et al. 2014)

## Project Structure

You will need to create a folder Structure like this for the scripts to run properly (they need the Outputs/ and figures/ folders to save the intermediate results):

```
PredictabilityMarineSpecies/
├── data/
├── scripts/
├── Outputs/
│   ├── inla-forecasts
│   └── inla-sim-testing/
└── figures/
```

## Data

### Abundance Data (2025-01-24-ATRC-data-for-review.csv)
- Long-term monitoring dataset from sites around Maria Island Marine Reserve and reference areas
- Contains annual abundance surveys from 1992-2021
- Key columns:
  - `Year`: Survey year
  - `Site_code`: Unique site identifier 
  - `Survey_id`: Unique survey identifier
  - `Region`: Location category (Maria-reserve, Reference-north, Reference-mainland, Reference-other)
  - Species abundance columns counts per survey:
    - `Trachinops caudimaculatus` (hulafish)
    - `Tosia australis` (Southern biscuit star)
    - `Haliotis rubra` (Blacklip abalone) 
    - `Jasus edwardsii` (Southern rock lobster)
  - Environmental principal components and lags:
    - `PC1`, `PC2`: First two environmental principal components
    - `PC1_1yrlag`, `PC2_1yrlag`: 1-year lagged PCs
    - `PC1_2yrlag`, `PC2_2yrlag`: 2-year lagged PCs
 - Site indicators
    - isite: numerical site code
    - iregion: numerical region code
    - Year2: copy of year, for INLA use

### Sea Surface Temperature Data (2025-01-24-SST-data-for-review.csv)
Sea surface temperature time series statistics:
Mean, min and max temperature per year, these statistics at 1 and 2 year lags. 

### NMRS Data (2025-01-24-NMRS-data-for-review.csv)
- Long-term monitoring dataset from the Maria Island NMRS
- Contains all data available in the period 1986 to 2023
- Key columns:
  - `Year`: Survey year
  - `Month`: Calendar month
  - `Nitrate_umolL`: Nitrate concentration
  - `Silicate_umolL`: Silicate concentration
  - `Salinity`: Water salinity
  - `Sample Depth`: Water sample depth in metres
  
### Environmental Data (2025-01-10_enviro_PC-timeseries-with-lags.csv)
- Yearly environmental time series from 1989-2023 to be used in model fitting
- Key columns:
  - `Year`: Year of data
  - `Month`: Month of GAM prediction (ie standardized to January)
  - Raw environmental variables:
    - `Nitrate_umolL`: Nitrate concentration
    - `Silicate_umolL`: Silicate concentration
    - `Salinity`: Water salinity
    - Temperature metrics (`mean_temp`, `max_temp`, `min_temp`)
  - Principal components:
    - `PC1`, `PC2`, `PC3`: Environmental principal components
    - `PC1_1yrlag`, `PC2_1yrlag`, `PC3_1yrlag`: 1-year lagged PCs
    - `PC1_2yrlag`, `PC2_2yrlag`, `PC3_2yrlag`: 2-year lagged PCs



## Scripts

### Data Processing and Analysis Scripts

1. `0_NMRS-PCA.R`
   - Performs Principal Component Analysis on environmental time-series data
   - Creates environmental covariates by combining NMRS and SST data
   - Generates PCA visualization plots
   - Outputs: processed environmental PC time-series with lags

2. `1_inla-forecasts.R`
   - Implements forecasting using INLA (Integrated Nested Laplace Approximation) models
   - Compares two approaches:
     - Legacy split: fixed training period
     - Modern split: rolling training period
   - Handles multiple species and model types (ar1, rw1, rw2)
   - Outputs: model forecasts saved in RDA format

3. `2_process-forecasts.R`
   - Processes and visualizes forecast results
   - Calculates Mean Absolute Scaled Error (MASE)
   - Creates comparison plots across different:
     - Time horizons
     - Species types (Increasing, Collapsing, Fluctuating, Stable)
     - Model types
   - Outputs: forecast accuracy visualization figures

4. `3_inla-forecasts-simtesting`
   - Sets up parameters for simulation testing 
   - Compares two approaches:
     - Legacy split: fixed training period
     - Modern split: rolling training period
   - Handles multiple parameter settings and model types (ar1, rw1, rw2)
   - Outputs: model forecasts saved in RDA format

4. `4_sim-testing-plots.R`
   - Creates visualization plots for simulation testing results
   - Combines results from AR1 and RW1/RW2 simulations
   - Generates comparative plots showing:
     - Forecast accuracy across different scenarios (25%, 50%, 75%, 100% K)
     - Performance of Legacy vs Modern split approaches
     - Error metrics across different forecast horizons
   - Outputs: Combined simulation testing plots and simplified forecast horizon plot

5. `5_parameter-plots.R`
   - Generates parameter estimation and model diagnostic plots
   - Analyzes hyperparameters (Rho, Precision) across different:
     - Species
     - Time periods
     - Split methods (Legacy vs Modern)
   - Creates visualizations for:
     - Parameter estimates with uncertainty ranges
     - Fixed effects comparisons
     - Model comparison plots
   - Outputs: Parameter plots and diagnostic visualizations

6. `6_inla-forecasts-simtestng-example-plots.R`
   - Generate example plots of simulation results 

7. `7_species-dynamics-categories.R`
   - Calculation of statistics characterizing time-series for case-study species 




### Supporting Files
- `functions-sim-predictions.R`: Helper functions for simulation predictions
- `inla-forecast-functions.R`: Core functions for INLA model fitting and forecasting

