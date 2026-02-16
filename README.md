# Predicting_e026_coexistence
All R scripts and datasets are provided for replicating all main text and supplementary figures for Catford et al. (2026), Multiple mechanisms are required to predict grass community composition.

Please cite as: Catford JA, Graham LJ, Shepherd HER, Hauser CE, Munro NT et al., (2026) Multiple mechanisms are required to predict grass community composition. Ecology letters. doi: TBC

Included in this repository are R scripts and data required for analyzing mechanistic predictions of community coexistence in e026 Cedar Creek. The code uses hierarchical multi-mechanism niche model of plant community assembly provided in: https://github.com/laurajanegraham/simulateCoexistence to predict the biomass of up to 5 grasses sown together across a soil N gradient. For full methods, see the above citation.

## Repository Structure

This repository contains two main folders:

- `R_scripts/`: Contains all R code required to produce modelled results and analyses included in the main manuscript.
  - `test_biomass_predictions.R`: Simulates all possible combinations of attributes in the model and outputs model estimates.
  - `output_results.R`: Takes model simulations produced in `test_biomass_predictions.R` and formats the data for subsequent analysis.
  - `Producing_Figure2.R`: Recreates Figure 2 in the main manuscript.
  - `Producing_Figure3.R`: Recreates Figure 3 in the main manuscript.
  - `Producing_Figure4.R`: Recreates Figure 4 in the main manuscript.
  - `Producing_supplementary_figures.R`: Recreates all supplementary figures in the manuscript.

- `Data/`: Contains all datasheets required for the analysis.
  - `non_zero_replicates_[0-11]_year6.R`: Simulation outputs produced in `output_results.R`. The number in each filename refers to the number of attributes switched on in the model simulations.
  - `files_combined.R`: Datasheet of simulated file outputs from `test_biomass_predictions.R`, including the attributes switched on in each simulation.
  - `switches.RDA`: List of all switch combinations and the number of attributes switched on in each combination.

