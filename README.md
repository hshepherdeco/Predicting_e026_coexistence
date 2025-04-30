# Predicting_e026_coexistence
All R scripts and datasets are provided for replicating all main text and supplementary figures for Catford et al., Mechanistic flexibility when predicting grassland community composition DOI: ...

## Repository Structure

This repository contains two main folders:

- `R_scripts/`: Contains all R code required to produce modelled results and analyses included in the main manuscript:
  - `xxxxx.R`: Simulates all possible combinations of attributes in the model and outputs model estimates.
  - `output_results.R`: Takes model simulations produced in `xxxxx.R` and formats the data for subsequent analysis.
  - `Producing_Figure2.R`: Recreates Figure 2 in the main manuscript.
  - `Producing_Figure3.R`: Recreates Figure 3 in the main manuscript.
  - `Producing_Figure4.R`: Recreates Figure 4 in the main manuscript.
  - `Producing_supplementary_figures.R`: Recreates all supplementary figures in the manuscript.

- `Data/`: Contains all datasheets required for the analysis:
  - `non_zero_replicates_[0-11]_year6.R`: Simulation outputs produced in `output_results.R`. The number in each filename refers to the number of attributes switched on in the model simulations.
  - `files_combined.R`: Datasheet of simulated file outputs from `xxxxx.R`, including the attributes switched on in each simulation.
  - `switches.RDA`: List of all switch combinations and the number of attributes switched on in each combination.

