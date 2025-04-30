# Predicting_e026_coexistence
All R scripts and datasets are provided for replicating all main text and supplementary figures for Catford et al.,... doi: ...

This repository contains two folders: /R_scripts and /Data

R_scripts: Contains all R code required to produce modelled results and analysis included in the main manuscript
            1. xxxxx.R: Simulate all possible combinations of attributes in the model and output model estimates.
            2. output_results.R: Take model simulations produced in xxxx.R and format the data for subsequent analysis. 
            3. Producing_Figure2.R: Recreate figure 2 in the main manuscript text.
            4. Producing_Figure3.R: Recreate figure 3 in the main manuscript text.
            5. Producing_Figure4.R: Recreate figure 4 in the main manuscript text.
            6. Producing_supplementary_figures.R: Recreate all supplementary figures in the main manuscript text.

Data: Contains all datasheets required for the analysis 
            1. non_zero_replicates_[0-11]_year6.R: Simulation outputs produced in output_results.R. Variable number in each file refers to the number of attributes switched               on in the model simulations
            2. files_combined.R: datasheet of simulated file ouputs from xxxx.R and the attributes switched on in each simulation
            3. switches.RDA: list of all switch combinations and the number of attributes switched on in each combination


