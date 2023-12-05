# Situating nature loss and conservation in global mental health

This repository is associated with the manuscript 'Situating nature loss and conservation in global mental health'. Specifically, it presents a re-analysis of the data published in 'Predicting the impacts of land management for sustainable development on depression risk in a Ugandan case study' (Pienkowski et al., 2022, https://doi.org/10.1038/s41598-022-14976-3). 

The repository contains the following files: 

- 'Anon_case_study.Rproj': An R Project.
- 'HH_locations': A folder containing rds file 'DF_anly.rds'. This rds file contains anonymised and spatially dislocated data collected by Pienkowski et al. (2022), used in the re-analysis. 
- 'functions.R': A script containing useful functions. 
- 'Processing.R': A script that processes the data contained in DF_anly.rds, including imputing missing data (0.2% of values) and extracting plausible values of latent variables, saved in ten datasets. 
- 'Analysis.R': A script that performs the statistical analysis. To run the entire script, users must download WDPA shapefile data (https://www.protectedplanet.net/en/thematic-areas/wdpa?tab=WDPA) and GRID3 settlement extends data (version 1.1., https://data.grid3.org/datasets/GRID3::grid3-uga-settlement-extents-v1-1/explore) for Uganda. These data need to be saved in the 'Spatial_data' folder and renamed following the names denoted in the script. 