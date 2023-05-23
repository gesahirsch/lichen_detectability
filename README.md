# LICHEN_DETECTABILITY
---------------------------------------------------------------------------------------------------

This repository contains the cleaned data and all the code used for the data analysis and the generation of the figures and numbers presented in the article von Hirschheydt, G., Kéry, M., Ekman, S., Stofer, S., Dietrich, M., Keller, C., Scheidegger, C. Occupancy models reveal limited detectability of lichens in a standardised large-scale monitoring.

All data and code currently underlies a strict copyright.


#  REPRODUCIBILITY OF RESULTS  ##
#-------------------------------#

Our results and figures should be reproducible (with small differences in the latter digits due to stochasticity of the MCMC sampler) with the provided code.
Below, we describe the structure and content of this repository.
The file "workflow.Rmd" guides the user through the R-files in the correct order so that you can:
- bundle the cleaned data into a data list readable for JAGS
- fit the multi-species occupancy model to the data and store the output
- assess the goodness-of-fit of the model to the data
- conduct a prior sensitivity analysis with 2 additional sets of priors
- extract the summary statistics reported in the manuscript and supplementary materials
- generate the figures shown in the manuscript and supplementary materials
IMPORTANT: The fitting of the occupancy models takes a long time: Each model file runs for approximately 1.5 days. We would therefore recommend to parallelize the model fitting process, e.g., by using a high-performance cluster (necessary singularity container is provided as .def and .sif file).
If you do not want to take the time to fit the occupancy models, but you would nevertheless like to explore the output of ours, please contact Gesa von Hirschheydt so the data package (ca. 1.35 GB) can be sent to you.


#  STRUCTURE OF THIS REPOSITORY  ##
#---------------------------------#

# FILE: lichen_detectability.Rproj
This file stores the information about the R-project (RStudio) of this repository. You can open the project by clicking on this file. Opening the project in this way will automatically define this repository as the working directory of the R-session.
If you are not using RStudio or you do not want to work with R-projects, you must make sure that you set this repository as the working directory before you run any of the R-code.

# FILE: workflow.Rmd
This file explains in which order the provided R-files must be executed in order to reproduce the results we presented in the manuscript.

# FILE: workflow.html
Visually more appealing version to explain in which order the provided R-files must be executed in order to reproduce the results we presented in the manuscript.

# FOLDER: 0_data
This folder contains the cleaned data used to generate the results. "Cleaned" means:
- coordinates have been removed to protect the local populations of endangered and protected lichen species
- observer names have been replaced with identification numbers ('pseudonymization') to protect their person
- some errors in the original data have been corrected (i.e., mis-identifications, erroneous survey dates or coordinates)

# FOLDER: 1_code
This folder contains 1 R-file to bundle the cleaned data into a list that can be used as model input in JAGS, 3 R-files to fit the multi-species occupancy models to the data (inkl. goodness-of-fit test and prior sensitivity analysis), and 1 R-file to generate the figures and extract the summary statistics reported in the manuscript or supplementary materials.
Note: Each R-file of the format 1***.R takes approximately 1.5 days to run.

# FOLDER: 2_output
This folder is currently empty (except a README file), but all model output generated when running the R-files in folder 1_code/ will be saved to this folder, i.e., model output, figures, and summary tables.

