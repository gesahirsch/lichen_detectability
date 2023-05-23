This folder contains all the code necessary to replicate the results and figures published in the manuscript and supplementary materials.
There are the following files:


#  0_bundle_data  #
#-----------------#
Combines the data from the folder 0_data/ into a list that can be submitted to JAGS for the occupancy model. It also reduces the data from the original 826 sites to only 416 by excluding all sites that lack epiphytic substrate.
Running this script will generate the file 0_data/data.RData.


#  1a_run_jags_model.R  #
#-----------------------#
Includes the JAGS-formulation of the multi-species occupancy model, initial values, MCMC settings, parameters to be traced, and saving of the output.
Running this script will generate the file 1_code/jags_model.txt and the file 2_output/output_jags_model.RData.


#  1b_run_jags_model_gof.R  #
#---------------------------#
Goodness-of-fit assessment of the multi-species occupancy model: The likelihood in this model is identical to the one in 1a_run_jags_model.R, but here, other parameters are traced to generate Bayesian p-values to assess model fit.
Includes the JAGS-formulation, initial values, MCMC settings, parameters to be traced, and saving of the output.
Running this script will generate the file 1_code/jags_model_gof.txt and the file 2_output/output_jags_model_gof.RData.
 

#  1c_run_jags_model_alternative_priors1.R  #
#-------------------------------------------#
For prior sensitivity analysis: The likelihood in this model is identical to the one in 1a_run_jags_model.R, but here, the prior for the regression coefficients is different from the one used in 1_run_jags_model.R.
Includes the JAGS-formulation, initial values, MCMC settings, parameters to be traced, and saving of the output.
Running this script will generate the file 1_code/jags_model_priors1.txt and the file 2_output/output_jags_model_priors1.RData.
 

#  1d_run_jags_model_alternative_priors2.R  #
#-------------------------------------------#
For prior sensitivity analysis: The likelihood in this model is identical to the one in 1a_run_jags_model.R, but here, the prior for the standard deviation of the random effects is different from the one used in 1_run_jags_model.R.
Includes the JAGS-formulation, initial values, MCMC settings, parameters to be traced, and saving of the output.
Running this script will generate the file 1_code/jags_model_priors2.txt and the file 2_output/output_jags_model_priors2.RData.
 

#  2_extract_results_and_make_figures.R  #
#----------------------------------------#
Generates the table, figures and summary statistics/numbers from the manuscript and supplementary materials.
Running this script will save 8 figures in .tif format and 2 tables in .csv format and place them in the folder 2_output/.