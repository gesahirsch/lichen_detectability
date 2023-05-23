This folder contains all the data necessary to replicate the results and figures published in the manuscript and supplementary materials.
There are the following files:


#  experience_long.csv  #
#-----------------------#
table indicating (0/1) whether the observer who conducted the respective survey had previous experience with the respective species, with 2 columns:
- CLNR (unique identifier that can link the detections to the respective sites; the content of this column matches that of the column sites$CLNR)
- species (name of the species)
- visit (factor with 2 levels specifying the first='A' or second='K' visit of the survey)


#  sites.csv  #
#-------------#
table of 826 sites x 14 columns describing the sites:
- CLNR (unique identifier that can link the sites to the underlying coordinates and lichen observations in the SwissLichens database)
- altitude (in meter above sea level)
- temperature (in degrees centigrade; mean annual average temperature)
- precipitation (in mm; mean annual amount of precipitation)
- solar (in percent; mean relative sunshine duration)
- observerA (unique identifier for the person who conducted the first survey 'A-Aufnahme'; NA where no field visit was made)
- observerK (unique identifier for the person who conducted the second survey 'Kontrollaufnahme'; NA where no survey was made)
- dateA (date when the first survey took place; NA where no field visit was made)
- dateK (date when the second survey took place; NA where no second survey was made)
- K (indicator: 1=second survey was conducted (called 'Kontrollaufnahme' in German in the original study); 0=no second survey was conducted)
- feldbesuch (indicator: 1=field visit was made; 0=no field visit was made, either due to certainty about the site's lack of epiphytic substrate or due to financial restrictions)
- epiphytic.substrate (indicator: 1=some substrate for epiphytic lichens is available; 0=no substrate for epiphytic lichens available)
- biogeo (factor with 5 levels specifying the biogeographic region this site falls into: Jura, Plateau=Mittelland, Prealps=Voralpen, Alps=Alpen, Southern Alps=Südalpen)
- forest (factor with 2 levels specifying whether a site was considered forest=Wald or non-forest=Nichtwald)
Note: we provide the data from all 826 sites contained in the original stratified random sample. The analysis uses only the data from the 416 sites that contained any potential substrate for corticolous lichens. This selection happens when running the file 1_code/0_bundle_data.R.


#  species_detections.csv  #
#--------------------------#
table with 3 columns where each row stands for a species detection during a visit at a site:
- CLNR (unique identifier that can link the detections to the respective sites; the content of this column matches that of the column sites$CLNR)
- species (name of the taxonomic unit/species)
- visit (factor with 2 levels specifying whether the detection occurred during the first='A' or second='K' visit)


#  species_indices.csv  #
#-----------------------#
table with a indices for all lichen species with 2 columns:
- species (names of the 373 lichen species assessed)
- conspicuousness (indicator: 1=species is considered conspicuous; 0=species is considered inconspicuous)
- identifiability (pseudo-continuous variable: 1=species can be identified with a handlens in the field; 0.5=species can be identified in the field with spot tests; -0.5=species requires a dissecting microscope for identification; -1=species requires microscopy or thin layer chromatography for identification)


#  substrate_long.csv  #
#----------------------#
table indicating (0/1) whether any suitable substrate is available for a given species at a given site, with 3 columns:
- CLNR (unique identifier that can link the detections to the respective sites; the content of this column matches that of the column sites$CLNR)
- species (name of the species)


#  experience_gain.csv  #
#-----------------------#
array (373 x 7 x 2) indicating (0/1) whether a given observer [ ,1:7, ] had previously had experience with a given species [1:373, , ] at the beginning [ , ,1] vs. the end [ , ,2] of the study period. This file is not required for the model fitting, but it is used to generate some of the numbers in the manuscript by running the file 1_code/2_extract_results_and_make_figures.R.



## ADDITIONAL FILE THAT WILL BE GENERATED INTO THIS FOLDER
#----------------------------------------------------------#

Running the file 1_code/0_bundle_data.R will generate one more file and place it into this folder:

#  data.RData  #
#--------------#
bundled list of data in the correct format for JAGS to use it as input for the occupancy model. It is generated from the above data files but includes only data for the 416 sites that had epiphytic substrate.