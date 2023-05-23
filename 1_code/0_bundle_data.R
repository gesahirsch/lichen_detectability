library(conflicted)
library(tidyverse) # reduce() & all dplyr %>% operations
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")


#  load data  ####

sites <- read.csv("0_data/sites.csv", sep=";")
species_detections <- read.csv('0_data/species_detections.csv', sep=';')
sub <- read.csv("0_data/substrate_long.csv", sep=";")
exper <- read.csv("0_data/experience_long.csv", sep=";")
species_indices <- read.csv("0_data/species_indices.csv", sep=";")



#  basics  ####

levspecies <- levels(as.factor(species_detections$species))
nspecies <- length(levspecies)
levsites <- levels(as.factor(sites$CLNR))
nsites <- length(levsites)
nvisits <- 2
nobservers <- 7



#  visits  ####

#  = indicator matrix to specify over how many visits the loop should run
visits <- matrix(data=NA, nrow=nsites, ncol=2)
visits[,1] <- 1
visits[,2] <- ifelse(sites$K==1, 2, 1)



#  experience  ####

exper$experience <- 1 # all rows represent a species encountered previously
head(exper)

# prepare 3D array
experience <- array(data=NA, dim=c(nsites, nvisits, nspecies))
dimnames(experience) <- list(levsites,
                                  c('A','K'),
                                  levspecies)
dim(experience)
experience[1:10,,1:4]

# ensure that order and levels of factors are kept
exper$CLNR <- factor(exper$CLNR, levels=levsites)
exper$species <- factor(exper$species, levels=levspecies)

experience[,1,] <- exper %>% 
  filter(visit=='A') %>% 
  spread(key=species, value='experience', drop=FALSE, fill=0) %>% 
  select(-c(CLNR,visit)) %>% 
  as.matrix()
experience[,2,] <- exper %>% 
  filter(visit=='K') %>% 
  spread(key=species, value='experience', drop=FALSE, fill=0) %>% 
  select(-c(CLNR,visit)) %>% 
  as.matrix()
experience[which(sites$K==0),2,] <- NA # no second survey was conducted
experience[1:10,,1:3]



#  dethist  ####

species_detections$detection <- 1 # all rows represent detected species
species_detections <- unique(species_detections)

# prepare 3D array
dethist <- array(data=NA, dim=c(nsites, nvisits, nspecies))
dimnames(dethist) <- list(levsites,
                          c('A','K'),
                          levspecies)
dim(dethist)
dethist[1:10,,1:4]

# ensure that order and levels of factors are kept
species_detections$CLNR <- factor(species_detections$CLNR, levels=levsites)
species_detections$species <- factor(species_detections$species, levels=levspecies)

dethist[,1,] <- species_detections %>% 
  filter(visit=='A') %>% 
  spread(key=species, value='detection', drop=FALSE, fill=0) %>% 
  select(-c(CLNR,visit)) %>% 
  as.matrix()
dethist[,2,] <- species_detections %>% 
  filter(visit=='K') %>% 
  spread(key=species, value='detection', drop=FALSE, fill=0) %>% 
  select(-c(CLNR,visit)) %>% 
  as.matrix()
dethist[which(sites$K==0),2,] <- NA # no second survey was conducted
dethist[1:10,,1:3]



#  substrate  ####

sub$substrate <- 1 # all rows represent suitable substrates

# ensure that order and levels of factors are kept
sub$CLNR <- factor(sub$CLNR, levels=levsites)
sub$species <- factor(sub$species, levels=levspecies)

substrate <- sub %>% 
  spread(key=species, value='substrate', drop=FALSE, fill=0) %>% 
  select(-CLNR) %>% 
  as.matrix()
substrate[which(sites$feldbesuch==0),] <- NA # no survey was conducted (= presence/absence of substrate not known)
dimnames(substrate) <- list(levsites,
                            levspecies)
dim(substrate)
substrate[1:10,1:2]



#  standardize continuous variables  ####

# check for missing values
apply(sites, 2, function(x) sum(is.na(x)))

standardized <- sites %>% 
  mutate_at(c('altitude','temperature','precipitation','solar'), ~(scale(.) %>% as.vector)) %>% 
  select(altitude, temperature, precipitation, solar)
head(standardized)
apply(standardized, 2, mean); apply(standardized, 2, sd)



#  selection of sites to include  ####

index.sites <- which(sites$epiphytic_substrate==1); length(index.sites)
nsites <- length(index.sites)
levsites <- levsites[index.sites]

# for goodness-of-fit test
help <- sites[index.sites,]
index <- which(help$K==1); length(index)



#  bundle and export data  ####

data <- NULL
str(data <- list(nsites=length(index.sites),
                 nvisits=nvisits,
                 nspecies=nspecies,
                 visits=visits[index.sites,],
                 index=index,    # only needed for goodness-of-fit assessment
                 
                 y=dethist[index.sites,,],
                 
                 altitude.st=standardized$altitude[index.sites],
                 altitude.st.2=(standardized$altitude[index.sites])^2,
                 precipitation.st=standardized$precipitation[index.sites],
                 precipitation.st.2=(standardized$precipitation[index.sites])^2,
                 solar.st=standardized$solar[index.sites],
                 solar.st.2=(standardized$solar[index.sites])^2,
                 substrate=substrate[index.sites,],
                 
                 observer=as.matrix(sites[index.sites,c('observerA','observerK')]),
                 nobservers=nobservers,
                 conspicuousness=species_indices$conspicuousness,
                 identifiability=species_indices$identifiability,
                 experience=experience[index.sites,,]
                 ))

rm(sites, standardized, sub, substrate, visits, species_detections,
   exper, experience, dethist, help, species_indices)

save.image("0_data/data.RData", version=2)
