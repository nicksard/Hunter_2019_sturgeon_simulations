#originally created on: August 21, 2017
#by: Nick Sard

#ABOUT: This script was written to create a input file (Input3.Par) for the simulation executable within COLONY using sturgeon-specific data

#loading in libraries
library(tidyverse)

#setting working directory and reading in data
setwd("C:/Users/sard/Google Drive/R/Projects/Hunter_2019_sturegon_simulations/Step 1 - Simulating data for colony/")

#to stop scientific notation
options(scipen=999)

#source scripts
source("Source/brd.mat.functions_10_08_2018.R")
source("Source/colony.poly.sim.create.R")

#using the actual number of alleles observed for sturgeon in the St. Clair-Detriot River system
alleles.per.locus <- c(13,16,13,8,5,5,7,10,11,16,11,14,4,23,8,15,14,12)
my.n.loci <- length(alleles.per.locus)

#########################
### setting variables ###
#########################

#Maximum expected number of offspring collected using specific gear
maxoff <- 50
minoff <- 50

#Minimum and maximum expected sex ratio - Based on data from Black Lake
min.sex.ratio <- 3
max.sex.ratio <- 3

#Minimum and maximum expected number of parents spawning in study area
min.parents <- 10
max.parents <- 500

#Minimum and maximum expected number of fertilized eggs per female
min.fert = 50000
max.fert = 700000

#setting minimum and maximum bounds for mean number of mate pairs
min.mates <- 3
max.mates <- 3

#creating object to store information about each simulation and parameter estimation
sim.info <- NULL

#now for the simulations
k <- 1
k <- NULL
for(k in 1:100){
  
  ############################
  ### setting output names ###
  ############################
  
  mating.str.name <- paste0("Output/sims/breeding_matrices/mating_matrix_sim_",k,".txt")
  sim.input.name <- paste0("Output/sims/input_par/input3_sim_",k,".Par") 
  colony.input.name <- paste0("Output/sims/colony_dat/colony2_sim_",k,'.Dat')
  
  ####################################################
  ### Clearing tmp varaibles before each simuation ###
  ####################################################
  n.mom <- NULL
  n.dad <- NULL
  n.par <- NULL
  n.off <- NULL
  sex.ratio <- NULL
  sex.ratio2 <- NULL
  picks <- NULL
  mat.str <- NULL
  ms1 <- NULL
  mat.sub <- NULL
  mp.lost <- NULL
  ped2 <-  NULL
  mat2 <- NULL
  ms2 <- NULL
  
  ######################################################
  ### Create a breeding matrix for entire study area ###
  ######################################################
  
  print(k)
  
  #defining the number of moms, dads, and rs distribution (poisson in this case) (20 and 5 mean, sd)
  
  #picking the number of parents from a uniform distrbution
  n.par <- round(runif(n = 1,min = min.parents, max = max.parents))
  n.par
  
  # assuming the the sex ratio of successful males to females is NOT 1:1
  sex.ratio <- round(runif(n = 1,min = min.sex.ratio,max = max.sex.ratio),1)
  sex.ratio
  
  #needed a way to pick of combination of males and females that would exactly (or as close as possible) to
  #the sex ratio that was choosen, here is my solution
  
  #enumerate all possible combinations of males and females, calculate sex ratio (sr) between them,
  #get the absoluate difference between their sr and the real sr (drawn above), and pick the one with smallest difference
  picks <- expand.grid(1:n.par,1:n.par) %>%
    mutate(sr =  Var1/Var2) %>%
    mutate(npar = Var1 + Var2) %>%
    mutate(diffs = abs(sr - sex.ratio)) %>%
    mutate(diffnp = abs(npar - n.par)) %>%
    filter(diffnp == min(diffnp)) %>%
    filter(diffs == min(diffs))
  picks
  
  #sometimes its possible to pick more than one pick - e.g. when the sex ratio is 2, so I put in an if statement
  #to randomly select one of the rows in that case
  if(nrow(picks)== 1){
    n.mom <- picks$Var2
    n.dad <- picks$Var1
    sex.ratio2 <- round(picks$sr,2)
    diff <- picks$diff
  } else {
    warning(paste("There are",nrow(picks),"sex ratios that work. Picking one randomly."))
    my.row <- as.numeric(sample(x = row.names(picks),size = 1))
    picks <- picks[my.row,]
    n.mom <- picks$Var2
    n.dad <- picks$Var1
    sex.ratio2 <- round(picks$sr,2)
    diff <- picks$diff
  }
  
  #using breding matrix function to create the breeding matrix
  mat.str <- brd.mat(moms = n.mom,dads = n.dad,min.mates = min.mates,max.mates = max.mates,min.fert = min.fert,max.fert = max.fert)
  head(mat.str)

  #getting summary stats
  ms1 <- mat.stats(mat = mat.str)
  ms1$type <- "Before"
  ms1$lost <- 0
  ms1$sim <- k
  ms1
  
  #######################################
  ### Subsamping full breeding matrix ###
  #######################################
  
  #randomly selecting the number of offspring that will be sampled in a given sampling effort
  n.off = round(runif(n = 1,min = minoff, max = maxoff))
  n.off
  
  #creating full pedigree with that breeding matrix, and sub-sampling each mate pair randomly
  mat.sub <- mat.sub.sample(mat = mat.str,noff = n.off)
  head(mat.sub)
  mat.sub
  
  #recording how many are lost
  mplost <- nrow(mat.sub[mat.sub$off1 == 0,])
  mplost
  
  #getting rid of the zeros
  mat.sub <- mat.sub[mat.sub$off1 != 0,]
  
  #converting to small "complete" genetic pedigree
  ped2 <- mat2ped(ped = mat.sub)
  head(ped2)
  ped2
  
  #converting back to a breeding matrix to calculate a range of stats
  mat2 <- ped2mat(ped = ped2)
  head(mat2)
  
  #getting numbers of mates and rs per sex
  ms2 <- mat.stats(mat = mat2)
  ms2$type <- "After"
  ms2$lost <- mplost
  ms2$sim <- k
  ms2
  
  #saving information associated with each breeding matrix for graphics
  sim.info1 <- rbind(ms1,ms2)
  sim.info <- rbind(sim.info,sim.info1)
  sim.info
  mat2
  
  ##########################
  ### Making Colony file ###
  ##########################
  
  #making the file
  colony.sim.create(mat2, update.alfs = 0, spp.type = 2, selfing.rate = 0, inbreeding = 0,
                    alf.dist = NA, ploidy = 8, fem.gamy = 0, mal.gamy = 0,
                    sib.prior = 0, known.alfs = 0, n.reps = 1, run.reps = 1, run.length = 3,
                    monitor = 0, windows.version = 0, likelihood.type = 1, likelihood.precision = 3,
                    prob.mom = 0.01, prob.dad = 0.01, n.mating.str = 1, nmarkers = my.n.loci, gt.err = 0,
                    mut.err = 0.01, pr.miss.gt = 0.001, marker.type = 0, n.alleles = alleles.per.locus, allele.dist = 2,
                    map.length = -1, inbrd.coef = 0, pat.sib.size = 1, mat.sib.size = 1,n.mom = NA, n.dad = NA,
                    convert.codom = 1, sib.scale = 0, codom_dropout =0.01, codom_false.allele =0.01,
                    delete.sup.files = T)
  
  #writing mating matrix to file and modifying some names
  write.table(x = mat2,file = mating.str.name,append = F,quote = F,sep = " ",row.names = F,col.names = F)
  file.rename(from = "input3.Par",to = sim.input.name)
  file.rename(from = "COLONY2_1.DAT",to = colony.input.name)
} #end of creating simulated datasets and colony input files


#writing each to file for evulation via graphics
write.table(x = sim.info,file = "Output/sims/simulation.information.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

#fin!