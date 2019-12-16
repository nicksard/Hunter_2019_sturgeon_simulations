#originally created on: Jan 14, 2019
#by: Nick Sard

#ABOUT: This script was written to analyze information from confusion matrices 
# to understand how full and half sib dyad inferences are correct/incorrect, and if incorrect, what type of 
#mistake was made - e.g. if a full-sib dyad was incorrectly identified, was it called a half-sib or unrelated dyad instead?
#and outputs (compares estiamtes for a number of paramters to real values) from power simulations 
#are used to determine how confusion correlations with an increase in true parents, as well as 
#the magnitude of the bias

#loading in libraries
library(tidyverse)

#loading my own functions
# - none - 

#setting working directory and reading in data
setwd("C:/Users/sard/Google Drive/R/Projects/Hunter_2019_sturegon_simulations/Step 3 - Post-processing")

#For each simulation, all information associated with their respective Best Configuration files are post.processed on a Linux individually
#I took all of that information and created a single file that contains all relavent information for each simulation

#reading in combined output file
df <- read.table(file = "Input/all.simulation.output.long.form.combined.txt",header = T,sep = "\t",stringsAsFactors = F)
head(df, n = 20)

################
### graphics ###
################

#confusion matrix of dyads
tiff(filename = "Output/confusion.matrix.50.figure.color.tiff",width = 11,height = 8,units = "in",res = 300)
ggplot(df %>% filter(noff == "50"), aes(x=real,y=val,color=type))+
  facet_grid(inferred~real.type)+
  geom_point(size=2,alpha=.5)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()+
  labs(x="Number of known parents",y="Proportion of dyads",color="Simulation type")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16,color="black"),
        strip.text = element_text(size=14))
dev.off()

#confusion matrix of dyads
tiff(filename = "Output/confusion.matrix.125.figure.color.tiff",width = 11,height = 8,units = "in",res = 300)
ggplot(df %>% filter(noff == "125"), aes(x=real,y=val,color=type))+
  facet_grid(inferred~real.type)+
  geom_point(size=2,alpha=.5)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()+
  labs(x="Number of known parents",y="Proportion of dyads",color="Simulation type")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16,color="black"),
        strip.text = element_text(size=14))
dev.off()

#confusion matrix of dyads
tiff(filename = "Output/confusion.matrix.750.figure.color.tiff",width = 11,height = 8,units = "in",res = 300)
ggplot(df %>% filter(noff == "750"), aes(x=real,y=val,color=type))+
  facet_grid(inferred~real.type)+
  geom_point(size=2,alpha=.5)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()+
  labs(x="Number of known parents",y="Proportion of dyads",color="Simulation type")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16,color="black"),
        strip.text = element_text(size=14))
dev.off()

head(df)
#the relationship between the real number of parents and the estimate based on the pedigree, as the real number increases
tiff(filename = "Output/bias.figure.color.tiff",width = 11,height = 8,units = "in",res = 300)
ggplot(df, aes(x=real,y=ratio,color= type))+
  facet_wrap(~noff2,scales = "free_x")+
  geom_point(size=2,alpha=.5)+
  scale_color_brewer(palette = "Set1")+
#  geom_smooth(method = "lm",se = F)+
  theme_bw()+
  labs(x="Number of known parents",y="Bias in number of parents estimation",color="Simulation type")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16,color="black"),
        strip.text = element_text(size=16))
dev.off()

#fin!