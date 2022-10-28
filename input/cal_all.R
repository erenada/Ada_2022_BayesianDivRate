### Calibrating all trees with Neognathae + Galliformes CALIBRATION

library(tidyverse)
library(ape)
library(phytools)
library(treeio)
library(ggtree)

neoaves_min <- 66.0
neoaves_max <- 99.6 
neoaves_outgroup <- c('StrCam')

#read trees


ALL50K_chronos <- root(read.tree("../RAxML_br/HKY85/RAxML_result.ALL50K"), outgroup = "StrCam", resolve.root = T)

gammaStat(ALL50K_chronos)

