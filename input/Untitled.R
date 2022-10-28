## calibration w multiple constraints


library(plyr)
library(tidyverse)
library(ape)
library(phytools)
library(jtools)
library(scales)
library(treeio)
library(ggtree)

## time constrains

neoaves_min <- 66.0
neoaves_max <- 99.6 
neoaves_outgroup <- c('StrCam')

#Read Trees

CDS_chronos <- read.nexus("../RAxML_br/CDS_tree")

root(CDS_chronos,outgroup = neoaves_outgroup)

ggplot(CDS_chronos, aes(x, y)) + geom_tree() + theme_tree()

## calibration

CDS_cal <- makeChronosCalib(CDS_chronos, interactive = TRUE)


#test_cal <- makeChronosCalib(CDS_chronos,age.min = neoaves_min,age.max = neoaves_max, node = 34)

CDS_edge_list <- compute.brlen(CDS_chronos)$edge.length
CDS_branch_list <- branching.times(compute.brlen(CDS_chronos))


#loop 

for(i in 1:10){
  CDS_chronos_iter <- chronos(CDS_chronos,calibration = CDS_cal)
  CDS_edge_list <- rbind(CDS_edge_list, CDS_chronos_iter$edge.length)
  CDS_branch_list <- rbind(CDS_branch_list, branching.times(CDS_chronos_iter))
}
