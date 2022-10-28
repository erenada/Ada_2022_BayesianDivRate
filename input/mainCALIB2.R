## unannotated ROOT + Neognathae CALIBRATION


library(tidyverse)
library(ape)
library(phytools)
library(treeio)
library(ggtree)

## time constrains

neoaves_min <- 66.0
neoaves_max <- 99.6 
neoaves_outgroup <- c('StrCam')

#read tree
unannotated_chronos <- root(read.tree("../RAxML_br/GTRG/RAxML_result.unannotated"), outgroup = "StrCam", resolve.root = F)

plot(unannotated_chronos)

#ggplot(unannotated_chronos, aes(x, y)) + geom_tree() + theme_tree()

#fix the root

unannotated_chronos$edge.length[which(!(unannotated_chronos$edge[,1] %in% unannotated_chronos$edge[,2]))] <- sum(unannotated_chronos$edge.length[which(!(unannotated_chronos$edge[,1] %in% unannotated_chronos$edge[,2]))])/2

#plot.phylo(unannotated_chronos, use.edge.length = F)

nodelabels()

#ggtree(unannotated_chronos, ladderize=T, branch.length="none")

#calibration start


unannotated_cal <- unannotated_cal <- add_row(makeChronosCalib(unannotated_chronos,age.min = neoaves_min, age.max = neoaves_max, soft.bounds=FALSE), node = 35, age.min = 66, age.max=70,soft.bounds = FALSE)


#get edge and branch lists

unannotated_edge_list <- compute.brlen(unannotated_chronos)$edge.length
unannotated_branch_list <- branching.times(compute.brlen(unannotated_chronos))

#repeat

for(i in 1:100){
  unannotated_chronos_iter <- chronos(unannotated_chronos,calibration = unannotated_cal, model = "clock")
  unannotated_edge_list <- rbind(unannotated_edge_list, unannotated_chronos_iter$edge.length)
  unannotated_branch_list <- rbind(unannotated_branch_list, branching.times(unannotated_chronos_iter))
}

plot(unannotated_chronos_iter)

#choosing median ages


unannotated_edge_list <- unannotated_edge_list[-1,]
unannotated_branch_list <- unannotated_branch_list[-1,]

unannotated_branches <- colnames(unannotated_branch_list)

unannotated_medians <- c()
for(i in 1:ncol(unannotated_branch_list)){
  unannotated_medians[colnames(unannotated_branch_list)[i]] <- median(unannotated_branch_list[,i])
}

# 
unannotated_timetree <- compute.brlen(unannotated_chronos)

is.ultrametric(unannotated_timetree)
unannotated_median_edge <- c()
for(i in 1:ncol(unannotated_edge_list)){
  unannotated_median_edge<-c(unannotated_median_edge,median(unannotated_edge_list[,i]))
}
unannotated_timetree$edge.length <- unannotated_median_edge

plot(unannotated_timetree)

write.tree(unannotated_timetree, file = "../TimeTrees/unannotated_CLCK_medianTT.tree")
