## pseudogene ROOT + Neognathae CALIBRATION


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
UTR_chronos <- root(read.tree("../RAxML_br/RAxML_result.UTR"), outgroup = "StrCam", resolve.root = F)

plot(UTR_chronos)

#ggplot(UTR_chronos, aes(x, y)) + geom_tree() + theme_tree()

#fix the root

#UTR_chronos$edge.length[which(!(UTR_chronos$edge[,1] %in% UTR_chronos$edge[,2]))] <- sum(UTR_chronos$edge.length[which(!(UTR_chronos$edge[,1] %in% UTR_chronos$edge[,2]))])/2

#plot.phylo(UTR_chronos, use.edge.length = F)

nodelabels()

#ggtree(UTR_chronos, ladderize=T, branch.length="none")

#calibration start


UTR_cal <- UTR_cal <- add_row(makeChronosCalib(UTR_chronos,age.min = neoaves_min, age.max = neoaves_max, soft.bounds=FALSE), node = 63, age.min = 66, age.max=70,
                                    soft.bounds = FALSE)


#get edge and branch lists

UTR_edge_list <- compute.brlen(UTR_chronos)$edge.length
UTR_branch_list <- branching.times(compute.brlen(UTR_chronos))

#repeat

for(i in 1:1000){
  UTR_chronos_iter <- chronos(UTR_chronos,calibration = UTR_cal, model = "clock")
  UTR_edge_list <- rbind(UTR_edge_list, UTR_chronos_iter$edge.length)
  UTR_branch_list <- rbind(UTR_branch_list, branching.times(UTR_chronos_iter))
}

plot(UTR_chronos_iter)

#choosing median ages


UTR_edge_list <- UTR_edge_list[-1,]
UTR_branch_list <- UTR_branch_list[-1,]

UTR_branches <- colnames(UTR_branch_list)

UTR_medians <- c()
for(i in 1:ncol(UTR_branch_list)){
  UTR_medians[colnames(UTR_branch_list)[i]] <- median(UTR_branch_list[,i])
}

# 
UTR_timetree <- compute.brlen(UTR_chronos)

is.ultrametric(UTR_timetree)
UTR_median_edge <- c()
for(i in 1:ncol(UTR_edge_list)){
  UTR_median_edge<-c(UTR_median_edge,median(UTR_edge_list[,i]))
}
UTR_timetree$edge.length <- UTR_median_edge

plot(UTR_timetree)

write.tree(UTR_timetree, file = "../TimeTrees/UTR_CLCK_medianTT.tree")

  
  