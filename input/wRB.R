
# calculating tree distances

library(phangorn)

#read GTRG tree

GTRG_tree1 <- read.tree("../RAxML_br/GTRG/RAxML_result.ALL50K")
GTRG_tree2 <- read.tree("../RAxML_br/GTRG/RAxML_result.CDS")
GTRG_tree3 <- read.tree("../RAxML_br/GTRG/RAxML_result.lnc_RNA")
GTRG_tree4 <- read.tree("../RAxML_br/GTRG/RAxML_result.other")
GTRG_tree5 <- read.tree("../RAxML_br/GTRG/RAxML_result.UTR")
GTRG_tree6 <- read.tree("../RAxML_br/GTRG/RAxML_result.unannotated")
GTRG_tree7 <- read.tree("../RAxML_br/GTRG/RAxML_result.other")
#GTRG_tree8 <- read.tree("../RAxML_br/GTRG/RAxML_result.pseudogene")

GTRG_cat_tree <- c(GTRG_tree1,GTRG_tree2,GTRG_tree3,GTRG_tree4,GTRG_tree5,GTRG_tree6)

#multiphylo
class(GTRG_cat_tree) = "multiPhylo"

GTRG_result_wrf <- wRF.dist(GTRG_cat_tree, normalize = T)

GTRG_result_wrf

apply(as.matrix(GTRG_result_wrf), 1, mean)



## calibrated 

time_GTRG_tree1 <- read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.ALL50K_timetree")
time_GTRG_tree2 <- read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.CDS_timetree")
time_GTRG_tree3 <- read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.lnc_RNA_timetree")
time_GTRG_tree4 <- read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.other_timetree")
time_GTRG_tree5 <- read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.UTR_timetree")
time_GTRG_tree6 <- read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.unannotated_timetree")
time_GTRG_tree7 <- read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.other_timetree")


time_GTRG_cat_tree <- c(GTRG_tree1,GTRG_tree2,GTRG_tree3,GTRG_tree4,GTRG_tree5,GTRG_tree6)

#multiphylo

class(time_GTRG_cat_tree) = "multiPhylo"

time_GTRG_result_wrf <- wRF.dist(time_GTRG_cat_tree, normalize = T)

time_GTRG_result_wrf


## Across the models

time_GTRG_tree <- read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.ALL50K_timetree")
time_HKY_tree2 <- read.tree("../TimeTrees/HKY85/HKY85_RAxML_result.ALL50K_timetree")
time_JC69_tree <- read.tree("../TimeTrees/JC69/JC69_RAxML_result.ALL50K_timetree")
time_K80_tree <- read.tree("../TimeTrees/K80/K80_RAxML_result.ALL50K_timetree")


models_cat <- c(time_GTRG_tree, time_HKY_tree2, time_JC69_tree,time_K80_tree)

class(models_cat)  = "multiPhylo"

time_models_result_wrf <- wRF.dist(models_cat, normalize = T)

time_models_result_wrf



