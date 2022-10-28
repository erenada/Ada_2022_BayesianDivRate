#calculating gamma stats and the lineage-throug-time plots

#libs
library(ape)
library(phytools)
library(phangorn)
library(ggtree)
library(ggplot2)
install.packages("TreeSim")
library(TreeSim)

#read tree

getwd()

base_dir <- as.character("../TimeTrees")

model_list <- as.character(list.files(base_dir))

#out_base_dir <- as.character("../TimeTrees")

#tree<-root(read.tree("../TimeTrees/GTRG/GTRG_RAxML_result.ALL50K_timetree"),outgroup = neoaves_outgroup, resolve.root = T)


for(model in model_list) {
  for(tree in list.files(paste(base_dir,model,sep = "/"))) {
    tree_file <- read.tree(paste(base_dir,model,tree,sep = "/"))
    print(paste("Is tree ultrametric?",model,tree,is.ultrametric(tree_file)))
    nnls <- nnls.tree(cophenetic(tree_file), tree_file, rooted = TRUE)
    tree_file <- nnls
    print(paste("Is tree ultrametric?",is.ultrametric(tree_file)))
    write.tree(tree_file, file=paste(base_dir,model,tree,sep="/"))
  }
}


#read all trees and print gamma stats

for(model in model_list) {
  for(tree in list.files(paste(base_dir,model,sep = "/"))) {
    tree_file <- assign(paste(tree), read.tree(paste(base_dir,model,tree, sep="/")))
    print(paste("gamma statistics for", tree, round(gammaStat(tree_file),digits = 3)))
  }
}

class(tree)

gammaStat(GTRG_RAxML_result.ALL50K_timetree)

GTRG_trees <- c(GTRG_RAxML_result.ALL50K_timetree, GTRG_RAxML_result.CDS_timetree,
                    GTRG_RAxML_result.intron_timetree, GTRG_RAxML_result.lnc_RNA_timetree,
                    GTRG_RAxML_result.lnc_RNA_timetree, GTRG_RAxML_result.other_timetree,
                    GTRG_RAxML_result.pseudogene_timetree, GTRG_RAxML_result.unannotated_timetree,
                    GTRG_RAxML_result.UTR_timetree)


class(GTRG_trees) = "multiPhylo"

mltt.plot(GTRG_trees,
         log = "y",
         dcol = FALSE,
         legend = FALSE,
         backward = T)



#yule.trees <- TreeSim::sim.bd.taxa(n=32, numbsim=100, lambda=0.1, mu=0, complete=T)
#bd.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=.9, complete=FALSE)



# depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(GTRG_trees,ape::branching.times)))
# max.depth <- sum(abs(depth.range)) #ape rescales depths
# plot(x=c(0, -1*max.depth), y=c(1, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
# colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
# list.of.both <- list(GTRG_trees, yule.trees)
# for (i in sequence(2)) {
#   tree.list <- list.of.both[[i]]
#   for (j in sequence(length(tree.list))) {
#     ape::ltt.lines(tree.list[[j]], col=colors[[i]])
#   }
# }

# legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)
