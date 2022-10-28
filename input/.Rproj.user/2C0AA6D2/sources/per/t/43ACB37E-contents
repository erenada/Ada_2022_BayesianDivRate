#looped calibration 

#loop

# base directory (model folders)

base_dir <- as.character("../RAxML_br")

model_list <- as.character(list.files(base_dir))

out_base_dir <- as.character("../TimeTrees")

neoaves_min <- 66.0
neoaves_max <- 99.6 
neoaves_outgroup <- c('StrCam')


for(model in model_list) {
  for(tree in list.files(paste(base_dir,model,sep = "/"))) {
  in_tree <- root(read.tree(paste(base_dir,model,tree, sep="/")),outgroup = neoaves_outgroup, resolve.root = T)
  cal <- add_row(makeChronosCalib(in_tree,age.min = neoaves_min, age.max = neoaves_max, soft.bounds=FALSE), node = c(35,34), age.min = c(66,50.03), age.max=c(70,50.03),soft.bounds = FALSE)
  edge_list <- compute.brlen(in_tree)$edge.length
  branch_list <- branching.times(compute.brlen(in_tree))
  
    for(i in 1:100){
    chronos_iter <- chronos(in_tree,calibration = cal, model = "clock")
    edge_list <- rbind(edge_list, chronos_iter$edge.length)
    branch_list <- rbind(branch_list, branching.times(chronos_iter))
    }
  
  edge_list <- edge_list[-1,]
  branch_list <- branch_list[-1,]
  branches <- colnames(branch_list)
  medians <- c()
  
  for(i in 1:ncol(branch_list)) {
    medians[colnames(branch_list)[i]] <- median(branch_list[,i])
  }
  
  timetree <- compute.brlen(in_tree)
  median_edge <- c()
  
  for(i in 1:ncol(edge_list)) {
    median_edge<-c(median_edge,median(edge_list[,i]))
  }
  
  timetree$edge.length <- median_edge
  
  tree_name <- paste(model,tree,c("timetree"),sep="_")
  
  out_model_dir <- paste(out_base_dir,model,sep = "/")
  
  write.tree(timetree, file = paste(out_model_dir,tree_name,sep="/"))
  }
}

#assign(paste(model,tree, sep="_"), read.tree(paste(base_dir,model,tree, sep="/")))


