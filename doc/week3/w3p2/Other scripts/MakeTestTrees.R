library(ape)
library(apTreeshape)
setwd("C:/TESTTREES/")

NUMTREES = 10000
ntips = 512

for (t in 1:NUMTREES) 
{
  message(t)
  curr_tree = rtreeshape(1,ntips,model="yule")[[1]]
  curr_tree_phylo = as.phylo(curr_tree)
  treefile_name = paste("testtree_yule_n", as.character(t),".newick",sep="")
  write.tree(curr_tree_phylo,treefile_name)
}

