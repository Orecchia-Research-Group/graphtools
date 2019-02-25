###############################################################################

###    FuzzyGraph -- 19 Feb 2016 -- Jacob Bor & Katia Oleinik, jbor@bu.edu ###	

###	wd_ART_fuzzyGraph_tasks_prelude, 4 April 2017

###############################################################################

# This code "consolidates" the output from the fuzzy matching.
# Taking a list of unique vertices (exact match IDs) 
# and a matrix of unique edges (pairs of exact match IDs),
# with corresponding weights (similarity scores),
# this code identifies all "clusters", i.e. new patient IDs.

sink("wd_prelude.txt", append=TRUE, split=TRUE)

# set working directory
setwd("/restricted/projectnb/salabs/RecordLinkage/HIV_cohort/alg_output_files")

# clear workspace
rm(list=ls())

# load libraries
library(igraph)
library(data.table)
library(foreign)
library(pryr)

## CHANGED FOR WEIGHTED DIAMETER, SINGLE THRESHOLD / DIAMETER PARAMETER
## CHANGED DEFAULT WEIGHT TO r_jimmy, i.e. column 7.
	#Get command line:
	argv <- commandArgs( TRUE )
	if ( length(argv) > 0 ){
	  threshold <- as.numeric( argv[1] )
	  edge.w <-  as.numeric( argv[2] )
	  ngroups <- as.numeric( argv[3] )
	  
	} else {
	  threshold <- 0.25   
	  edge.w <- 3
	  ngroups <- 50
	}

# turn off scientific notation
options(scipen = 999)

#Get the task ID
#itask <- as.numeric(Sys.getenv("SGE_TASK_ID"))


## create list of unique vertices, format as a list, and add to graph
ref_ids <- read.dta("HIV_EMIDplus_vertices.dta")
vertices <- as.character(as.vector(ref_ids[,1])) ##ADDED "as.character"
vertices <- list(name = vertices)

rm(ref_ids)


## create list of unique edges, format as a matrix, and add to graph
fuzzy.matches <- read.dta("HIV_all_edges_R.dta")
edge.matrix <- as.matrix(fuzzy.matches[,c(1,2)])
## Add weights. CHANGED FOR WEIGHTED DIAMETER, NEGATIVE LOG-WEIGHT, CHANGED FOR r_jimmy
p = exp(-11.13 + 0.366*fuzzy.matches[,edge.w]) / (1 + exp(-11.13 + 0.366*fuzzy.matches[,edge.w]))
w = -log(p)
mem_used()

rm(fuzzy.matches, p)


# create graph
G <- graph.empty(n = 0, directed = F)
G <- add.vertices(G, length(vertices$name), attr=vertices)
G <- add.edges(G, as.character(as.numeric(t(edge.matrix))))
G <- set_edge_attr(G,"weight", value = w)

mem_used()

rm(vertices, edge.matrix, w)


# Create empty graphs to receive output
G.out.A <- graph.empty(n = 0, directed = F)
G.out.B <- graph.empty(n = 0, directed = F)



### Eliminate very low-ranked edges and break up large graph components

# Delete edges below the threshold, p = 0.25
G <- delete.edges(G,which(E(G)$weight > -log(0.25)))

# Identify clusters in G
cl <- clusters(G)


# Identify and break up large clusters
big.cls <- which(cl$csize>=5000)
G.small <- delete.vertices(G,which(cl$membership %in% big.cls))
G.big <- delete.vertices(G,which(!(cl$membership %in% big.cls)))

mem_used()

#save.image(file = "debug prelude.RData")



#rm(G)
G.big.exists <- 1
loop.threshold <- 0.25
while (G.big.exists ==1) {
  loop.threshold <- loop.threshold + 0.05
  G.big <- delete.edges(G.big,which(E(G.big)$weight > -log(loop.threshold)))
  cl.big <- clusters(G.big)
  big.cls <- which(cl.big$csize>=5000)
  G.big.small <- delete.vertices(G.big,which(cl.big$membership %in% big.cls))
  G.big <- delete.vertices(G.big,which(!(cl.big$membership %in% big.cls)))
  G.small <- union(G.small, G.big.small)
  E(G.small)$weight <- pmin(E(G.small)$weight_1, E(G.small)$weight_2, na.rm=TRUE)
  G.small <- remove.edge.attribute(G.small, "weight_1")
  G.small <- remove.edge.attribute(G.small, "weight_2")
  G.big.exists <- length(V(G.big)) > 0
}

print(paste("   The last look threshold was ",loop.threshold))


#check same V, fewer E
print(paste("   The number of vertices in G.small is ",length(V(G))))
print(paste("   This should equal the number in G, ",length(V(G.small))))
print(paste("   The number of edges in G.small is ",nrow(get.edgelist(G))))
print(paste("   This should be less than the number of edges in G, ",nrow(get.edgelist(G.small))))


# G.old <- G # this is just for trouble shooting; drop this 
G <- G.small
rm(G.small)

rm(cl.big, G.big, G.big.exists, G.big.small, big.cls, loop.threshold)


# Identify new clusters in G (now revised to eliminate problematic edges)
cl <- clusters(G)
print(paste("   The number of clusters is ", cl$no))
print(paste("   The max cluster size is ", max(cl$csize)))

# OUTER LOOP: LOOP OVER PARTITIONS OF THE FULL GRAPH

temp.clust.ids <- unique(cl$membership)
group.size <- ceiling ( length( temp.clust.ids ) / ngroups )

#randomize clusters
set.seed(6819)
sample.clust.ids <- sample ( temp.clust.ids, length( temp.clust.ids ) )

rm(temp.clust.ids)

save.image(file = "wd_HIV_fuzzyGraph_tasks_prep.RData")


# 27,303,301 clusters, all with <5000 nodes
# 62,808,665 vertices
# 83,026,988 edges, reduced from 305,576,121 original edges

