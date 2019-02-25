
## I USE WEIGHTED DIAMETER HERE, WITH WEIGHTS EQUAL TO  p = exp(-11.45 + 0.366*r_jimmy) / (1 + exp(-11.45 + 0.366*r_jimmy))
## IN MAKING THIS REVISION, I NEED TO: 
## (A) REPLACE THE DIAMETER CONSTRAINT WITH A WEIGHTED DIAMETER CONSTRAINT, NOTING THAT HIGHER WEIGHTED DIAMETER IS GOOD (NOT BAD)
## (B) REVISE THE "BIG CLUSTER" / "SMALL CLUSTER" CONSTRAINTS SO THAT I DON'T ALWAYS RULE OUT (OR IN) BIG (SMALL) CLUSTERS
## (C) CALCULATE THE NEW WEIGHTS, AND CALCULATE WEIGHTED DIAMETER (USING LOG TRANSFORMATION)


#1) additional degrees of separation imply further distance; but higher weights indicate stronger edges and thus shorter distance.Therefore, we need to invert the rank ordering somehow (e.g., put a minus sign in front)
#2) from an interpretability standpoint, we thought we might translate the weights into pseudo-probabilities; and that these probabilities would thenbe multiplied together. I.e., distance between two nodes say 1 and 3 is the product of the edges: w12*w23. This implies a multiplicative model but R uses addition. We could solve this by implement on the log scale: exp(log(w12)+log(w23))
#So, this suggests the following approach: convert existing weights into probabilities. Create the weights as log-probabilities. Let R add them up. Then use log(threshold-probability) as the threshold.


###############################################################################

###    FuzzyGraph -- 19 Feb 2016 -- Jacob Bor & Katia Oleinik, jbor@bu.edu ###	diameter = WEIGHTED

###############################################################################

# This code "consolidates" the output from the fuzzy matching.
# Taking a list of unique vertices (exact match IDs) 
# and a matrix of unique edges (pairs of exact match IDs),
# with corresponding weights (similarity scores),
# this code identifies all "clusters", i.e. new patient IDs.
# We execute the process iteratively, starting with a low but
# plausible threshold, identifying clusters, then applying simple
# "rules" about acceptable network structure, e.g. weighted diameter (probability)
# mulst be greater than some threshold. If not, drop the lowest scoring edge. 

  # Input parameters that could be changed in future versions:
  # -- Weighted diameter (probability) of acceptable cluster (currently >=0.5)


# To run:
#   in terminal: Rscript fuzzyGraph_2_19feb2016.R > fuzzyGraph_2.log &
#   qsub: qsub...
#   start_tasks.sh to optimize



# set working directory
setwd("/restricted/projectnb/salabs/RecordLinkage/HIV_cohort/alg_output_files")

# clear workspace
rm(list=ls())

# load libraries
library(igraph)
library(data.table)
library(foreign)
library(pryr)

# load "prep" workspace with full graph already loaded (moved from below)
load(file="wd_HIV_fuzzyGraph_tasks_prep.RData")


#Get command line:
argv <- commandArgs( TRUE )
if ( length(argv) > 0 ){
  threshold <- as.numeric( argv[1] )
  edge.w <-  as.numeric( argv[2] )
  ngroups <- as.numeric( argv[3] )
  
} else {
  threshold <- 0.5   
  edge.w <- 3
  ngroups <- 50
}

# turn off scientific notation
options(scipen = 999)

#Get the task ID
itask <- as.numeric(Sys.getenv("SGE_TASK_ID"))


# ### COMMENTED OUT FOR VERSION WITH PRELUDE / TASKS

  # ## create list of unique vertices, format as a list, and add to graph
  # ref_ids <- read.dta("/restricted/projectnb/salabs/RecordLinkage/Feb2016/Full275M/Full275_EMIDplus_CD4VLworkup.dta")
  # vertices <- as.vector(ref_ids[,1])
  # vertices <- list(name = vertices)
  # rm(ref_ids)
  # 
  # 
  # ## create list of unique edges, format as a matrix, and add to graph
  # fuzzy.matches <- read.dta("/restricted/projectnb/salabs/RecordLinkage/Feb2016/Full275M/Full275_CD4VLworkup_edges_R3.dta")
  # edge.matrix <- as.matrix(fuzzy.matches[,c(1,2)])
  # ## Add weights. CHANGED FOR WEIGHTED DIAMETER, NEGATIVE LOG-WEIGHT, CHANGED FOR r_jimmy
  # p = exp(-11.45 + 0.366*fuzzy.matches[,edge.w]) / (1 + exp(-11.45 + 0.366*fuzzy.matches[,edge.w]))
  # 
  # mem_used()
  # 
  # rm(fuzzy.matches)
  # 
  # 
  # # create graph
  # G <- graph.empty(n = 0, directed = F)
  # G <- add.vertices(G, length(vertices$name), attr=vertices)
  # G <- add.edges(G, as.character(as.numeric(t(edge.matrix))))
  # G <- set_edge_attr(G,"weight", value = -log(p))
  # 
  # mem_used()
  # 
  # rm(vertices, edge.matrix)
  # 
  # 
  # # Delete edges below the threshold, and replace G
  # G <- delete.edges(G,which(E(G)$weight > -log(threshold)))
  # 
  # # Create empty graphs to receive output
  # G.out.A <- graph.empty(n = 0, directed = F)
  # G.out.B <- graph.empty(n = 0, directed = F)
  # 
  # # Identify clusters in G
  # cl <- clusters(G)
  # 

# OUTER LOOP: LOOP OVER PARTITIONS OF THE FULL GRAPH (itask = 1/ngroups is called from qsub)
temp.clust.ids <- unique(cl$membership)
group.size <- ceiling ( length( temp.clust.ids ) / ngroups )

#randomize clusters
set.seed(6819)
sample.clust.ids <- sample ( temp.clust.ids, length( temp.clust.ids ) )

# sample clusters for partition for this task 
i.start <- (itask - 1) * group.size + 1
i.end <- min( itask * group.size, length(sample.clust.ids) )

G.i.group <- delete.vertices(G,which(!(cl$membership %in% (sample.clust.ids[i.start:i.end]))))
G <- G.i.group  
rm(G.i.group) #cleanup

#generate V and E totals for the itask sample of the full graph.
total.n.vertices <- length(V(G))
total.n.edges <- nrow(get.edgelist(G))

sum.edges.deleted <- 0


mem_used()

  
# INNER LOOP: EACH ITERATION, DROP THE LOWEST VALUED EDGES FROM EACH IMPLAUSIBLE CLUSTER  

h = 0  

len.V.G <- length(V(G))

while (len.V.G > 0) {

  h <- h + 1 
  print(paste("----- Loop = ", h, "-----"))
  
  # Identify clusters in G
  cl <- clusters(G)
  
  # extract small clusters, size <=2, save in temp graph
  small.size <- 2
  small.cls <- which(cl$csize <= small.size)
  G.out.A.temp <- delete.vertices(G,which(!(cl$membership %in% small.cls)))
  
  # append G.out.A.temp to G.out.A, (i.e., the final output graph that will collect the final set of clusters [which we can also think of as the final set of vertices and edges])
  v.G.out.A.temp <- list(name=V(G.out.A.temp)$name)
  e.G.out.A.temp <- get.edgelist(G.out.A.temp)
  w.G.out.A.temp <- as.vector(E(G.out.A.temp)$weight)
  G.out.A <- add.vertices(G.out.A,length(v.G.out.A.temp$name),attr=v.G.out.A.temp)
  G.out.A <- add.edges(G.out.A,as.character(as.numeric(t(e.G.out.A.temp))),weight=w.G.out.A.temp)
  
  #extract remaining clusters, and put in graph "G.out.remainder"
  #small.big.cls <- c(small.cls, big.cls)
  G.out.remainder <- delete.vertices(G,which(cl$membership %in% small.cls))
  
  # identify new.cl cluster list
  new.cl <- clusters(G.out.remainder)
  
  # define number of vertices and edges in G.out.remainder 
  v.G.out.remainder <- length(V(G.out.remainder))
  e.G.out.remainder <- nrow(get.edgelist(G.out.remainder))
  
  # create output holders and counters for the loop that follows
  ## NOTE: SUBSTITUTED 0 FOR 0L TO TRY TO SOLVE PRECISION PROBLEM
  v.out.a <- rep(0,total.n.vertices)
  e.out.a <- matrix(0,nrow=total.n.edges,ncol=3)
  v.out.b <- rep(0,total.n.vertices)
  e.out.b <- matrix(0,nrow=total.n.edges,ncol=3)
  
  # start at 0 - we do not have any edges or vertices yet
  v.counter.a <- 0
  e.counter.a <- 0
  v.counter.b <- 0
  e.counter.b <- 0
  edges.deleted <- 0

  ## decompose G.out.remainder into components; calculate WEIGHTED diameter of each component; and move PLAUSIBLE components to A, IMPLAUSIBLE components to B
  
  print("Start loop over cluster strata")
  
  # loop over strata of clusters by cluster size (solves non-linear processing time issue with big graphs)  
  for (j in sort(unique(new.cl$csize))){

    #print("")
    print(paste("Cluster size = ",j))
      
    j.cls <- which(new.cl$csize==j)
    G.j <- delete.vertices(G.out.remainder,which(!(new.cl$membership %in% j.cls)))
    
    #print("")
    #print(paste("   number of vertices in G.j = ", length(V(G.j))))
    #G.j.n.edges <- nrow(get.edgelist(G.j))
    #print(paste("   number of edges in G.j = ", G.j.n.edges))
    
    if (length(V(G.j)) > 0) {
      
      j.components <- decompose(G.j)
      
      v.in.j.components <- 0 # number of vertices in j.components, calculated as a check
      e.in.j.components <- 0 # number of edges in j.components, calculated as a check
      
      for (k in 1:length(j.components)) {

        v.in.j.components <-  v.in.j.components + length(V(j.components[[k]]))
        e.in.j.components <- e.in.j.components + nrow(get.edgelist(j.components[[k]]))
        
        ## Is cluster plausible? <100 vertices? weighted diameter < threshold?
        ## Note: threshold cluster size 100 rarely binds and is only used to avoid calculating diameter of large clusters.
        plausible <- j < 100
        if (plausible==1) {
          plausible <- (diameter(j.components[[k]]) <= -log(threshold)) ##weighted: sum (log_p), removed ", weights=NA"
          }
          
        if (plausible==1) {
          e.mat.a <- get.edgelist(j.components[[k]])
          v.vec.a <- list(name=unique(as.vector(e.mat.a)))$name
          e.mat.a <- cbind(e.mat.a,E(j.components[[k]])$weight)
          
          v.out.a[ (v.counter.a+1) : ( v.counter.a + length(v.vec.a) )] <- v.vec.a
          e.out.a[ (e.counter.a+1) : ( e.counter.a + nrow(e.mat.a) ), ] <- e.mat.a
          
          v.counter.a <- v.counter.a + length(v.vec.a)  
          e.counter.a <- e.counter.a + nrow(e.mat.a)  
          
        }  else {
          
          e.mat.b <- matrix(get.edgelist(j.components[[k]]),ncol=2)
          v.vec.b <- list(name=unique(as.vector(e.mat.b)))$name
          e.mat.b <- matrix(cbind(e.mat.b,E(j.components[[k]])$weight),ncol=3)
          n.edge.b <- nrow(e.mat.b)
          
          ##DROP LOWEST 10% OF SCORED EDGES FROM e.mat.b, sort descending by weight, and then by random number to break ties
          e.mat.b <- e.mat.b[order(e.mat.b[,3],runif(nrow(e.mat.b)),decreasing=TRUE),]   #sort low to high
          e.mat.b <- matrix(e.mat.b[ (ceiling(nrow(e.mat.b)/10)+1) : nrow(e.mat.b),],ncol=3)	    #drop lowest 10%, rounding up
          n.edge.b.trimmed <- nrow(e.mat.b)
          
          v.out.b[ (v.counter.b+1) : ( v.counter.b + length(v.vec.b) )] <- v.vec.b
          e.out.b[ (e.counter.b+1) : ( e.counter.b + nrow(e.mat.b) ), ] <- matrix(e.mat.b,ncol=3)
          
          v.counter.b <- v.counter.b + length(v.vec.b)  
          e.counter.b <- e.counter.b + n.edge.b.trimmed  
          edges.deleted <- edges.deleted + (n.edge.b - n.edge.b.trimmed)
          
        } # close diameter if statement
      } # close loop over k subcomponents of j.components
    }  # close if statement that graph G.j is not emptygraph
  }  # close loop over strata of clusters of size j 

  print(" ")
  print("Summary of loop h")
  print(" ")
  print(paste("   We have now completed the G.out.remainder loop for round h = ",h))
  print(paste("   Cumulative number of vertices put into A is ",v.counter.a))
  print(paste("   Cumulative number of vertices put into B is ",v.counter.b))
  print(paste("   vertices A + B= ",v.counter.a + v.counter.b, " of total n of vertices:", v.G.out.remainder))
  
  print(" ")
  print(paste("   Cumulative number of edges put into A is ",e.counter.a))
  print(paste("   Cumulative number of edges put into B is ",e.counter.b))
  print(paste("   Edges A + B = ",e.counter.a + e.counter.b, " of total n of edges:", e.G.out.remainder))
  print(paste("   Edges A + B + total deleted = ",e.counter.a + e.counter.b + edges.deleted, " of total n of edges:", e.G.out.remainder))
  
  #katia: save all the data here, so we can start debugging from here:
  #save.image(file = paste0("debug/debug_task_",itask,".RData"))
  #load(file="debug_task_1.RData")
    
  # trim output matrix/vectors from G.out.remainder loop
  e.out.a <- e.out.a[1 : e.counter.a, ]
  v.out.a <- v.out.a[1 : v.counter.a ]
  v.list.a <- list(name = v.out.a)
  
  e.out.b <- matrix(e.out.b[1 : e.counter.b, ],ncol=3)
  v.out.b <- v.out.b[1 : v.counter.b ]
  v.list.b <- list(name=v.out.b)
  
  # append to G.out.A, G.out.B
  if (sum(as.numeric(v.list.a$name))>0) {
    G.out.A <- add.vertices(G.out.A, length(v.list.a$name), attr=v.list.a)
  }
  if (sum(as.numeric(e.out.a)) >0){
    G.out.A <- add.edges(G.out.A, as.character(as.numeric(t(e.out.a[,1:2]))),weight=e.out.a[,3])
  }

  G.out.B <- graph.empty(n=0,directed=F)
  if (sum(as.numeric(v.list.b$name))>0) {
    G.out.B <- add.vertices(G.out.B, length(v.list.b$name), attr=v.list.b)
  }
  if (sum(as.numeric(e.out.b)) >0){
    G.out.B <- add.edges(G.out.B, as.character(as.numeric(t(e.out.b[,1:2]))),weight=e.out.b[,3])
  }
  
  sum.edges.deleted <- sum.edges.deleted+edges.deleted
  
  # Summary
  print(" ")
  print(paste("   Summary of aggregate progress from loops h = 1 through h = ", h))
  print(" ")
  print(paste("   G.out.A contains n = ",length(V(G.out.A)), " vertices."))
  print(paste("   G.out.B contains n = ",length(V(G.out.B)), " vertices."))
  print(paste("   The sum of vertices in A + B", length(V(G.out.A)) + length(V(G.out.B)), " should equal ",total.n.vertices))
  print(paste("   The sum of edges in A + B + deleted", nrow(get.edgelist(G.out.A)) + nrow(get.edgelist(G.out.B)) + sum.edges.deleted, " should equal ",total.n.edges))
  print(paste("   The total number of clusters A + B is", no.clusters(G.out.A) + no.clusters(G.out.B)))
  print(" ")
  # re-define G.out.B as G, to loop over again
  G <- G.out.B
  
  # added to break out of loop if we reach the end (i.e. no more implausible graphs to break up)
  len.V.G <- length(V(G))
  #if (length(V(G))==0) break

}  # END OF ROUND "h" LOOP

# final summary of G.out.A
#no.clusters(G.out.A)
#length(V(G.out.A))

save(G.out.A, file=paste0("wd_G.out.final.R/G.out.final.A_",itask,".Rdata"))
#save(G.out.B, file=paste0("wd_G.out.final.R/G.out.final.B_",itask,".Rdata"))


mem_used()
     
# Export to Stata
cl.A  <- clusters(G.out.A)
#cl.B  <- clusters(G.out.B)
df.G.out.A <- data.frame(ref_id=V(G.out.A)$name, new_id=cl.A$membership)
#df.G.out.B <- data.frame(ref_id=V(G.out.B)$name, new_id=cl.B$membership)

# SAVE WITH "i"
# G_out_final_A_d`d'_t`t_low'_`t_high'_`w'_1.dta
write.dta(df.G.out.A, file=paste0( "wd_G_out_final_dta/G_out_final_A_wd",
                                   threshold,"_",
                                   edge.w,"_",
                                   itask,".dta"), 
				   convert.factors = "numeric")
#write.dta(df.G.out.B, file=paste0( "wd_G_out_final_dta/G_out_final_B_wd",
#                                   threshold,"_",
#                                   edge.w,"_",
#                                   itask,".dta"), 
#				    convert.factors = "numeric")

#sink()
