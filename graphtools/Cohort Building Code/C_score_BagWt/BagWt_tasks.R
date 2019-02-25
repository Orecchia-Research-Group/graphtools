## Weight Optimization

## Summary: This script conducts 1 bootstrap (by Ref Episode ID) of sim scores to identify optimal weights
## Input: Data file with Ref Ep ID, Ep ID, Manual Match Result, and the 6 raw sim scores
## Output: One bootstrapped and optimized set of weights

## NOTE: this script does not bootstrap by itself. It is coded to be run by a .qsub process that assigns a task ID to each core running the code
## .qsub file = "weightJobv2.qsub"

setwd("/restricted/projectnb/salabs/RecordLinkage/HIV_cohort/alg_code/JimmyBagWt/")

# Load Packages to read DTA file
library(readstata13)

# Pull in the task number (which is a unix environmental variable assigned to each of 500 parallel processes in the .qsub file)
id <-as.numeric(Sys.getenv("SGE_TASK_ID"))

#### Load Data file
raw.data <- read.dta13("/restricted/projectnb/salabs/RecordLinkage/HIV_cohort/alg_input_files/Training_Data_forBAGGING.dta")

#### Calculate the number of ref episodes to bootstrap (currently total number of ref episodes)
samp.size <- length(unique(raw.data$ref_episode_no))

#### Create a list of all the unique ref episodes for bootstrap drawing
ref.eps <- unique(raw.data$ref_episode_no)

#### Create Custom Function "quasi.c"
#### Input: Vector of Sim Score weights (length 6)
#### (Uses the data frame "dd", required to be already loaded, which contains the matching dataset)
#### Output: Quasi C Score for dataset dd, given the input weights

quasi.c <- function(inp.weights, dd){

  ## Input vector will be inp.weights, length 6 vector
  ## Example: inp.weights <- c(.5,.5,.5,.5,.5,.5)
  ## Calculate exp(weights) to ensure that weights are non-negative
  exp.inp.weights <- exp(inp.weights)
  ##exp.inp.weights <- inp.weights

  ## Calculate the weighted Tot Sim Score from the 6 components of Tot Sim Score
  dd$w.tot.sim.score <- exp.inp.weights[1] * dd$DOB_sim_score +
    exp.inp.weights[2] * dd$first_sim_score +
    exp.inp.weights[3] * dd$last_sim_score +
    exp.inp.weights[4] * dd$gender_sim_score +
    exp.inp.weights[5] * dd$prov_sim_score +
    exp.inp.weights[6] * dd$facility_sim_score

  ## n_true_matches == total # of man matched positive values (with non-missing sim score)
  n.matches <- sum(dd$r_manmatch)

  ## Get the percentiles (100 bins) of the tot_sim_score within the man-match positives
  bins <- quantile(dd[dd$r_manmatch == 1,]$w.tot.sim.score,seq(0,.99,.01))

  ## Get the count of the number that are above each percentile value (True Positives) within the man-match positives
  n.tp <- sapply(bins, function(x) sum(dd[dd$w.tot.sim.score > x,]$r_manmatch))

  ## Calculate the sensitivy for each percentile as the above count divided by total true matches
  sens.bins <- n.tp/n.matches

  ## Take the mean value of r_manmatch (man-match) for all values above the percentile value as the PPV
  ppv.bins <- sapply(bins, function(x) mean(dd[dd$w.tot.sim.score > x,]$r_manmatch))

  ## Calculate the quasi-C-statistic (ROC curve area) as the sum of all of the PPV values
  quasi.c.stat <- sum(ppv.bins)

  ## Output of function is the quasi-C-statistic (area under curve)
  ## Note: is negative so that optim will minimize this output
  return(-quasi.c.stat)
}



## Set seed for parallel computation (independent random seeds for each parallel processor)
set.seed( id )
  
  #### Grab the ref episode numbers to bootstrap
  draws <- sample(1:samp.size,samp.size, replace=TRUE)

  #### Use the random draws to generate a new data frame (dd) with each selected set of ref episodes
  dd <- NULL
  for(j in 1:samp.size){
   dd <- rbind(dd, raw.data[raw.data$ref_episode_no == ref.eps[draws[j]],])
  }

  
  #### Optimize, to find vector of weights that minimizes the (neg) quasi-c statistic
  opt.out <- optim(par=c(0,0,0,0,0,0), quasi.c, , dd=dd, control=list(trace=TRUE))

  print("Parameters from latest loop are ")
  print(opt.out$par)
  

## Save output from loop to an RDS file
write( exp(opt.out$par), file=paste0("optResults_", sprintf("%03d",id),".txt"), sep=",")


