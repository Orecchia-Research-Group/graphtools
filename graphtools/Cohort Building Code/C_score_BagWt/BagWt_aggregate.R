## This script takes the results from the parallelized bootstrap process (which creates 500 text files) and aggregates them into final results

## Input: 500 individual results from Bootstrap runs of the algorithm
##    - Located in the path below and named optResults_[001-500].txt
## Output: Two files
##    - "all_bootstrapped_weights.csv" collects all files above and creates one csv file
##    - "bootstrap_aggregation_results.csv" contains the final results, i.e. a simple average of bootstrapped weights for each category

## Set working directory
#setwd("/restricted/projectnb/salabs/RecordLinkage/jpotter/OptOutput/")

setwd("/restricted/projectnb/salabs/RecordLinkage/HIV_cohort/alg_code/JimmyBagWt/")


## Create empty matrix for weights
optres <- as.data.frame(matrix(NA,nrow=500,ncol=6))

## Loop 1-500, reading in each result, formatting and adding to data frame
for(i in 1:500){
  
  ## Format the path for each result file - formatC() used to pad zeroes for formatting
  result.path <- paste("optResults_",formatC(i,width=3,format="d",flag="0"),".txt",sep="")
  
  ## Read in test result
  test.result <- read.csv(result.path,header=FALSE)
  
  ## Place result in matrix
  # note: the text file outputs from BagWt_tasks.R have 2 rows for some reason, which is why we have this wonky code here to read it in. Could fix this at some point (in BagWt_tasks.R), but it works.
  optres[i,] <- as.matrix(c(test.result[1,],test.result[2,1]))
}

## Assign names to the columns
names(optres) <- c("DOB", "first", "last", "gender", "province", "facility")

## Save weights together as a .csv file
write.csv(optres, file="all_bootstrapped_weights.csv")

## Calculate the average weights (which is the bootstrap aggregation procedure)
average.weights <- apply(optres,2,mean)

## Assign names to the columns
names(average.weights) <- c("DOB", "first", "last", "gender", "province", "facility")

## Results of Bootstrap Aggregation
average.weights

## Save results to a file
write.csv(average.weights, file="bootstrap_aggregation_results.csv")
