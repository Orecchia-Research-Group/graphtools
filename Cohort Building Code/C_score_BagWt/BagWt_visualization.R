library(ggplot2)
library(GGally)
library(RColorBrewer)

## Set working directory
##setwd("/restricted/projectnb/salabs/RecordLinkage/jpotter/OptOutput/")
setwd("/restricted/projectnb/salabs/RecordLinkage/HIV_cohort/alg_code/JimmyBagWt/")


# Save weights to file
read.csv(optres, file="all_bootstrapped_weights.csv")

average.weights <- apply(optres,2,mean)
names(average.weights) <- c("DOB", "first", "last", "gender", "province", "facility")
average.weights

optres <- rbind(optres, average.weights)

optres$type <- "Single Bootstrap Run"
optres[501,7] <- "Mean"
optres$type <- as.factor(optres$type)
names(optres) <- c("DOB", "first", "last", "gender", "province", "facility", "type")

pairs(optres[,1:6])

## Using ggplot2 and GGally packages
ggpairs(data=optres, mapping=aes(color=type), columns=1:6, diag=NULL, upper=list(continuous="points"))

## Make a color scale for use in the scatterplots
myColors <- brewer.pal(3,"Set1")
names(myColors) <- levels(optres$type)
colScale <- scale_colour_manual(name="type",values=myColors)

## Make a custom function tuned for the scatterplot we want to duplicate into matrix format
#customplot <- function(data,mapping, ..., low="#132B43", high="56B1F7"){
#  ggplot(data=data,mapping=mapping) + geom_point() + scale_colour_grey() +
#    scale_alpha_discrete(range=c(.99,.6))
#}
#
#ggpairs(data=optres, mapping=aes(color=type, alpha=type), columns=1:6,diag=NULL,
#        upper=list(continuous=customplot),
#        lower=list(continuous=customplot))


## Make a custom function tuned for the scatterplot we want to duplicate into matrix format
customplot <- function(data,mapping, ..., low="#132B43", high="56B1F7"){
  ggplot(data=data,mapping=mapping) + geom_point() + scale_colour_grey() +
    scale_alpha_discrete(range=c(.99,.05)) +
    #scale_x_continuous(limits=c(0.2,2)) +
    #scale_y_continuous(limits=c(0.2,2)) +
    theme_bw()
}

ggpairs(data=optres, mapping=aes(alpha=type), columns=1:6,
        #diag=NULL,
        upper=list(continuous=customplot),
        lower=list(continuous=customplot))

## This visual works but needs to be photoshopped to remove extra
## plots (lower half) and move axes from bottom and left to proper positions
ggsave("bootstrap_visual_v2_needs_photoshopping.png")
