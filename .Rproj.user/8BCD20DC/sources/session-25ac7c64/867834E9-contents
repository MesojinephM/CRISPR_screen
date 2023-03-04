#### MARIE MORIARTY SENIOR CAPSTONE PROJECT
#### TEMPORARY STORAGE LOCATION FOR UNUSED/NON-WORKING CODE CHUNKS

### Exploratory plots
# Compare positive and negative group scores
ggplot(mageck, aes(x = neg.score, 
                   y = pos.score, 
                   alpha = 0.05)) +
  geom_point()

ggplot(mageck, aes(x = neg.score)) +
  geom_histogram(bins = 500)

ggplot(mageck, aes(x = pos.score)) +
  geom_histogram(bins = 500)

cor(mageck$neg.score, mageck$pos.score)



ggplot(mageck, aes(x = neg.p_value)) +
  geom_histogram(bins = 500)

ggplot(mageck, aes(x = pos.p_value)) +
  geom_histogram(bins = 500)

cor(mageck$neg.p_value, mageck$pos.p_value)


histogram <- function(var){
  ggplot(mageck, aes(x = var)) +
    geom_histogram(bins = 500) +
    labs(main = var)
}

# Positive histograms
histogram(mageck$pos.score)
histogram(mageck$pos.p_value)
histogram(mageck$pos.fdr)
histogram(mageck$pos.lfc)

# Negative histograms
histogram(mageck$neg.score)
histogram(mageck$neg.p_value)
histogram(mageck$neg.fdr)
histogram(mageck$neg.lfc)

### Separate more significant negative and positive p-values

# I'm trying to split the data set into two based on whether 
# the positive or negative p-value is more significant (smaller). 
# I know I can use a loop to go element by element, but I'm wondering if 
# there is a more efficient way. 

##### Code chunk
# pos.p = NULL
# neg.p = NULL
# # equal.p = NULL
# 
# sort_fun <- function(x){
#   if (x$pos.p_value < x$neg.p_value){
#     pos.p <- c(pos.p, x$pos.p_value)
#   }
#   else{
#     neg.p <- c(neg.p, x$neg.p_value)
#   }
# }
# 
# apply(X = mageck,
#       MARGIN = ,
#       FUN = sort_fun)



# Alternative set using GSEABase
library(GSEABase)
gmt2 <- getGmt("m2.cp.v2022.1.Mm.symbols.gmt")

### How to separate the mageck data frame according to pathways? 
### Do I try to new objects combining the pathways and the mageck results 
### or do I find a way to compare the two sets.