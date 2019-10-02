#load related, the package we will be using. 
library(related)
library(tidyverse)
library(ggpubr)

#load in the properly formatted .txt file containing all the genotypes; also load in the parent edge list, which should have a column of dyad names called "dyads". 
Genotypes <- readgenotypedata("INSERTGENOTYPESHERE.txt")
Edges <- read_csv("INSERTPARENTEDGESHERE.csv")

#generate an output file with relatedness estimator values for all possible pairs of individuals. In the coancestry() command, you will need to specify which estimators you want to be generated, and when you subset the data later in this script, you will need to specify which relatedness indicator columns to include, depending on which you generate at this step. 
Output <- coancestry(Genotypes$gdata, wang = 2)

#we need a column in the output file that contains the name of the dyad, in the format "Individual1ID_Individual2ID" so that we can compare this to the parent edge list, and assign a status of MATED or UNMATED to each dyad in the relatedness output. 
RelatednessValues <- unite(Output$relatedness, "dyads", c(ind1.id, ind2.id), remove = FALSE)

#now we create a column in the output that tells us if a dyad in the output file is present in the dyads column of the edge list (so basically a column that tells us if a particular dyad is known to have reproduced together or not).
RelatednessValues$mated <- RelatednessValues$dyads %in% Edges$dyads

#we also want a column in the output that tells us if a dyad is mixed sex or single sex (i.e. if the row contains the relatedness value for a pair of individuals that were potential mates, or the relatedness value for two males/two females.)
RelatednessValues$MixedSex <- ifelse(grepl("M", RelatednessValues$dyads) & grepl("F", RelatednessValues$dyads), "YesMixedSex", "NoSingleSex")

# Now we need to extract all the rows where MixedSex = YesMixedSex, so that we can compare the relatedness values for individuals that mated (i.e. where Mated = TRUE) to those of individuals who could have but did not mate (i.e. where Mated = FALSE). If you generated more relatedness indicators earlier in this script, you will need to add them to the subsetting here. 
RelatednessOfPotentialMates <- subset(RelatednessValues, MixedSex == "YesMixedSex", select = c("dyads", "ind1.id", "ind2.id", "wang", "mated"))

# We will also extract all the rows where mated = TRUE, so that we can compare the relatedness values for individuals that mated (i.e. where Mated = TRUE) to a mock population created by randomly sampling all potential mated pairs. 
RelatednessOfRealMatedPairs <- subset(RelatednessValues, mated == TRUE, select = c("dyads", "ind1.id", "ind2.id", "wang", "mated"))

# Bootstrapping: We need to randomly sample XXX pairs (change this number to match the number of actual mated pairs in your dataset) from all possible male-female dyads, do that a large number of times, and  calculate mean relatedness all those model populations, and turn that numeric vector into a dataframe, so you can make a distribution of means.  
DistributionOfModelMeans <- replicate(n = 10000, expr = mean(sample(RelatednessOfPotentialMates$wang, XXX, replace = TRUE)))
DistributionOfModelMeans <- as.data.frame(DistributionOfModelMeans)

# Then calculate mean relatedness of actual mated pairs.
MeanOfRealPairs <- mean(RelatednessOfRealMatedPairs$wang)

# Now compare the real mean to the bootstrapped means- is it more than one standard deviation from the model population mean?
Bootstrapping <- ggplot(data = DistributionOfModelMeans, mapping = aes(x = DistributionOfModelMeans)) + geom_histogram(mapping = aes(y = (..count../sum(..count..))), binwidth = 0.008) + geom_vline(xintercept = MeanOfRealPairs, size = 1, colour = "#FF3721", linetype = "dashed") + ggtitle("Null Distribution of Mean Relatedness in POPULATION") + labs(x = "Mean Relatedness of Resample", y = "Proportion of Resamples") + theme_bw()
plot(Bootstrapping)







