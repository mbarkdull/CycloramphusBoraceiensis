# load the packages we need:
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(related)
  library(tidygraph)

# Generate relatedness estimators for all dyads:
  # Load in the properly formatted .txt file containing all the genotypes; also load in the parent edge list, which should have a column of dyad names called "dyads". 
    Genotypes <- readgenotypedata("INSERTGENOTYPESHERE.txt")
    Edges <- read_csv("INSERTPARENTEDGESHERE.csv")

# Generate an output file with relatedness estimator values for all possible pairs of individuals. In the coancestry() command, you will need to specify which estimators you want to be generated, and when you subset the data later in this script, you will need to specify which relatedness indicator columns to include, depending on which you generate at this step. 
  Output <- coancestry(Genotypes$gdata, wang = 2)

# We need a column in the output file that contains the name of the dyad, in the format "Individual1ID_Individual2ID" so that we can compare this to the parent edge list, and assign a status of MATED or UNMATED to each dyad in the relatedness output. 
  RelatednessValues <- unite(Output$relatedness, "dyads", c(ind1.id, ind2.id), remove = FALSE)

# We also want a column in the output that tells us if a dyad is mixed sex or single sex (i.e. if the row contains the relatedness value for a pair of individuals that were potential mates, or the relatedness value for two males/two females.)
  RelatednessValues$MixedSex <- ifelse(grepl("M", RelatednessValues$dyads) & grepl("F", RelatednessValues$dyads), "YesMixedSex", "NoSingleSex")

# Now, we first need to identify individuals with more than one mate, because we are interested in looking at the relatedness of indviduals that share a mate. 
  # Load in the node and edge data, and construct a network:
    Nodes <- read_csv("INSERTPARENTNODESHERE.csv")
    Edges <- read_csv("INSERTPARENTEDGESHERE.csv")
    Network <- tbl_graph(nodes = Nodes, edges = Edges)

# Now get the degree of each individual, and subset list of individuals to only those with more than one mate (i.e. degree > 1), by their vertex number.
  Nodes$degree <-degree(Network)
  Nodes$Vertex <- V(Network)
  IndividualsWithMultipleMates <- subset(Nodes, Nodes$degree > 1, select = c(Vertex, IndividualID))

# Use adjacent_vertices to list the individuals that have a mate in common (i.e. vertices adjacent to the vertices identified as having degree greater than 1) 
  IndividualsWithSharedMates1 <- adjacent_vertices(Network, IndividualsWithMultipleMates$Vertex, mode = c("all"))

# Turn that list into a dataframe, then convert the vertex numbers into the individual IDs:
  # Turn the list into a dataframe called IndividualsWithSharedMates2:
    IndividualsWithSharedMates2 <- plyr::ldply(IndividualsWithSharedMates1, rbind)
  # Now convert the vertex names to the actual individuals' names. 
   IndividualsWithSharedMates3 <- IndividualsWithSharedMates2
  IndividualsWithSharedMates3 <- as.data.frame(sapply(IndividualsWithSharedMates2, plyr::mapvalues, from = Nodes$Vertex, to = Nodes$IndividualID))

#Need to generate all permutations of row elements:
  # Produce a list of dataframes with all possible combinations of columns (1/2, 1/3, 1/4, 2/3, 2/4, 3/4)
    IndividualsSharingAMate1 <- combn(IndividualsWithSharedMates3, 2, simplify=FALSE)
  # Rename the column names in those dataframes so that they are all the same, and will bind together properly:
    ColumnNames <- c("1", "2")
    IndividualsSharingAMate2 <- lapply(IndividualsSharingAMate1, setNames, ColumnNames)
  # Bind that list of dataframes into a single dataframe:
    IndividualsSharingAMate3 <- plyr::ldply(IndividualsSharingAMate2, rbind)
  # Remove rows with NAs, i.e. rows that aren't a dyad (would results from individuals who had an odd number of mates instead of an even number)
    IndividualsSharingAMate4 <- na.omit(IndividualsSharingAMate3)
  # Now, since we made combinations of columns with combn, not permutations, we need to make another dataframe with the same columns, but in order 2,1 instead of 1,2 (in case the individuals are in the opposite order in RelatednessValues)
    IndividualsSharingAMate5 <- data.frame(IndividualsSharingAMate4$"2", IndividualsSharingAMate4$"1")
  # Then, again, rename the columns so we can join up our two sets and get one dataframe with all the permutations:
    names(IndividualsSharingAMate5) <- ColumnNames
  # And join the two dataframes to get all possible permutations:
    IndividualsSharingAMate6 <- bind_rows(IndividualsSharingAMate4, IndividualsSharingAMate5)

# Turn that dataframe into a set of dyad names, so we can search it against the dyad names in RelatednessValues:
  IndividualsSharingAMate7 <- unite(IndividualsSharingAMate6, "dyads", c("1", "2"), remove = FALSE)

# Search RelatednessValues$dyads against IndividualsSharingAMate6$dyads so that we can subset out the relatedness values of individuals that shared a mate:
  RelatednessValues$SharedMate <- RelatednessValues$dyads %in% IndividualsSharingAMate7$dyads
  # Now subset out the rows where SharedMate = TRUE, so that we can calculate mean relatedness for individuals that shared a mate. 
    RelatednessOfSharedMates <- subset(RelatednessValues, SharedMate == TRUE, select = c("dyads", "ind1.id", "ind2.id", "wang"))

# Calculate mean relatedness for those individuals that shared a mate:
  MeanOfSharedMates <- mean(RelatednessOfSharedMates$wang)

# Now create a bootstrapped distribution of mean relatedness for randomly sampled same-sex dyads, to test for departure from random mating:
  # Make a column in the output that tells us if a dyad is mixed sex or single sex (i.e. if the row contains the relatedness value for a pair of individuals that were potential mates, or the relatedness value for two males/two females.)
    RelatednessValues$MixedSex <- ifelse(grepl("M", RelatednessValues$dyads) & grepl("F", RelatednessValues$dyads), "YesMixedSex", "NoSingleSex")
  # Subset out only the same-sex dyads: 
    RelatednessOfPotentialSharers <- subset(RelatednessValues, MixedSex == "NoSingleSex", select = c("dyads", "ind1.id", "ind2.id", "wang"))
  # Bootstrapping: We need to randomly sample XXX pairs (number will depend on number of real mating-sharing dyads) from all possible single-sex dyads, do that a large number of times, and  calculate mean relatedness all those model populations, and turn that numeric vector into a dataframe, so you can make a distribution of means.  
    DistributionOfModelMeans <- replicate(n = 10000, expr = mean(sample(RelatednessOfPotentialSharers$wang, 43, replace = TRUE)))
    DistributionOfModelMeans <- as.data.frame(DistributionOfModelMeans)

# Now compare the real mean to the bootstrapped means- is it more than one standard deviation from the model population mean?
  Bootstrapping <- ggplot(data = DistributionOfModelMeans, mapping = aes(x = DistributionOfModelMeans)) + geom_histogram(mapping = aes(y = (..count../sum(..count..))), binwidth = 0.008) + geom_vline(xintercept = MeanOfSharedMates, size = 1, colour = "#FF3721", linetype = "dashed") + ggtitle("Null Distribution of Mean Relatedness in Picinguaba") + labs(x = "Mean Relatedness of Resample", y = "Proportion of Resamples") + theme_bw()
  plot(Bootstrapping)

