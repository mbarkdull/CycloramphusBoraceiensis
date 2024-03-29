# This code allows you to create a diagram of a mating system, represented as two rings of nodes (one for parents of one sex, and the second for parents of the second sex), connected by edge links that correspond to mating between the parents. I am using input data generated by parentage analysis conducted using the program Colony (Jones and Wang 2010).
# First, load the packages you'll need. 
  library(tidyverse)
  library(tidygraph)
  library(igraph)
  library(ggraph)

# Now, read in your node and edge data, and construct the graph object with ggraph. 
  Nodes <- read_csv("INSERTPARENTNODESHERE.csv")
  Edges <- read_csv("INSERTPARENTEDGESHERE.csv")
  Graph <- tbl_graph(nodes = Nodes, edges = Edges)

# Add a column for the degree of each indvidual in the graph (be sure to change the name of the output file when you run this code):
Degrees <- degree(Graph)
Nodes$degree <- Degrees
Nodes$Vertex <- V(g)
write_csv(Nodes, path = "INSERTPARENTNODESHERE.csv")

# Next, create an object that represents your graph with a bipartite layout. 
  g <- Graph
  V(g)$type <- bipartite_mapping(g)$type

# Based on that bipartite graph object, we will find x and y coordinates for each node that arrange the nodes into concentric circles. 
  # First, use create_layout to make an object, "coords", that describes the current, bipartite layout. coords contains columns for the x and y coordinates of each node in the linear bipartite graph.
    coords <- create_layout(g, layout = "bipartite") %>%
  # Next, we will use some trig to change the x and y coordinates from those of the linear bipartite graph to x and y coordinates that describe a circular graph. First, select the x and y columns. 
    select(x, y) %>%
  # This line creates a column, "theta", which basically divides the unit circle into slices based on the number of nodes we need to place along the edge of the circle, and assigns a slice to each node. 
    mutate(theta = x / (max(x) + 1) * 2 * pi, 
   # The next lines create r values, the radii of the circles the nodes will be placed on, and find the x and y coordinates on a Cartesian plane that will place each node along the edge of the circles. 
    r = y + 1, 
    x = r * cos(theta), 
    y = r * sin(theta))

# Now we create an entirely new graph layout, one that gives the node positions based on those circular coordinates we just created in coords. 
  my_graph <- create_layout(g, "manual", node.position = coords)

# I also want to make sure that the labels for each node are rotated in a logical way, so that there is no label overlapping and all labels are easy to read. 
  # To do this, first make an object called "label_data" based on "coords".
    label_data <- coords
  # Add a new column to label_data, called "angle", by converting the theta column from radians to degrees. 
    label_data$angle <- (label_data$theta)*(180/pi)
  # If we were to just use this angle value, then the labels in the third and fourth quadrants of the circle would be upside down. So we transform the angle values in these regions to flip them right side up. 
    label_data$plottingangle<-ifelse(label_data$angle < 270 & label_data$angle > 90, label_data$angle - 180, label_data$angle)

# I still need to figure out how to move the label text so that it isn't centered right on the dot for each node; my current plan is to do something similar to how I created label angle values, but then applying it to hjust. 
# > label_data$hjust <- ifelse(label_data$angle > 90 & label_data$angle < 180 & label_data$r > 1.5, 1, 2)

# Finally, we just plot the graph object that has the layout we create, specifying the angle of rotation for the labels in the aesthetic of geom_node_text. Tada!
  # Potentially might want to remove the filter() from geom_node_text(), depending on the context the figure will be used in. 
    ggraph(my_graph) + geom_edge_link0(edge_colour = "gray73") + geom_node_point(mapping = aes(color = LifeStage, size = degree, alpha = Inferred)) +  geom_node_text(aes(filter = degree > 1, label = IndividualID, angle = label_data$plottingangle), size = 3) + scale_shape_manual(values = c(0, 19)) +  coord_fixed() + theme_graph() + theme(legend.position = "none") + ggtitle("Mating System in Caraguatatuba")
    ParentNetworkGraph <- ggraph(my_graph) + geom_edge_link0(edge_colour = "gray73") + geom_node_point(mapping = aes(color = LifeStage, size = degree, alpha = Inferred)) +  geom_node_text(aes(filter = degree > 1, label = IndividualID, angle = label_data$plottingangle), size = 3) + scale_shape_manual(values = c(0, 19)) +  coord_fixed() + theme_graph() + theme(legend.position = "none") + ggtitle("Mating System in Caraguatatuba")

# Save your plot, if desired:
  ggsave("INSERTFILENAME.tiff")
    
# This graph changes the transparency of the edges depending on the probability value assigned to them by Colony:
  # I'd like to scale the transparency to match the bins of probability values we will mention in the text (>90%, >75%, etc.).
    ggraph(my_graph) + geom_edge_link0(mapping = aes(alpha = Probability), edge_colour = "gray73") + geom_node_point(mapping = aes(color = LifeStage, size = degree, alpha = Inferred)) +  geom_node_text(aes(filter = degree > 1, label = IndividualID, angle = label_data$plottingangle), size = 3) + scale_shape_manual(values = c(0, 19)) +  coord_fixed() + theme_graph() + theme(legend.position = "none") + ggtitle("Mating System in Caraguatatuba")
