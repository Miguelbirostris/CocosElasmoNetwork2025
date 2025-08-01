#Code to reproduce a single network plot. Much of it is visualization and values relative to one another. see other scripts for more meaningful results (Like network strength)

#This code has not been optimized for speed

#Created by Miguel de Jesus Gomez Garcia
#Created: 01-August-2025
#Last edited: 01-August-2025
#First uploaded to GitHub on 01-August-2025


#Packages -----------------------------------------------------------------

library('igraph')
library('networkD3')
library(tidyverse)
library(raster)
library(Matrix)
library(ggplot2)
library(ggraph)
library(tidygraph)

## Prepare data ---------------------------------------------------------------------

#read formated data
Cocos_Edit<-read.csv("cocos_df.csv")[,-1]


#fix date format
Cocos_Edit$Date<-as.POSIXct(Cocos_Edit$Date,format="%Y-%m-%d")


sum(is.na(Cocos_Edit))

Cocos_Edit <- Cocos_Edit %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))%>%
  group_by(Date, Year, Month, Day) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = T)))

sum(is.na(Cocos_Edit))

# Arrange species alphabetically
species_columns <- sort(names(Cocos_Edit)[-(1:4)]) # Exclude Date, Year, Month, Day
Cocos_Edit <- Cocos_Edit %>%
  dplyr::select(Date, Year, Month, Day, all_of(species_columns))

# # Transform data to presence/absence (0 or 1)
#
# Cocos_Edit <- Cocos_Edit %>%
#   mutate(across(all_of(species_columns), ~ ifelse(.x > 0, 1, 0)))

head(Cocos_Edit)


Years<-as.character(unique(Cocos_Edit$Year))

#Create Seasonal grouping

Cocos_Season <- Cocos_Edit %>%
  mutate(Season = ifelse(Month %in% 6:11, "Wet", "Dry"))



# By year-season ----------------------------------------------------------

# Create a list of split data frames by Year and Season
Cocos_Season_Split <- split(Cocos_Season, list(Cocos_Season$Year, Cocos_Season$Season), drop = TRUE)

# Loop through and assign each one with a custom name
for (name in names(Cocos_Season_Split)) {
  # Clean the name: remove extra dots from split()
  clean_name <- gsub("\\.", "_", name)
  assign(paste0("Cocos_Season_", clean_name), Cocos_Season_Split[[name]])
}



# Example 2014 --------------------------------------------------------------------


##Single year, 2014. data grouped by day. Filtered to only wet season 2014 co-occurrences
data_year <- Cocos_Season_2014_Wet

# Melt data to long format to count species co-occurrences
data_long <- data_year %>%
  pivot_longer(cols = 5:17, names_to = "Species", values_to = "Abundance")

# Extract unique species
unique_species <- unique(data_long$Species)
unique_species <- species_columns

# Create pairs of species that co-occur
species_pairs <- expand.grid(unique_species, unique_species)
species_pairs <- species_pairs[species_pairs$Var1 != species_pairs$Var2, ]

# Count the number of months each pair of species co-occurs
co_occurrences <- species_pairs %>%
  rowwise() %>%
  mutate(count = sum(data_year[[as.character(Var1)]] > 0 & data_year[[as.character(Var2)]] > 0)) %>%
  filter(count > 0)

# Create a numeric ID for each species (node)
nodes <- data.frame(id = seq_along(unique_species) - 1, label = unique_species)

# Prepare edges data frame
edges <- co_occurrences %>%
  rename(from = Var1, to = Var2) %>%
  dplyr::select(from, to, count)

# Map species to numeric IDs in edges
edges <- edges %>%
  left_join(nodes, by = c("from" = "label")) %>%
  rename(from_id = id) %>%
  left_join(nodes, by = c("to" = "label")) %>%
  rename(to_id = id)

# Finalize edges and nodes
gedges <- edges %>%
  dplyr::select(from_id, to_id, count) %>%
  rename(from = from_id, to = to_id, value = count)

gedges$value<-sqrt(gedges$value)  # This improvis visualization a lot, by making edge lines thinner and easier to see

# Aggregate abundance by species
species_abundance <- data_long %>%
  group_by(Species) %>%
  summarize(Abundance = sum(Abundance, na.rm = TRUE))


g_nodes <- nodes %>%
  rename(name = label)

# Merge the abundance data with the node data
g_nodes <- g_nodes %>%
  left_join(species_abundance, by = c("name" = "Species"))

# Set default abundance to avoid NA values (optional, if any)
g_nodes$Abundance[is.na(g_nodes$Abundance)] <- 1  # or another value

# Normalize node sizes based on abundance (you can adjust this scaling factor)
g_nodes$size <- log(g_nodes$Abundance + 1) * 10  # Scaling and transforming abundance for visualization

g_nodes$size <- as.numeric(as.factor(g_nodes$Abundance)) # Scaling and transforming abundance for visualization #Factor scale

# Create the network graph
Species_net <- graph_from_data_frame(gedges, directed = TRUE, vertices = g_nodes)

# Set vertex color based on species groups for more distinction

V(Species_net)$col_values <- round( as.numeric(as.factor((V(Species_net)$name)), 2)) * 100
# Colour vertecies
colours <- colorRampPalette(c("gray50","tomato" ,"gold"))(4000)
V(Species_net)$color <- colours[ V(Species_net)$col_values ]
V(Species_net)$size<- g_nodes$size


# Use transparent edges for better visibility
E(Species_net)$width <- gedges$value
E(Species_net)$color <- adjustcolor("gray40", alpha.f = 0.5)



# tidygraph ---------------------------------------------------------------



# Convert to tidygraph format
Species_tidy_net <- as_tbl_graph(Species_net) # Generate the plot using ggraph with improvements

ggraph(Species_tidy_net, layout = "dh" ) +
  # Customize edge links with transparency and variable width
  geom_edge_link(aes(width = value),
                 alpha = 0.5,
                 color = "gray40",
                 show.legend = F) +
  # Add the white halo around nodes (slightly larger)
  geom_node_point(aes(size = g_nodes$size[match(name, g_nodes$name)] + 1), # Adjust this value for halo size
                  color = "white",
                  show.legend = FALSE) +
  # Customize node appearance (colored nodes)
  geom_node_point(aes(size = g_nodes$size[match(name, g_nodes$name)],
                      color = V(Species_net)$color),
                  show.legend = FALSE) +
  # Add node labels with repulsion to avoid overlapping
  geom_node_text(aes(label = name),
                 repel = TRUE, size = 25) +
  # Adjust scales for node size and color
  scale_size(range = c(5, 60)) +
  scale_color_identity() +

  # Set theme and title
  theme_void() +
  ggtitle(" ")

# Save the plot as a PNG image
ggsave(filename = "Network_Plot_recolor2014_wetc.jpg", width = 18, height = 15, units = "in", dpi = 300)




