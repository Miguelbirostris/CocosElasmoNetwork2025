#This script runs several network statistics for a community of elasmobranchs in  Cocos Island, Costa Rica.

#These analysis are described in Gomez-Garcia et al (in press) "Inter-specific relationships and their ecological role in an oceanic elasmobranch community"

#Note that some loops create several figures and tables at once in the work space, overwriting previous itterations. Run sequentially saving results outside the workspace as needed.

#This code has not been optimized for speed

#Created by Miguel de Jesus Gomez Garcia
#Created: 07-October-2024
#Last edited: 31-July-2025
#First uploaded to GitHub on 08-Apr-2025





#Network analysis -----------------------------------------------------------------

library('igraph')
library('networkD3')
library(tidyverse)
library(raster)
library(Matrix)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(mgcv)

## Prepare data ---------------------------------------------------------------------
#
# write.csv(Cocos_Edit,"cocos_df.csv")

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


Cocos_Edit_Year <- split(Cocos_Edit, Cocos_Edit$Year)
#BySeason ---------------------------------------------------------

#**PLEASE NOTE** that this will overwrite the previously created plots and tables.Save them elsewere before running the next part.


#Create Seasonal grouping

Cocos_Season <- Cocos_Edit %>%
  mutate(Season = ifelse(Month %in% 6:11, "Wet", "Dry"))


# Define the range of years to process
years <- 1993:2019



# Unique combinations of Year and Season
year_seasons <- Cocos_Season %>%
  ungroup() %>%                       # Remove any existing groupings
  dplyr::select(Year, Season) %>%     # Select only the Year and Season columns
  distinct()

# Create lists to store results for each year
co_occurrences_list <- list()
gedges_list <- list()
g_nodes_list <- list()
Species_net_list <- list()


# Initialize data frames to store results for each year
initialize_results <- function() {
  list(
    density_df = data.frame(Year = integer(), Density = numeric(), stringsAsFactors = FALSE),
    reciprocity_df = data.frame(Year = integer(), Reciprocity = numeric(), stringsAsFactors = FALSE),
    dyad_df = data.frame(Year = integer(), Mutual = integer(), Asymmetric = integer(), Null = integer(), stringsAsFactors = FALSE),
    calculated_reciprocity_df = data.frame(Year = integer(), Calculated_Reciprocity = numeric(), stringsAsFactors = FALSE),
    transitivity_df = data.frame(Year = integer(), Global_Transitivity = numeric(), Local_Transitivity = numeric(), stringsAsFactors = FALSE),
    triad_df = data.frame(Year = integer(), Triad_Census = I(list()), stringsAsFactors = FALSE),
    diameter_df = data.frame(Year = integer(), Diameter = numeric(), stringsAsFactors = FALSE),
    mean_distance_df = data.frame(Year = integer(), Mean_Distance_Undirected = numeric(), Mean_Distance_Directed = numeric(), stringsAsFactors = FALSE),
    degree_df = data.frame(Year = integer(), Degree = I(list()), stringsAsFactors = FALSE),
    strength_df = data.frame(Year = integer(), Strength = I(list()), stringsAsFactors = FALSE),
    degree_centrality_df = data.frame(Year = integer(), Degree_Centrality = I(list()), stringsAsFactors = FALSE),
    closeness_centrality_df = data.frame(Year = integer(), Closeness_Centrality = I(list()), stringsAsFactors = FALSE),
    eigenvector_centrality_df = data.frame(Year = integer(), Eigenvector_Centrality = I(list()), stringsAsFactors = FALSE),
    betweenness_centrality_df = data.frame(Year = integer(), Betweenness_Centrality = I(list()), stringsAsFactors = FALSE),
    edge_betweenness_df = data.frame(Year = integer(), Edge_Betweenness = I(list()), stringsAsFactors = FALSE),
    hub_scores_df = data.frame(Year = integer(), Hub_Scores = I(list()), stringsAsFactors = FALSE)
  )
}

export_results <- function(results) {
  for (name in names(results)) {
    file_name <- paste0(gsub("_df", "_Stats", name), ".txt")
    write.table(results[[name]], file = file_name, sep = "\t", row.names = FALSE)
  }
}


results <- initialize_results()


# Loop through each Year and Season
for (row in 1:nrow(year_seasons)) {
  year <- year_seasons$Year[row]
  season <- year_seasons$Season[row]

  # Filter data for the current year and season
  data_season <- Cocos_Season %>%
    filter(Year == year & Season == season) %>%
    summarise(across(Blacktips:Whitetips, ~ sum(.x, na.rm = TRUE)))

  # Melt data to long format to count species co-occurrences
  data_long <- data_season %>%
    pivot_longer(cols = 2:14, names_to = "Species", values_to = "Abundance")

  # Extract unique species
  unique_species <- unique(data_long$Species)
  unique_species <- species_columns

  # Create pairs of species that co-occur
  species_pairs <- expand.grid(unique_species, unique_species)
  species_pairs <- species_pairs[species_pairs$Var1 != species_pairs$Var2, ]

  # Count the total number of individuals co-occurring for each species pair
  co_occurrences <- species_pairs %>%
    rowwise() %>%
    mutate(count = sum(data_season[[as.character(Var1)]] > 0 & data_season[[as.character(Var2)]] > 0)) %>%
    filter(count > 0)

  # Prepare edges data frame with abundance-based weights
  edges <- co_occurrences %>%
    rename(from = Var1, to = Var2) %>%
    dplyr::select(from, to, count)

  # Assign weights based on the abundance counts
  edges$value <- edges$count

  # Create a numeric ID for each species (node)
  nodes <- data.frame(id = seq_along(unique_species) - 1, label = unique_species)

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

  gedges$value<-sqrt(gedges$value)  # This improves visualization a lot, by making edge lines thinner and easier to see

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



  # Create the network graph
  Species_net <- graph_from_data_frame(gedges, directed = TRUE, vertices = g_nodes)

  # Calculate the degree of each node
  node_degrees <- degree(Species_net, mode = "all")

  #Degree histogram plot
  png(filename = paste0("Degree_Histogram_", year,season, ".png"), width = 1000, height = 800)

  par(mar = c(7, 8, 4, 4) + 0.1,    # Margins: bottom, left, top, right
      mgp = c(5, 2, 0))           # mgp controls the margin for axis titles, labels, and tick marks


  hist(node_degrees,main="", xlab = "Node Degree", ylab= "Interaction Frequency",
       cex.axis = 4, cex.lab = 4)

  dev.off()

  #Node strength
  node_strength <- strength(Species_net, mode = "all", weights = E(Species_net)$value)

  #Strength histogram plot
  png(filename = paste0("Strength_Histogram_", year,season, ".png"), width = 1000, height = 800)

  par(mar = c(7, 8, 4, 4) + 0.1,    # Margins: bottom, left, top, right
      mgp = c(5, 2, 0))           # mgp controls the margin for axis titles, labels, and tick marks

  hist(node_strength,main="", xlab = "Node Strength", ylab= "Interaction Frequency",
       cex.axis = 4, cex.lab = 4)
  dev.off()

  # Store results in lists
  co_occurrences_list[[as.character(year)]] <- co_occurrences
  gedges_list[[as.character(year)]] <- gedges
  g_nodes_list[[as.character(year)]] <- g_nodes

  # Create the network graph
  Species_net <- graph_from_data_frame(gedges, directed = TRUE, vertices = g_nodes)
  Species_net_list[[as.character(year)]] <- Species_net


  g_nodes$size <- node_degrees
  #g_nodes$size <- as.numeric(as.factor(node_degrees[match(g_nodes$name, names(node_degrees))]))


  V(Species_net)$col_values <- round( as.numeric(as.factor((V(Species_net)$name)), 2)) * 100
  # Colour vertecies
  colours <- colorRampPalette(c("gray50","tomato" ,"gold"))(4000)
  V(Species_net)$color <- colours[ V(Species_net)$col_values ]
  V(Species_net)$size<- g_nodes$Abundance


  E(Species_net)$width <-gedges$value #In network plot, make edge line thickness represent the Strenght of the interaction

  # Convert to tidygraph format
  Species_tidy_net <- as_tbl_graph(Species_net) # Generate the plot using ggraph with improvements
  plot<-ggraph(Species_tidy_net, layout = "fr") +
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
                   repel = TRUE, size = 12) +
    # Adjust scales for node size and color
    scale_size(range = c(1, 20)) +
    scale_color_identity() +
    # Color identity will apply your predefined colors

    # Set theme and title
    theme_void() +
    ggtitle(paste("Species Co-occurrence Network (", year,season, ")", sep = ""))




  ggsave(paste("Species_Network_", year,season, ".jpg", sep = ""), plot, width = 10, height = 8)

  # Calculate network statistics
  density <- edge_density(Species_net, loops = FALSE)
  reciprocity_value <- reciprocity(Species_net)
  dyad <- dyad_census(Species_net)
  calculated_reciprocity <- 2 * dyad$mut / ecount(Species_net)
  global_transitivity <- transitivity(Species_net, type = "global")
  local_transitivity <- transitivity(Species_net, type = "local")
  triad <- triad_census(Species_net)
  diameter_value <- diameter(Species_net, directed = FALSE, weights = NA)
  mean_dist_undirected <- mean_distance(Species_net, directed = FALSE)
  mean_dist_directed <- mean_distance(Species_net, directed = TRUE)
  degree_value <- node_degrees
  strength_value <- node_strength
  degree_centrality <- degree(Species_net)
  closeness_centrality <- closeness(Species_net, mode = "all", weights = NA)
  eigen_centrality_value <- eigen_centrality(Species_net, directed = TRUE, weights = NA)$vector
  betweenness_centrality <- betweenness(Species_net, directed = FALSE, weights = NA)
  edge_betweenness_value <- edge_betweenness(Species_net, directed = FALSE, weights = NA)
  hs <- hub_score(Species_net, weights = NA)$vector
  as <- authority_score(Species_net, weights = NA)$vector

  # Append results to the respective data frames within the results list
  results$density_df <- rbind(results$density_df, data.frame(Year = year, Season = season, Density = density))
  results$reciprocity_df <- rbind(results$reciprocity_df, data.frame(Year = year, Season = season, Reciprocity = reciprocity_value))
  results$dyad_df <- rbind(results$dyad_df, data.frame(Year = year, Season = season, Mutual = dyad$mut, Asymmetric = dyad$asym, Null = dyad$null))
  results$calculated_reciprocity_df <- rbind(results$calculated_reciprocity_df, data.frame(Year = year, Season = season, Calculated_Reciprocity = calculated_reciprocity))
  results$transitivity_df <- rbind(results$transitivity_df, data.frame(Year = year, Season = season, Global_Transitivity = global_transitivity, Local_Transitivity = local_transitivity))
  results$triad_df <- rbind(results$triad_df, data.frame(Year = year, Season = season, Triad_Census = I(list(triad))))
  results$diameter_df <- rbind(results$diameter_df, data.frame(Year = year, Season = season, Diameter = diameter_value))
  results$mean_distance_df <- rbind(results$mean_distance_df, data.frame(Year = year, Season = season, Mean_Distance_Undirected = mean_dist_undirected, Mean_Distance_Directed = mean_dist_directed))
  results$degree_df <- rbind(results$degree_df, data.frame(Year = year, Season = season, Degree = I(list(degree_value))))
  results$strength_df <- rbind(results$strength_df, data.frame(Year = year, Season = season, Strength = I(list(strength_value))))
  results$degree_centrality_df <- rbind(results$degree_centrality_df, data.frame(Year = year, Season = season, Degree_Centrality = I(list(degree_centrality))))
  results$closeness_centrality_df <- rbind(results$closeness_centrality_df, data.frame(Year = year, Season = season, Closeness_Centrality = I(list(closeness_centrality))))
  results$eigenvector_centrality_df <- rbind(results$eigenvector_centrality_df, data.frame(Year = year, Season = season, Eigenvector_Centrality = I(list(eigen_centrality_value))))
  results$betweenness_centrality_df <- rbind(results$betweenness_centrality_df, data.frame(Year = year, Season = season, Betweenness_Centrality = I(list(betweenness_centrality))))
  results$edge_betweenness_df <- rbind(results$edge_betweenness_df, data.frame(Year = year, Season = season, Edge_Betweenness = I(list(edge_betweenness_value))))
  results$hub_scores_df <- rbind(results$hub_scores_df, data.frame(Year = year, Season = season, Hub_Scores = I(list(hs))))


}

results_season<-results


## Export Tables -----------------------------------------------------------

export_results(results_season)


## Degree / strength plots -------------------------------------------------------------



# Calculate the mean degree for each year
degree_summary <- results_season$degree_df %>%
  rowwise() %>%
  mutate(Mean_Degree = mean(unlist(Degree), na.rm = TRUE),
         SE_Degree = sd(unlist(Degree), na.rm = TRUE) / sqrt(length(unlist(Degree)))) %>%
  dplyr::select(Year,Season, Mean_Degree, SE_Degree)


# Plot Mean Degree over time with separate lines for dry and wet seasons
ggplot() +
  # Dry season line and points (solid line)
  geom_line(data = filter(degree_summary, Season == "Dry"),
            aes(x = Year, y = Mean_Degree, linetype = "Dry"), size = 1.5) +
  geom_point(data = filter(degree_summary, Season == "Dry"),
             aes(x = Year, y = Mean_Degree), color = "#7A0403FF", size = 3) +
  geom_errorbar(data = filter(degree_summary, Season == "Dry"),
                aes(x = Year, ymin = Mean_Degree - SE_Degree, ymax = Mean_Degree + SE_Degree),
                width = 0.2, color = "black") +

  # Wet season line and points (dashed line)
  geom_line(data = filter(degree_summary, Season == "Wet"),
            aes(x = Year, y = Mean_Degree, linetype = "Wet"), size = 1.5) +
  geom_point(data = filter(degree_summary, Season == "Wet"),
             aes(x = Year, y = Mean_Degree), color = "#1AE4B6FF", size = 3) +
  geom_errorbar(data = filter(degree_summary, Season == "Wet"),
                aes(x = Year, ymin = Mean_Degree - SE_Degree, ymax = Mean_Degree + SE_Degree),
                width = 0.2, color = "black") +

  labs(title = "Network Degree Over Time by Season",
       x = "Year",
       y = "Mean Degree") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"))




ggsave(filename = "Degree_Over_Season.png",
       plot = last_plot(),
       width = 10,
       height = 8,
       units = "in",
       dpi = 300)


# Calculate the mean and standard error of strength for each season
strength_summary_season <- results_season$strength_df %>%
  rowwise() %>%
  mutate(Mean_Strength = mean(unlist(Strength), na.rm = TRUE),
         SE_Strength = sd(unlist(Strength), na.rm = TRUE) / sqrt(length(unlist(Strength)))) %>%
  dplyr::select(Year, Season, Mean_Strength, SE_Strength)

# Plot the evolution of mean strength over time by season with error bars
ggplot() +
  # Dry season line and points (solid line)
  geom_line(data = filter(strength_summary_season, Season == "Dry"),
            aes(x = Year, y = Mean_Strength, linetype = "Dry"), size = 1.5) +
  geom_point(data = filter(strength_summary_season, Season == "Dry"),
             aes(x = Year, y = Mean_Strength), fill = "firebrick",color= "white", size = 3, shape = 21) +
  geom_errorbar(data = filter(strength_summary_season, Season == "Dry"),
                aes(x = Year, ymin = Mean_Strength - SE_Strength, ymax = Mean_Strength + SE_Strength),
                width = 0.2, color = "black") +

  # Wet season line and points (dashed line)
  geom_line(data = filter(strength_summary_season, Season == "Wet"),
            aes(x = Year, y = Mean_Strength, linetype = "Wet"), size = 1.5) +
  geom_point(data = filter(strength_summary_season, Season == "Wet"),
             aes(x = Year, y = Mean_Strength), fill = "firebrick", color= "white",size = 3, shape = 23) +
  geom_errorbar(data = filter(strength_summary_season, Season == "Wet"),
                aes(x = Year, ymin = Mean_Strength - SE_Strength, ymax = Mean_Strength + SE_Strength),
                width = 0.2, color = "black") +
  scale_y_continuous(limits = c(0, 160))+ #custo ylimits

  labs(title = "Network Strength Over Time by Season",
       x = "Year",
       y = "Mean Strength") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"))


ggsave(filename = "Strength_Over_Season.png",
       plot = last_plot(),
       width = 15,
       height = 8,
       units = "in",
       dpi = 300)


## Signifficance Models ----------------------------------------------------


#all year/seasons DF
strength_summary_season
strength_summary_season$Season<-as.factor(strength_summary_season$Season)

shapiro.test(strength_summary_season$Mean_Strength) #Normal

#LM
lm_all<-lm(Mean_Strength ~ Year*Season, data = strength_summary_season)
summary(lm_all)

# Boxplot for visualization
ggplot(strength_summary_season, aes(x = factor(Year), y = Mean_Strength)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Mean Network Strength by Year", x = "Year", y = "Mean Strength")+
  facet_wrap(~Season)


#Dry season DF
strength_summary_Dry<-strength_summary_season[strength_summary_season$Season=="Dry",]
shapiro.test(strength_summary_Dry$Mean_Strength)#Normal


lm_dry<-lm(Mean_Strength ~ Year, data = strength_summary_Dry)
summary(lm_dry)

gam_dry <- gam(Mean_Strength ~ s(Year), data = strength_summary_Dry)
summary(gam_dry)

# Plot the smooth trend
plot(gam_dry, shade = TRUE, main = "Dry Season",ylab="Fitted Strength")


#Wet season DF
strength_summary_wet<-strength_summary_season[strength_summary_season$Season=="Wet",]
shapiro.test(strength_summary_wet$Mean_Strength) #Normal


lm_wet<-lm(Mean_Strength ~ Year, data = strength_summary_wet)
summary(lm_wet)

gam_wet <- gam(Mean_Strength ~ s(Year), data = strength_summary_wet)
summary(gam_wet)

# Plot the smooth trend
plot(gam_wet, shade = TRUE, main = "Wet Season",ylab="Fitted Strength")


# Independent t-test for Season (Wet vs. Dry)
t_test_res <- t.test(Mean_Strength ~ Season, data = strength_summary_season)
print(t_test_res)

# Boxplot for season comparison
ggplot(strength_summary_season, aes(x = Season, y = Mean_Strength, fill = Season)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Mean Network Strength by Season", x = "Season", y = "Mean Strength")

#BySites + Season Long ---------------------------------------------------------



## ArrangeData -------------------------------------------------------------



Cocos_Edit <- read.csv("cocosseason_df.csv")

#fix date
Cocos_Edit$Date<-as.POSIXct(Cocos_Edit$Date,format="%Y-%m-%d")





sum(is.na(Cocos_Edit))

Cocos_Edit <- Cocos_Edit %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))%>%
  group_by(Date, Year, Month, Day, Site) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = T)))

sum(is.na(Cocos_Edit))

# Arrange species alphabetically
species_columns <- sort(names(Cocos_Edit)[-(1:5)]) # Exclude Date, Year, Month, Day
Cocos_Edit <- Cocos_Edit %>%
  dplyr::select(Date, Year, Month, Day, Site,all_of(species_columns))

# # Transform data to presence/absence (0 or 1)
#
# Cocos_Edit <- Cocos_Edit %>%
#   mutate(across(all_of(species_columns), ~ ifelse(.x > 0, 1, 0)))

head(Cocos_Edit)


Years<-as.character(unique(Cocos_Edit$Year))


Cocos_Edit_Year <- split(Cocos_Edit, Cocos_Edit$Year)


for (year in Years) {
  assign(paste0("Cocos_Edit_", year), Cocos_Edit_Year[[year]])
}



#**PLEASE NOTE** that this will overwrite the previously created plots and tables.Save them elsewere before running the next part.


#Create Seasonal grouping

Cocos_Season <- Cocos_Edit %>%
  mutate(Season = ifelse(Month %in% 6:11, "Wet", "Dry"))


# Define the range of years to process
years <- 1993:2019

#Merge close proximity (Manuelita) Sites

unique(Cocos_Season$Site)

Cocos_Season <- Cocos_Season %>%
  mutate(Site = if_else(str_detect(Site, "^Manuelita[A-F]$"), "Manuelita", Site))

unique(Cocos_Season$Site)

# Unique combinations of Year and Season and site
year_seasons_site <- Cocos_Season %>%
  ungroup() %>%                       # Remove any existing groupings
  dplyr::select(Year, Season,Site) %>%     # Select only the Year and Season columns
  distinct()

# Create lists to store results for each year
co_occurrences_list <- list()
gedges_list <- list()
g_nodes_list <- list()
Species_net_list <- list()

# Initialize data frames to store results for each year
# Initialize data frames to store results for each year
initialize_results <- function() {
  list(
    density_df = data.frame(Year = integer(), Density = numeric(), stringsAsFactors = FALSE),
    reciprocity_df = data.frame(Year = integer(), Reciprocity = numeric(), stringsAsFactors = FALSE),
    dyad_df = data.frame(Year = integer(), Mutual = integer(), Asymmetric = integer(), Null = integer(), stringsAsFactors = FALSE),
    calculated_reciprocity_df = data.frame(Year = integer(), Calculated_Reciprocity = numeric(), stringsAsFactors = FALSE),
    transitivity_df = data.frame(Year = integer(), Global_Transitivity = numeric(), Local_Transitivity = numeric(), stringsAsFactors = FALSE),
    triad_df = data.frame(Year = integer(), Triad_Census = I(list()), stringsAsFactors = FALSE),
    diameter_df = data.frame(Year = integer(), Diameter = numeric(), stringsAsFactors = FALSE),
    mean_distance_df = data.frame(Year = integer(), Mean_Distance_Undirected = numeric(), Mean_Distance_Directed = numeric(), stringsAsFactors = FALSE),
    degree_df = data.frame(Year = integer(), Degree = I(list()), stringsAsFactors = FALSE),
    strength_df = data.frame(Year = integer(), Strength = I(list()), stringsAsFactors = FALSE),
    degree_centrality_df = data.frame(Year = integer(), Degree_Centrality = I(list()), stringsAsFactors = FALSE),
    closeness_centrality_df = data.frame(Year = integer(), Closeness_Centrality = I(list()), stringsAsFactors = FALSE),
    eigenvector_centrality_df = data.frame(Year = integer(), Eigenvector_Centrality = I(list()), stringsAsFactors = FALSE),
    betweenness_centrality_df = data.frame(Year = integer(), Betweenness_Centrality = I(list()), stringsAsFactors = FALSE),
    edge_betweenness_df = data.frame(Year = integer(), Edge_Betweenness = I(list()), stringsAsFactors = FALSE),
    hub_scores_df = data.frame(Year = integer(), Hub_Scores = I(list()), stringsAsFactors = FALSE)
  )
}

export_results <- function(results) {
  for (name in names(results)) {
    file_name <- paste0(gsub("_df", "_Stats", name), ".txt")
    write.table(results[[name]], file = file_name, sep = "\t", row.names = FALSE)
  }
}


results <- initialize_results()




for (row in 1:nrow(year_seasons_site)) {
  year <- year_seasons_site$Year[row]
  season <- year_seasons_site$Season[row]
  site <- year_seasons_site$Site[row]

  # Filter data for the current combination
  data_season <- Cocos_Season %>%
    filter(Year == year & Season == season & Site == site) %>%
    summarise(across(Blacktips:Whitetips, ~ sum(.x, na.rm = TRUE)))

  # Reshape data into long format for species co-occurrence analysis
  data_long <- data_season %>%
    pivot_longer(cols = 5:17, names_to = "Species", values_to = "Abundance")


  # Extract unique species
  unique_species <- unique(data_long$Species)
  unique_species <- species_columns

  # Create pairs of species that co-occur
  species_pairs <- expand.grid(unique_species, unique_species)
  species_pairs <- species_pairs[species_pairs$Var1 != species_pairs$Var2, ]

  # Count the total number of individuals co-occurring for each species pair
  co_occurrences <- species_pairs %>%
    rowwise() %>%
    mutate(count = sum(data_season[[as.character(Var1)]] > 0 & data_season[[as.character(Var2)]] > 0)) %>%
    filter(count > 0)

  # Prepare edges data frame with abundance-based weights
  edges <- co_occurrences %>%
    rename(from = Var1, to = Var2) %>%
    dplyr::select(from, to, count)

  # Assign weights based on the abundance counts
  edges$value <- edges$count

  # Create a numeric ID for each species (node)
  nodes <- data.frame(id = seq_along(unique_species) - 1, label = unique_species)

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

  gedges$value<-sqrt(gedges$value)  # This improves visualization a lot, by making edge lines thinner and easier to see

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



  # Create the network graph
  Species_net <- graph_from_data_frame(gedges, directed = TRUE, vertices = g_nodes)

  # Calculate the degree of each node
  node_degrees <- degree(Species_net, mode = "all")

  #Degree histogram plot
  png(filename = paste0("Degree_Histogram_", year,season,site, ".png"), width = 1000, height = 800)

  par(mar = c(7, 8, 4, 4) + 0.1,    # Margins: bottom, left, top, right
      mgp = c(5, 2, 0))           # mgp controls the margin for axis titles, labels, and tick marks


  hist(node_degrees,main="", xlab = "Node Degree", ylab= "Interaction Frequency",
       cex.axis = 4, cex.lab = 4)

  dev.off()

  #Node strength
  node_strength <- strength(Species_net, mode = "all", weights = E(Species_net)$value)

  #Strength histogram plot
  png(filename = paste0("Strength_Histogram_", year,season,site, ".png"), width = 1000, height = 800)

  par(mar = c(7, 8, 4, 4) + 0.1,    # Margins: bottom, left, top, right
      mgp = c(5, 2, 0))           # mgp controls the margin for axis titles, labels, and tick marks

  hist(node_strength,main="", xlab = "Node Strength", ylab= "Interaction Frequency",
       cex.axis = 4, cex.lab = 4)
  dev.off()

  # Store results in lists
  co_occurrences_list[[as.character(year)]] <- co_occurrences
  gedges_list[[as.character(year)]] <- gedges
  g_nodes_list[[as.character(year)]] <- g_nodes

  # Create the network graph
  Species_net <- graph_from_data_frame(gedges, directed = TRUE, vertices = g_nodes)
  Species_net_list[[as.character(year)]] <- Species_net


  g_nodes$size <- node_degrees
  #g_nodes$size <- as.numeric(as.factor(node_degrees[match(g_nodes$name, names(node_degrees))]))


  V(Species_net)$col_values <- round( as.numeric(as.factor((V(Species_net)$name)), 2)) * 100
  # Colour vertecies
  colours <- colorRampPalette(c("gray50","tomato" ,"gold"))(4000)
  V(Species_net)$color <- colours[ V(Species_net)$col_values ]
  V(Species_net)$size<- g_nodes$Abundance


  E(Species_net)$width <-gedges$value #In network plot, make edge line thickness represent the Strenght of the interaction

  # Convert to tidygraph format
  Species_tidy_net <- as_tbl_graph(Species_net) # Generate the plot using ggraph with improvements
  plot<-ggraph(Species_tidy_net, layout = "fr") +
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
                   repel = TRUE, size = 12) +
    # Adjust scales for node size and color
    scale_size(range = c(1, 20)) +
    scale_color_identity() +
    # Color identity will apply your predefined colors

    # Set theme and title
    theme_void() +
    ggtitle(paste("Species Co-occurrence Network (", year,season, site,")", sep = ""))




  ggsave(paste("Species_Network_", year,season,site, ".jpg", sep = ""), plot, width = 10, height = 8)

  # Calculate network statistics
  density <- edge_density(Species_net, loops = FALSE)
  reciprocity_value <- reciprocity(Species_net)
  dyad <- dyad_census(Species_net)
  calculated_reciprocity <- 2 * dyad$mut / ecount(Species_net)
  global_transitivity <- transitivity(Species_net, type = "global")
  local_transitivity <- transitivity(Species_net, type = "local")
  triad <- triad_census(Species_net)
  diameter_value <- diameter(Species_net, directed = FALSE, weights = NA)
  mean_dist_undirected <- mean_distance(Species_net, directed = FALSE)
  mean_dist_directed <- mean_distance(Species_net, directed = TRUE)
  degree_value <- node_degrees
  strength_value <- node_strength
  degree_centrality <- degree(Species_net)
  closeness_centrality <- closeness(Species_net, mode = "all", weights = NA)
  eigen_centrality_value <- eigen_centrality(Species_net, directed = TRUE, weights = NA)$vector
  betweenness_centrality <- betweenness(Species_net, directed = FALSE, weights = NA)
  edge_betweenness_value <- edge_betweenness(Species_net, directed = FALSE, weights = NA)
  hs <- hub_score(Species_net, weights = NA)$vector
  as <- authority_score(Species_net, weights = NA)$vector

  # Append results to the respective data frames within the results list
  results$density_df <- rbind(results$density_df, data.frame(Year = year, Season = season, Site=site, Density = density))
  results$reciprocity_df <- rbind(results$reciprocity_df, data.frame(Year = year, Season = season, Site=site, Reciprocity = reciprocity_value))
  results$dyad_df <- rbind(results$dyad_df, data.frame(Year = year, Season = season, Site=site, Mutual = dyad$mut, Asymmetric = dyad$asym, Null = dyad$null))
  results$calculated_reciprocity_df <- rbind(results$calculated_reciprocity_df, data.frame(Year = year, Season = season, Site=site, Calculated_Reciprocity = calculated_reciprocity))
  results$transitivity_df <- rbind(results$transitivity_df, data.frame(Year = year, Season = season, Site=site, Global_Transitivity = global_transitivity, Local_Transitivity = local_transitivity))
  results$triad_df <- rbind(results$triad_df, data.frame(Year = year, Season = season, Site=site, Triad_Census = I(list(triad))))
  results$diameter_df <- rbind(results$diameter_df, data.frame(Year = year, Season = season, Site=site, Diameter = diameter_value))
  results$mean_distance_df <- rbind(results$mean_distance_df, data.frame(Year = year, Season = season, Site=site, Mean_Distance_Undirected = mean_dist_undirected, Mean_Distance_Directed = mean_dist_directed))
  results$degree_df <- rbind(results$degree_df, data.frame(Year = year, Season = season, Site=site, Degree = I(list(degree_value))))
  results$strength_df <- rbind(results$strength_df, data.frame(Year = year, Season = season, Site=site, Strength = I(list(strength_value))))
  results$degree_centrality_df <- rbind(results$degree_centrality_df, data.frame(Year = year, Season = season, Site=site, Degree_Centrality = I(list(degree_centrality))))
  results$closeness_centrality_df <- rbind(results$closeness_centrality_df, data.frame(Year = year, Season = season, Site=site, Closeness_Centrality = I(list(closeness_centrality))))
  results$eigenvector_centrality_df <- rbind(results$eigenvector_centrality_df, data.frame(Year = year, Season = season, Site=site, Eigenvector_Centrality = I(list(eigen_centrality_value))))
  results$betweenness_centrality_df <- rbind(results$betweenness_centrality_df, data.frame(Year = year, Season = season, Site=site, Betweenness_Centrality = I(list(betweenness_centrality))))
  results$edge_betweenness_df <- rbind(results$edge_betweenness_df, data.frame(Year = year, Season = season, Site=site, Edge_Betweenness = I(list(edge_betweenness_value))))
  results$hub_scores_df <- rbind(results$hub_scores_df, data.frame(Year = year, Season = season, Site=site, Hub_Scores = I(list(hs))))


}

results_season<-results


## Export Tables -----------------------------------------------------------

export_results(results_season)


# Calculate the mean degree for each year
degree_summary <- results_season$degree_df %>%
  rowwise() %>%
  mutate(Mean_Degree = mean(unlist(Degree), na.rm = TRUE),
         SE_Degree = sd(unlist(Degree), na.rm = TRUE) / sqrt(length(unlist(Degree)))) %>%
  dplyr::select(Year,Season,Site, Mean_Degree, SE_Degree)


# Plot Mean Degree over time with separate lines for dry and wet seasons
ggplot() +
  # Dry season line and points (solid line)
  geom_line(data = filter(degree_summary, Season == "Dry"),
            aes(x = Year, y = Mean_Degree, linetype = "Dry"), size = 1.5) +
  geom_point(data = filter(degree_summary, Season == "Dry"),
             aes(x = Year, y = Mean_Degree), color = "#7A0403FF", size = 3) +
  geom_errorbar(data = filter(degree_summary, Season == "Dry"),
                aes(x = Year, ymin = Mean_Degree - SE_Degree, ymax = Mean_Degree + SE_Degree),
                width = 0.2, color = "black") +

  # Wet season line and points (dashed line)
  geom_line(data = filter(degree_summary, Season == "Wet"),
            aes(x = Year, y = Mean_Degree, linetype = "Wet"), size = 1.5) +
  geom_point(data = filter(degree_summary, Season == "Wet"),
             aes(x = Year, y = Mean_Degree), color = "#1AE4B6FF", size = 3) +
  geom_errorbar(data = filter(degree_summary, Season == "Wet"),
                aes(x = Year, ymin = Mean_Degree - SE_Degree, ymax = Mean_Degree + SE_Degree),
                width = 0.2, color = "black") +

  labs(title = "Network Degree Over Time by Season",
       x = "Year",
       y = "Mean Degree") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"))+
  facet_wrap(~Site)




ggsave(filename = "Degree_Over_Season_site.png",
       plot = last_plot(),
       width = 30,
       height = 24,
       units = "in",
       dpi = 900)


# Calculate the mean and standard error of strength for each season
strength_summary_season <- results_season$strength_df %>%
  rowwise() %>%
  mutate(Mean_Strength = mean(unlist(Strength), na.rm = TRUE),
         SE_Strength = sd(unlist(Strength), na.rm = TRUE) / sqrt(length(unlist(Strength)))) %>%
  dplyr::select(Year, Season, Site, Mean_Strength, SE_Strength)

# Plot the evolution of mean strength over time by season with error bars
ggplot() +
  # Dry season line and points (solid line)
  geom_line(data = filter(strength_summary_season, Season == "Dry"),
            aes(x = Year, y = Mean_Strength, linetype = "Dry"), size = 1.5) +
  geom_point(data = filter(strength_summary_season, Season == "Dry"),
             aes(x = Year, y = Mean_Strength), fill = "firebrick",color= "white", size = 3, shape = 21) +
  geom_errorbar(data = filter(strength_summary_season, Season == "Dry"),
                aes(x = Year, ymin = Mean_Strength - SE_Strength, ymax = Mean_Strength + SE_Strength),
                width = 0.2, color = "black") +

  # Wet season line and points (dashed line)
  geom_line(data = filter(strength_summary_season, Season == "Wet"),
            aes(x = Year, y = Mean_Strength, linetype = "Wet"), size = 1.5) +
  geom_point(data = filter(strength_summary_season, Season == "Wet"),
             aes(x = Year, y = Mean_Strength), fill = "#1AE4B6FF", color= "white",size = 3, shape = 23) +
  geom_errorbar(data = filter(strength_summary_season, Season == "Wet"),
                aes(x = Year, ymin = Mean_Strength - SE_Strength, ymax = Mean_Strength + SE_Strength),
                width = 0.2, color = "black") +
  #scale_y_continuous(limits = c(0, 160))+ #custo ylimits

  labs(title = " ",
       x = "Year",
       y = "Mean Strength") +
  theme_classic(base_size = 40) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 80),
    axis.title.y = element_text(size = 80),
    axis.text = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black")
  ) +
  facet_wrap(~Site, scales = "free_y")


ggsave(filename = "Strength_Over_Season_Site.jpg",
       plot = last_plot(),
       width = 35,
       height = 25,
       units = "in",
       dpi = 600)
