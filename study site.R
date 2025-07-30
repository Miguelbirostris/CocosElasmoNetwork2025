library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(sf)
library(cowplot)

# Get world map data
world <- ne_countries(scale = "large", returnclass = "sf")

# Define the coordinates of Cocos Island
cocos_island <- data.frame(
  lon = -87.071667,
  lat = 5.518889
)


# Plot the map
ggplot(data = world) +
  geom_sf() +
  geom_point(data = cocos_island, aes(x = lon, y = lat), color = "red", size = 3) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-95, -80), ylim = c(-5, 15), expand = FALSE) +
  labs(title = "Study Site: Cocos Island, Costa Rica",
       x = "Longitude", y = "Latitude") +
  theme_minimal()




# Two ---------------------------------------------------------------------


# Plot of Costa Rica with Cocos Island marked
main_map <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = cocos_island, aes(x = lon, y = lat), color = "red", size = 3) +
  coord_sf(xlim = c(-90, -82), ylim = c(4, 12), expand = FALSE) +
  labs(title = "Study Site: Cocos Island, Costa Rica",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# Zoomed-in plot for Cocos Island
inset_map <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = cocos_island, aes(x = lon, y = lat), color = "red", size = 3) +
  coord_sf(xlim = c(-87, -85), ylim = c(4.5, 6.5), expand = T) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        plot.title = element_blank())
# Combine the main map and the inset
final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.65, y = 0.55, width = 0.3, height = 0.3)

# Display the final map
print(final_map)



# three -------------------------------------------------------------------


library(sf)
library(ggplot2)
library(rnaturalearth)
library(cowplot)

# Specify the path to the shapefile (excluding the file extension)
shapefile_path <- "CocosIslandShape.shp"

# Read the shapefile
cocos_shapefile <- st_read(shapefile_path)

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot the base map with Cocos Island polygon, no axis labels, only title
main_map <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  coord_sf(xlim = c(-90, -82), ylim = c(4, 12), expand = FALSE) +
  labs(title = " ",x = "Longitude", y = "Latitude") +
  theme_minimal()+ 
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  )



#+
  # theme(axis.title = element_blank(),
  #       axis.text = element_blank(),
  #       axis.ticks = element_blank())

# Create a zoomed-in inset for Cocos Island with a complete box and spaced coordinates
inset_map <- ggplot() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  coord_sf(xlim = c(-87.10, -87.02), ylim = c(5.49, 5.57), expand = FALSE) +
  scale_x_continuous(breaks = seq(-87.10, -87.02, by = 0.03)) +
  scale_y_continuous(breaks = seq(5.49, 5.57, by = 0.03)) +
 # labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


# Combine the main map and the inset
final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.25, y = 0.12, width = 0.5, height = 0.5) 

# Display the final map
print(final_map)

# Save the final map
ggsave("cocos_island_study_site_with_inset.jpg", plot = final_map, width = 10, height = 8)



# Four --------------------------------------------------------------------

library(sf)
library(ggplot2)
library(rnaturalearth)
library(cowplot)
library(ggspatial)

# Specify the path to the shapefile
shapefile_path <- "CocosIslandShape.shp"

# Read the shapefile
cocos_shapefile <- st_read(shapefile_path)

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Coordinates for Cocos Island (N 05°31'08", W87°04'18")
cocos_coords <- data.frame(
  lon = -87.0717,  # Converted W87°04'18" to decimal degrees
  lat = 5.5189     # Converted N 05°31'08" to decimal degrees
)

costa_rica_coords <- data.frame(
  lon = -83.8,     # Approximate longitude for central Costa Rica
  lat = 9.75       # Approximate latitude for central Costa Rica
)

# Plot the base map with Cocos Island polygon, no axis labels, only title
main_map <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  geom_point(data = cocos_coords, aes(x = lon, y = lat), color = "black", size = 3) +  # Point at Cocos location
  geom_text(data = costa_rica_coords, aes(x = lon, y = lat, label = "CR"), size = 5, fontface = "bold", color = "black") +  # Add "CR"
  coord_sf(xlim = c(-95, -82), ylim = c(4, 12), expand = FALSE) +  # Zoom out for smaller Costa Rica
  labs(title = " ", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(size=20),
        plot.margin = margin(10, 25, 10, 10))+
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.1, "in"), pad_y = unit(0.1, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  ) +
  ggspatial::annotation_scale(
    location = "bl",          # Place in the bottom left
    width_hint = 0.4,         # Adjust the width of the scale bar
    bar_cols = c("grey60", "white"),
    text_cex = 0.7
  )

# Create a zoomed-in inset for Cocos Island
inset_map <- ggplot() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  coord_sf(xlim = c(-87.10, -87.02), ylim = c(5.49, 5.57), expand = FALSE) +
  scale_x_continuous(breaks = seq(-87.10, -87.02, by = 0.03)) +
  scale_y_continuous(breaks = seq(5.49, 5.57, by = 0.03)) +
  theme_classic() +
  theme(
    text = element_text(size=20),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Firebrick border around the plot
    axis.line = element_blank(),     # Removes default axis lines (so they don't overlap)
    axis.ticks.length = unit(5, "pt"), # Adjust the length of the axis ticks
    axis.ticks = element_line(color = "black", size = 0.5) # Keeps the ticks, not the lines
  )



# Combine the main map and the inset, moving the inset to the upper left
final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.12, y = 0.36, width = 0.4, height = 0.5)  # Adjust x, y, width, height for positioning

# Display the final map
print(final_map)

ggsave("cocos_island_study_site_with_inset.jpg", plot = final_map, width = 10, height = 8)


# Five --------------------------------------------------------------------
library(sf)
library(ggplot2)
library(rnaturalearth)
library(cowplot)
library(ggspatial)

# Specify the path to the shapefile
shapefile_path <- "CocosIslandShape.shp"

# Read the shapefile
cocos_shapefile <- st_read(shapefile_path)

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Coordinates for Cocos Island (N 05°31'08", W87°04'18")
cocos_coords <- data.frame(
  lon = -87.0717,  # Converted W87°04'18" to decimal degrees
  lat = 5.5189     # Converted N 05°31'08" to decimal degrees
)

costa_rica_coords <- data.frame(
  lon = -84.8,     # Approximate longitude for central Costa Rica
  lat = 9.75       # Approximate latitude for central Costa Rica
)



# Updated main map with Pacific Ocean and Cocos label
main_map <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  geom_point(data = cocos_coords, aes(x = lon, y = lat), color = "black", size = 3) +  # Point at Cocos location
  geom_text(data = costa_rica_coords, aes(x = lon, y = lat, label = ""), size = 5, fontface = "bold", color = "black") +  # Add "CR"
  annotate("text", x = -86.8, y = 5.3, label = "Cocos", size = 4, fontface = "bold", color = "black") + # Label Cocos
  annotate("text", x = -84.5, y = 6.5, label = "Pacific Ocean",fontface = "bold", size = 7, color = "black") + # Label Pacific Ocean
  annotate("text", x = -82, y = 10.4, label = "Caribbean Sea", size = 5, fontface = "bold", color = "black") + # Label Caribbean Sea
    coord_sf(xlim = c(-93, -81), ylim = c(4, 12), expand = FALSE) +  # Zoom out for smaller Costa Rica
  labs(title = " ", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(size=20),
        plot.margin = margin(10, 25, 10, 10)) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  ) +
  ggspatial::annotation_scale(
    location = "bl",          # Place in the bottom left
    width_hint = 0.4,         # Adjust the width of the scale bar
    bar_cols = c("grey60", "white"),
    text_cex = 0.7
  )

# Updated inset map Cocos labels
inset_map <- ggplot() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  coord_sf(xlim = c(-87.10, -87.02), ylim = c(5.49, 5.57), expand = FALSE) +
  scale_x_continuous(breaks = seq(-87.10, -87.02, by = 0.03)) +
  scale_y_continuous(breaks = seq(5.49, 5.57, by = 0.03)) +
  annotate("text", x = -87.06, y = 5.53, label = "Cocos", size = 15, fontface = "bold", color = "black") + # Label Cocos
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    text = element_text(size=20),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Firebrick border around the plot
    axis.line = element_blank(),     # Removes default axis lines (so they don't overlap)
    axis.title = element_blank(),
    axis.ticks.length = unit(5, "pt"), # Adjust the length of the axis ticks
    axis.ticks = element_line(color = "black", size = 0.5) # Keeps the ticks, not the lines
  )

# Combine the main map and the inset, moving the inset to the upper left
final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.12, y = 0.36, width = 0.4, height = 0.5)  # Adjust x, y, width, height for positioning

# Display the final map
print(final_map)

# Save the updated map
ggsave("cocos_island_study_site_with_labels.jpg", plot = final_map, width = 10, height = 8)

# Save just main
ggsave("cocos_island_study_main map.png", plot = main_map, width = 10, height = 8)


# Save just main
ggsave("cocos_island_study_insetmap.png", plot = inset_map, width = 10, height = 8)

# Dive sites --------------------------------------------------------------------
library(sf)
library(ggplot2)
library(rnaturalearth)
library(cowplot)
library(ggspatial)

# Specify the path to the shapefile
shapefile_path <- "CocosIslandShape.shp"

# Read the shapefile
cocos_shapefile <- st_read(shapefile_path)

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Coordinates for Cocos Island (N 05°31'08", W87°04'18")
cocos_coords <- data.frame(
  lon = -87.0717,  # Converted W87°04'18" to decimal degrees
  lat = 5.5189     # Converted N 05°31'08" to decimal degrees
)

costa_rica_coords <- data.frame(
  lon = -84.8,     # Approximate longitude for central Costa Rica
  lat = 9.75       # Approximate latitude for central Costa Rica
)



# Updated main map with Pacific Ocean and Cocos label
main_map <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  geom_point(data = cocos_coords, aes(x = lon, y = lat), color = "black", size = 2) +  # Point at Cocos location
  geom_text(data = costa_rica_coords, aes(x = lon, y = lat, label = ""), size = 5, fontface = "bold", color = "black") +  # Add "CR"
  annotate("text", x = -86.8, y = 5.3, label = "Cocos", size = 4, fontface = "bold", color = "black") + # Label Cocos
  annotate("text", x = -84.5, y = 6.5, label = "Pacific Ocean",fontface = "bold", size = 7, color = "black") + # Label Pacific Ocean
  annotate("text", x = -82, y = 10.4, label = "Caribbean Sea", size = 5, fontface = "bold", color = "black") + # Label Caribbean Sea
  coord_sf(xlim = c(-93, -81), ylim = c(4, 12), expand = FALSE) +  # Zoom out for smaller Costa Rica
  labs(title = " ", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(size=20),
        plot.margin = margin(10, 25, 10, 10)) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.05, "in"), pad_y = unit(0.05, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  ) +
  ggspatial::annotation_scale(
    location = "bl",          # Place in the bottom left
    width_hint = 0.4,         # Adjust the width of the scale bar
    bar_cols = c("grey60", "white"),
    text_cex = 0.7
  )


#Create site coords

sites_coords<-data.frame(site=c("Alcyone", "AmigosG","AmigosP","Chatham", "DirtyRock", "ManuelitaA","Pajara", "PMaria","Sharkfin","Silverado","SRock","Ulloa","Viking"),
                         x=c(-87.03,-87.35,-87.10,-87.042,-87.08,-87.046,-87.055,-87.09,-87.08,-87.026,-87.05,-87.29,-87.065),
                         y=c(5.51,5.7,5.50,5.551,5.55,5.56,5.554,5.53,5.49,5.54,5.508,5.54,5.552))

#Calculate distance between sites
# Convert to sf object with geographic coordinates (WGS84)
sites_sf <- st_as_sf(sites_coords, coords = c("x", "y"), crs = 4326)

# Transform to projected CRS for accurate distance (UTM Zone 16N or 17N depending on location)
# Cocos Island is around 87°W, UTM Zone 16N should work well
sites_utm <- st_transform(sites_sf, crs = 32616)  # EPSG 32616 = WGS 84 / UTM zone 16N

# Calculate pairwise distances in meters
dist_matrix_m <- st_distance(sites_utm)

# Convert to kilometers and make a clean distance matrix
dist_matrix_km <- as.matrix(dist_matrix_m) / 1000
rownames(dist_matrix_km) <- sites_coords$site
colnames(dist_matrix_km) <- sites_coords$site

# Print distance matrix
RoundMatrix<-round(dist_matrix_km, 2)
print(RoundMatrix)

#Minimum
min(as.numeric(RoundMatrix)[as.numeric(RoundMatrix)>0])

#Maximum
max(as.numeric(RoundMatrix)[as.numeric(RoundMatrix)>0])

# Updated inset map Cocos labels
inset_map <- ggplot() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  coord_sf(xlim = c(-87.11, -87.02), ylim = c(5.48, 5.57), expand = FALSE) +
  scale_x_continuous(breaks = seq(-87.10, -87.02, by = 0.03)) +
  scale_y_continuous(breaks = seq(5.49, 5.57, by = 0.03)) +
  annotate("text", x = -87.06, y = 5.533, label = "Cocos", size = 13, fontface = "bold", color = "black") + # Label Cocos
  geom_point(data = sites_coords, aes(x = x, y = y), color = "black", size = 3)+
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    text = element_text(size=20),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Firebrick border around the plot
    axis.line = element_blank(),     # Removes default axis lines (so they don't overlap)
    axis.title = element_blank(),
    axis.ticks.length = unit(5, "pt"), # Adjust the length of the axis ticks
    axis.ticks = element_line(color = "black", size = 0.5) # Keeps the ticks, not the lines
  )+ #scale bar
  ggspatial::annotation_scale(
    location = "bl",          # Place in the bottom left
    width_hint = 0.4,         # Adjust the width of the scale bar
    bar_cols = c("grey30", "white"),
    text_cex = 0.7
  )
print(inset_map)
# Combine the main map and the inset, moving the inset to the upper left
final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.12, y = 0.36, width = 0.4, height = 0.5)  # Adjust x, y, width, height for positioning

# Display the final map
print(final_map)

# Save the updated map
ggsave("cocos_island_study_site_with_labels2.jpg", plot = final_map, width = 10, height = 8)

# Save just main
ggsave("cocos_island_study_main map.png", plot = main_map, width = 10, height = 8)

# Updated inset map Cocos labels with site annotations
inset_map_zoom <- ggplot() +
  geom_sf(data = cocos_shapefile, fill = "grey", color = "black") +
  coord_sf(xlim = c(-87.11, -87.02), ylim = c(5.48, 5.57), expand = FALSE) +
  scale_x_continuous(breaks = seq(-87.10, -87.02, by = 0.03)) +
  scale_y_continuous(breaks = seq(5.49, 5.57, by = 0.03)) +
  annotate("text", x = -87.06, y = 5.533, label = "Cocos", size = 13, fontface = "bold", color = "black") + # Label Cocos
  geom_point(data = sites_coords, aes(x = x, y = y), color = "black", size = 3) +  # Add points
  geom_text(data = sites_coords, aes(x = x, y = y, label = site), 
            hjust = 1, vjust = -0.5, size = 2, color = "black") +  # Add site labels
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    text = element_text(size=20),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Border around the plot
    axis.line = element_blank(),     # Removes default axis lines
    axis.title = element_blank(),
    axis.ticks.length = unit(5, "pt"), # Adjust the length of the axis ticks
    axis.ticks = element_line(color = "black", size = 0.5) # Keeps the ticks, not the lines
  ) + 
  ggspatial::annotation_scale(
    location = "bl",          # Place in the bottom left
    width_hint = 0.4,         # Adjust the width of the scale bar
    bar_cols = c("grey30", "white"),
    text_cex = 0.7
  )
print(inset_map_zoom)


# Save just inset
ggsave("cocos_island_study_insetmap.jpg", plot = inset_map_zoom, width = 10, height = 8)
