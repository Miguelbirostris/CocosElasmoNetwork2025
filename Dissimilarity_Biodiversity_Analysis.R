#This script runs several biodiversity analysis fore a community of elasmobranchs in  Cocos Island, Costa Rica. 

#These analysis are described in Gomez-Garcia et al (in press) "Inter-specific relationships and their ecological role in an oceanic elasmobranch community" 

#This code has not been optimized for speed

#Created by Miguel de Jesus Gomez Garcia
#Created: 07-October-2024
#Last edited: 08-Apr-2025
#First uploaded to GitHub on 08-Apr-2025



# Libraries
library(tidyverse)
library(ggplot2)
library(viridis)
library(coda)
library(Metrics)
library(vegan)
library(dplyr)
library(vegan)
library(reshape2)
library(foreach)

# Load dataset
data_19 <- read.csv("data.csv")

# Limit printing for ease of visualization
options(max.print = 500)


# Anosim-Dissimilarity ----------------------------------------------------



## Yearly analysis ------------------------------------------------------------------

# Prepare species data for ANOSIM
species_data_Anosim <- data_19 %>%
  select(Year, YearWeek, Species, MeanCount) %>% #Group variables
  spread(key = Species, value = MeanCount, fill = 0) %>%  #replace missing entries with 0s. Can use pivot wider_too
  mutate(across(starts_with("Species"), sqrt, .names = "sqrt_{col}")) #sqrt transformation reduces the effect of extreme values 

# Remove YearWeek column and ensure numeric columns
species_matrix <- species_data_Anosim %>%
  select(-YearWeek, -Year) %>%
  select_if(is.numeric)

# Calculate Bray-Curtis dissimilarity matrix
dissimilarity_matrix <- vegdist(species_matrix, method = "jaccard")

# Perform ANOSIM analysis
anosim_result <- anosim(dissimilarity_matrix, grouping = species_data_Anosim$Year)

# Export ANOSIM summary and plot
print(summary(anosim_result), "AnosimSummary.csv")
#ggsave("AnosimPlot.jpg", plot(anosim_result))

# Extract the Bray-Curtis dissimilarity matrix
dissimilarity_matrix <- as.matrix(vegdist(species_matrix, method = "bray"))

# Create a dataframe to store the results
paired_dissimilarities <- data.frame(
  Year1 = integer(),
  Year2 = integer(),
  MeanDissimilarity = numeric()
)

# Compute mean dissimilarities between consecutive years through looping

unique_years <- sort(unique(species_data_Anosim$Year))
for (i in 1:(length(unique_years) - 1)) {
  year1 <- unique_years[i]
  year2 <- unique_years[i + 1]
  
  indices_year1 <- which(species_data_Anosim$Year == year1)
  indices_year2 <- which(species_data_Anosim$Year == year2)
  
  dissimilarities <- dissimilarity_matrix[indices_year1, indices_year2]
  mean_dissimilarity <- mean(dissimilarities)
  
  paired_dissimilarities <- rbind(paired_dissimilarities, data.frame(
    Year1 = year1,
    Year2 = year2,
    MeanDissimilarity = mean_dissimilarity
  ))
}

#Export results
write_csv(paired_dissimilarities, "PairedDissimilarities.csv")

# Plot mean dissimilarities between consecutive years

DisplotYear<- ggplot(paired_dissimilarities, aes(x = Year2, y = MeanDissimilarity)) +
  geom_line() +
  geom_point() +
  labs(title = "Mean Bray-Curtis Dissimilarity Between Consecutive Years",
       x = "Year",
       y = "Mean Dissimilarity") +
  theme_minimal()+
  scale_x_continuous(breaks=seq(1993,2020, by=3))

DisplotYear <- ggplot(paired_dissimilarities, aes(x = Year2, y = MeanDissimilarity)) +
  geom_line(color = "black", size = 1) +  # Thicker black line
  geom_point(size = 3, shape = 21, fill = "white", color = "firebrick", stroke = 1) +  # Points with white fill and firebrick border
  scale_x_continuous(name = "Year", breaks = c(seq(1993, 2018, by = 7), 2019)) +  # Adjust breaks to show years as per your example
  scale_y_continuous(limits = c(min(paired_dissimilarities$MeanDissimilarity) - 0.01, max(paired_dissimilarities$MeanDissimilarity) + 0.01), expand = c(0, 0)) +  # Adjust y-axis
  labs(title = "Mean Yearly Dissimilarity",
       x = "Year",
       y = "Mean Dissimilarity") +
  theme_classic(base_size = 14) +  # Increase base font size for readability
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),  # Bold and centered title
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkgray"),  # Centered subtitle (if needed)
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(color = "black")  # Black axis text
    # Optional: Adjust grid lines if needed
    # panel.grid.minor = element_blank(),  # Remove minor grid lines for a cleaner look
    # panel.grid.major = element_line(color = "gray80")  # Lighten major grid lines
  )

# Display Dissimilarity plot
print(DisplotYear)

# Save Dissimilarity plot
ggsave(filename = "MeanDissimilarityByYear2.jpg", plot = DisplotYear, width = 12, height = 8, dpi = 300)
  
# Perform SIMPER analysis

simper_result <- simper(species_matrix, species_data_Anosim$Year)
simper_summary <- summary(simper_result)

# Initialize an empty data frame to store the results
simper_top_species <- data.frame(Pair = character(), 
                                 Species = character(), 
                                 Contribution = numeric(), 
                                 stringsAsFactors = FALSE)

# Loop over each pairwise comparison in the simper summary
for (pair in names(simper_summary)) {
  # Extract the data for the current pair
  pair_data <- simper_summary[[pair]]
  
  # Sort by cumulative contribution (if not already sorted)
  pair_data <- pair_data[order(pair_data$cumsum), ]
  
  # Get the top 3 species based on the cumsum contribution
  top_species <- head(pair_data, 3)
  top_species$rank<-as.numeric(as.factor(top_species$cumsum))
  
  # Create a temporary data frame with the pair, species names, and cumsum values
  temp_df <- data.frame(Pair = rep(pair, 3), 
                        Rank =top_species$rank ,
                        Species = rownames(top_species), 
                        Contribution = 1-top_species$cumsum)
  
  # Append the results to the final data frame
  simper_top_species <- rbind(simper_top_species, temp_df)
}


#Check top contributing species.
unique(simper_top_species$Species)

#Details of top contributing species

summary(simper_top_species)
str(simper_top_species)
summary(simper_top_species[simper_top_species$Species=="Hammerheads",])
summary(simper_top_species[simper_top_species$Species=="Whitetips",])
summary(simper_top_species[simper_top_species$Species=="MarbledRays",])

# Count the number of instances each species ended in first, second, and third rank
rank_counts <- simper_top_species %>%
  filter(Rank %in% 1:3) %>%  # Filter to only include ranks 1, 2, and 3
  group_by(Species, Rank) %>%  # Group by Species and Rank
  summarise(Count = n()) %>%  # Count the number of instances
  arrange(Species, Rank)  # Arrange by Species and Rank

# View the result
print(rank_counts)

# Summarize each species
summary_hammerheads <- summary(simper_top_species[simper_top_species$Species == "Hammerheads",])
summary_whitetips <- summary(simper_top_species[simper_top_species$Species == "Whitetips",])
summary_marbledrays <- summary(simper_top_species[simper_top_species$Species == "MarbledRays",])

# Combine summaries into a single data frame
summary_combined <- rbind(
  data.frame(Species = "Hammerheads", summary_hammerheads),
  data.frame(Species = "Whitetips", summary_whitetips),
  data.frame(Species = "MarbledRays", summary_marbledrays)
)

# Save the summary to a CSV file
write.csv(summary_combined, "simper_top_species_summary.csv", row.names = FALSE)


# Reorder the Species factor within each Pair based on Contribution
simper_top_species <- simper_top_species %>%
  group_by(Pair) %>%
  arrange(desc(Contribution), .by_group = TRUE)%>%
  ungroup() 


#Add other species cathegory

simper_top_species_updated <- simper_top_species %>%
  group_by(Pair) %>%
  mutate(Contribution_Rank_1_to_3 = sum(Contribution[Rank %in% 1:3], na.rm = TRUE)) %>% # Calculate sum of contributions for Rank 1-3
  filter(Rank %in% 1:3) %>% # Keep only Rank 1-3 for further processing
  summarise(
    Rank = 4,
    Species = "Other",
    Contribution = 1 - Contribution_Rank_1_to_3
  ) %>%
  bind_rows(simper_top_species) %>% # Combine back with original data
  arrange(Pair, Rank)%>%
  unique() #Remove unnecesary duplicates

#fix decimal errors if needed
simper_top_species_updated$Contribution[simper_top_species_updated$Contribution< 0] <-0

# View the updated dataset
head(simper_top_species_updated)

#Paired dissimilarities plot
pair_plot<-  ggplot(simper_top_species_updated, aes(x = Pair, y = Contribution, fill = factor(Species, levels = unique(Species)))) +
  geom_bar(stat = "identity") +
  labs(title = "Top Species Contributions to Dissimilarity",
       x = "Comparison Pair",
       y = "Cumulative Contribution",
       fill = "Species") +
  scale_x_discrete(breaks = c("1993_1994","2001_2002", "2018_2019"),
                   labels = c("1993_1994","2001_2002", "2018_2019")) +
  scale_fill_viridis_d(option = "mako", direction = -1, begin=0, end = 0.85)+
  theme_classic() +
  theme(axis.text.x = element_text( vjust = 0.5, hjust=0.5))

#Print plot
print(pair_plot)

#Save plot
ggsave(filename = "TopPairPlot2.jpg", plot = pair_plot, width = 12, height = 8, dpi = 300)


## Monthly analysis ------------------------------------------------------------------

# Add Month column
data_19$Month <- format(as.Date(data_19$refdate), "%m")

# Prepare species data for ANOSIM by Month
species_data_Anosim_month <- data_19 %>%
  select(Month, YearWeek, Species, MeanCount) %>%
  spread(key = Species, value = MeanCount, fill = 0) %>%
  mutate(across(starts_with("Species"), sqrt, .names = "sqrt_{col}"))

# Remove YearWeek column and ensure numeric columns
species_matrix_month <- species_data_Anosim_month %>%
  select(-YearWeek, -Month) %>%
  select_if(is.numeric)

# Calculate Bray-Curtis dissimilarity matrix for months
dissimilarity_matrix_month <- as.matrix(vegdist(species_matrix_month, method = "jaccard"))

# Perform ANOSIM analysis by Month
anosim_result_month <- anosim(dissimilarity_matrix_month, grouping = species_data_Anosim_month$Month)

# Export ANOSIM summary and plot
#write_csv(summary(anosim_result_month), "AnosimMonthSummary.csv")
#ggsave("AnosimMonthPlot.jpg", plot(anosim_result_month))

# Compute mean dissimilarities between consecutive months


paired_dissimilarities_month <- data.frame(
  Month1 = character(),
  Month2 = character(),
  MeanDissimilarity = numeric(),
  stringsAsFactors = FALSE
)

unique_months <- sort(unique(species_data_Anosim_month$Month))

for (i in 1:(length(unique_months) - 1)) {
  month1 <- unique_months[i]
  month2 <- unique_months[i + 1]
  
  indices_month1 <- which(species_data_Anosim_month$Month == month1)
  indices_month2 <- which(species_data_Anosim_month$Month == month2)
  
  dissimilarities_month <- dissimilarity_matrix_month[indices_month1, indices_month2]
  mean_dissimilarity <- mean(dissimilarities_month)
  
  paired_dissimilarities_month <- rbind(paired_dissimilarities_month, data.frame(
    Month1 = month1,
    Month2 = month2,
    MeanDissimilarity = mean_dissimilarity
  ))
}

# Export ANOSIM summary and plot
print(paired_dissimilarities_month)


kruskal.test(MeanDissimilarity ~ Month1 , data = paired_dissimilarities_month)

# Plot mean dissimilarities between consecutive months
DisplotMonth <- ggplot(paired_dissimilarities_month, aes(x = as.numeric(Month2), y = MeanDissimilarity)) +
  geom_line(color = "black", size = 1) +  # Thicker line with a distinct color
  geom_point(size = 3, shape = 21, fill = "white", color = "firebrick", stroke = 1) +  # Add data points with a clean style
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(limits = c(min(paired_dissimilarities_month$MeanDissimilarity) - 0.025, max(paired_dissimilarities_month$MeanDissimilarity) + 0.025), expand = c(0, 0)) +  # Adjust y-axis
  labs(
    title = "Mean Bray-Curtis Dissimilarity Between Consecutive Months",
    x = "Month",
    y = "Mean Dissimilarity"
  ) +
  theme_classic(base_size = 25) +  # Increase base font size for readability
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkgray"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )
  
print(DisplotMonth)


ggsave("MeanDissimilarityByMonth2.jpg", plot=DisplotMonth,width = 12, height = 8, dpi = 300)


# Perform SIMPER analysis by Month
simper_result_month <- simper(species_matrix_month, species_data_Anosim_month$Month)
simper_summary_month <- summary(simper_result_month)
write_csv(as.data.frame(simper_summary_month$species), "SIMPER_Month_Species.csv")
write_csv(as.data.frame(simper_summary_month$groups), "SIMPER_Month_Groups.csv")


# Save the entire global environment to an .RData file
save.image("AnosimCocosPremise.RData")



# Biodiversity ------------------------------------------------------------

# Install and load the vegan package if not already installed
# install.packages("vegan")


# Prepare the data: pivot to get a species matrix (rows = unique YearWeek, columns = species)
species_matrix <- dcast(data_19, refdate ~ Species, value.var = "MeanCount", fun.aggregate = sum, fill = 0)

# Calculate the Shannon index for each unique YearWeek
shannon_index <- diversity(species_matrix[,-1], index = "shannon")  # Exclude the YearWeek column

simpson_index <- diversity(species_matrix[,-1], index = "simpson")  # Simpson index

# Create a data frame with YearWeek and corresponding Shannon Index
shannon_data <- data.frame(refdate = species_matrix$refdate, ShannonIndex = shannon_index)

# Check the result
head(shannon_data)


# Plotting the Shannon Index time series
ShannonPlot <- ggplot(shannon_data, 
                      aes(x = as.numeric(substr(refdate, 1, 4)) + (as.numeric(substr(refdate, 6, 7)) - 1) / 12, #x axis as continuous factorial year/month 
                          y = ShannonIndex, group = 1)) + #y axis, shanon index
  geom_line() +
  #geom_point() +
  labs(title = "Shannon Diversity Index Over Time (Weekly Data)",
       x = "Year",
       y = "Shannon Index") +
  theme_minimal()

# Display the plot
print(ShannonPlot)

# Save the plot
ggsave("ShannonIndex_TimeSeries.jpg", plot = ShannonPlot)



# Extract Year and Month from YearWeek

#By year
shannon_data$Year <- as.numeric(substr(shannon_data$refdate, 1, 4))


#By month
shannon_data$Month <- as.numeric(substr(shannon_data$refdate, 6, 7))

## Aggregate the Shannon Index by Month and year
shannon_monthly_Year <- aggregate(ShannonIndex ~ Year + Month, data = shannon_data, mean)


# Aggregate the Shannon Index by Year
shannon_monthly <- aggregate(ShannonIndex ~ Month, data = shannon_data, mean)


# Enhanced Shannon Index Plot by Month
ShannonPlotMonth <- ggplot(shannon_monthly, aes(x = Month, y = ShannonIndex)) +
  geom_line(color = "black", size = 1) +  # Thicker line with a distinct color
  geom_point(size = 3, shape = 21, fill = "white", color = "firebrick", stroke = 1) +  # Add data points with a clean style
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(limits = c(0.9, max(shannon_monthly$ShannonIndex) +0.025), expand = c(0, 0)) +  # Adjust y-axis
  labs(title = "Monthly Variation in Shannon Diversity Index",
       x = "Month",
       y = "Mean Shannon Index") +
  theme_classic(base_size = 14) +  # Increase base font size for readability
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkgray"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    #panel.grid.minor = element_blank(),  # Remove minor grid lines for a cleaner look
    #    panel.grid.major = element_line(color = "gray80")  # Lighten major grid lines
  )

# Display the plot
ShannonPlotMonth

# Save the plot
ggsave("ShannonIndexByMonth.jpg", plot = ShannonPlotMonth)

# Aggregate the Shannon Index by Year
shannon_yearly <- aggregate(ShannonIndex ~ Year, data = shannon_data, mean)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(ShannonIndex ~ Year, data = shannon_data)
kruskal_result

# Create a box plot for Shannon Index grouped by year

ShannonBoxPlotYear <- ggplot(shannon_monthly_Year, aes(x = factor(Year), y = ShannonIndex)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.color = "red", outlier.shape = 1) +
  scale_x_discrete(name = "Year", breaks = c(seq(min(shannon_monthly_Year$Year), max(shannon_monthly_Year$Year-2), by = 7),2019)) +
  labs(title = "Shannon Diversity Index by Year",
       y = "Shannon Index") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Display the plot
print(ShannonBoxPlotYear)

# Save the plot
ggsave("ShannonIndexByYear_BoxPlot_CustomBreaks.jpg", plot = ShannonBoxPlotYear)


# Enhanced Shannon Index Plot by Year
ShannonPlotYear <- ggplot(shannon_yearly, aes(x = Year, y = ShannonIndex)) +
  geom_line(color = "black", size = 1) +  # Thicker line with a distinct color
  geom_point(size = 3, shape = 21, fill = "white", color = "firebrick", stroke = 1) +  # Add data points with a clean style
  scale_x_continuous(name = "Year", breaks = c(seq(min(shannon_yearly$Year), max(shannon_yearly$Year-2), by = 7),2019)) +  # Ensure all years are shown
  scale_y_continuous(limits = c(0.75, max(shannon_yearly$ShannonIndex) + 0.025), expand = c(0, 0)) +  # Adjust y-axis
  labs(title = "Yearly Variation in Shannon Diversity Index",
       x = "Year",
       y = "Mean Shannon Index") +
  theme_classic(base_size = 14) +  # Increase base font size for readability
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkgray"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    # panel.grid.minor = element_blank(),  # Remove minor grid lines for a cleaner look
    # panel.grid.major = element_line(color = "gray80")  # Lighten major grid lines
  )

# Display the plot
ShannonPlotYear

# Save the plot
ggsave("ShannonIndexByYear_LinePlot.jpg", plot = ShannonPlotYear,width = 12, height = 8, dpi = 300)



## HillTransform -----------------------------------------------------------

#Hill transform
hill_index <- exp(shannon_index)  # Exclude the YearWeek column

# Create a data frame with YearWeek and corresponding hill Index
hill_data <- data.frame(refdate = species_matrix$refdate, hillIndex = hill_index)

# Check the result
head(hill_data)


# Plotting the hill Index time series
hillPlot <- ggplot(hill_data, 
                   aes(x = as.numeric(substr(refdate, 1, 4)) + (as.numeric(substr(refdate, 6, 7)) - 1) / 12, #x axis as continuous factorial year/month 
                       y = hillIndex, group = 1)) + #y axis, shanon index
  geom_line() +
  #geom_point() +
  labs(title = "Dill Diversity Index Over Time (Weekly Data)",
       x = "Year",
       y = "Hill Index") +
  theme_minimal()

hillPlot
# Extract Year and Month from YearWeek

# Save the plot
ggsave("HillIndex_TimeSeries.jpg", plot = hillPlot)

#By year
hill_data$Year <- as.numeric(substr(hill_data$refdate, 1, 4))


#By month
hill_data$Month <- as.numeric(substr(hill_data$refdate, 6, 7))

# Plotting the hill Index as a box plot grouped by month0

hillBoxPlot <- ggplot(hill_data, aes(x = factor(Month), y = hillIndex)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.color = "red", outlier.shape = 1) +
  scale_x_discrete(labels = month.abb) +
  labs(title = "hill Diversity Index by Month",
       x = "Month",
       y = "hill Index") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Display the plot
hillBoxPlot

# Aggregate the hill Index by Year
hill_yearly <- aggregate(hillIndex ~ Year, data = hill_data, mean)

# Create a box plot for hill Index grouped by year
hillPlotYear <- ggplot(hill_yearly, aes(x = Year, y = hillIndex)) +
  geom_line(color = "black", size = 1) +  # Thicker line with a distinct color
  geom_point(size = 3, shape = 21, fill = "white", color = "firebrick", stroke = 1) +  # Add data points with a clean style
  scale_x_continuous(name = "Year", breaks = c(seq(min(hill_yearly$Year), max(hill_yearly$Year-2), by = 7),2019)) +  # Ensure all years are shown
  scale_y_continuous(limits = c(2, max(hill_yearly$hillIndex) + 0.025), expand = c(0, 0)) +  # Adjust y-axis
  labs(title = "Yearly Variation in Hill Diversity",
       x = "Year",
       y = "Mean Hill Index") +
  theme_classic(base_size = 14) +  # Increase base font size for readability
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkgray"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    # panel.grid.minor = element_blank(),  # Remove minor grid lines for a cleaner look
    # panel.grid.major = element_line(color = "gray80")  # Lighten major grid lines
  )

# Display the plot
hillPlotYear

ggsave("HillByYear_LinePlot.jpg", plot = hillPlotYear,width = 12, height = 8, dpi = 300)

#data summary and tests

summary(shannon_yearly)
summary(shannon_monthly)
summary(shannon_data)


summary(lm(hillIndex ~ Year, data = hill_data))

kruskal.test(ShannonIndex ~ Month, data = shannon_data)


#By month
summary(lm(hillIndex ~ Month, data = hill_data))


# Aggregate the hill Index by Year
hill_monthly <- aggregate(hillIndex ~ Month, data = hill_data, mean)


# Enhanced hill Index Plot by Month
hillPlotMonth <- ggplot(hill_monthly, aes(x = Month, y = hillIndex)) +
  geom_line(color = "black", size = 1) +  # Thicker line with a distinct color
  geom_point(size = 3, shape = 21, fill = "white", color = "firebrick", stroke = 1) +  # Add data points with a clean style
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(limits = c(2.5, max(hill_monthly$hillIndex) +0.025), expand = c(0, 0)) +  # Adjust y-axis
  labs(title = "Monthly Variation in hill Diversity Index",
       x = "Month",
       y = "Mean hill Index") +
  theme_classic(base_size = 14) +  # Increase base font size for readability
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14, color = "darkgray"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    #panel.grid.minor = element_blank(),  # Remove minor grid lines for a cleaner look
    #    panel.grid.major = element_line(color = "gray80")  # Lighten major grid lines
  )

# Display the plot
hillPlotMonth
