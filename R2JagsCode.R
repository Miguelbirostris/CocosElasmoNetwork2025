#This script runs interactive-species hierarchical models for several  species of elasmobranchs in  Cocos Island, Costa Rica. 

#This model is desribed in Gomez-Garcia et al (in press) "Inter-specific relationships and their ecological role in an oceanic elasmobranch community" 

#Currently, the script uses a negative binomial distribution expressed as a Poisson model with diffuse priors for each covariate for abundant species, and a logistic binomial model following the same structure for rare species (see below and paper)

#This code has not been optimized for speed

#Created by Miguel de Jesus Gomez Garcia
#Created: 07-October-2024
#Last edited: 08-Apr-2025
#First uploaded to GitHub on 08-Apr-2025



#Load Libraries

library(tidyverse)
library(ggplot2)
library(viridis)
library(rjags)
library(R2jags)
library(coda)
library(Metrics)
library(zoo)


# Prepare data ------------------------------------------------------------

#Load Data
data_19 <- read.csv("data.csv")

#You may need this extra bit to remove the first column for some spreadsheets [,-1]

#format date
data_19$refdate <- as.Date(data_19$refdate, format="%Y-%m-%d")

#prevent huge console outputs
options(max.print=500)

#Data formating
data_19 <- data_19 %>%
  mutate(
    week = ifelse(week < 10, paste0("0", week), as.character(week)),
    Date = as.Date(refdate),
    YearWeek = paste(year(Date), week, sep = "-")
  ) %>%
  arrange(refdate)

# Get the unique species
unique_species <- unique(data_19$Species)

# Create a list to store data frames for each species
species_data_list <- list()

# Extract data for each species
for (species in unique_species) {
  species_data <- filter(data_19, Species == species)
  species_data_list[[species]] <- species_data
}

# Add counts of other species as covariates and replace NA with 0
for (species in unique_species) {
  for (other_species in unique_species) {
    if (species != other_species) {
      species_data_list[[species]] <- species_data_list[[species]] %>%
        left_join(
          select(species_data_list[[other_species]], YearWeek, MeanCount),
          by = "YearWeek",
          suffix = c("", paste0("_", other_species))
        ) %>%
        select(-matches(paste0("Species_", other_species)), -matches(paste0("week_", other_species)), -matches(paste0("Year_", other_species))) %>%
        rename_at(vars(paste0("MeanCount_", other_species)), ~paste0(other_species, "_Count")) %>%
        # Replace NA values with 0
        mutate(across(ends_with("_Count"), ~replace_na(., 0)))
    }
  }
}

# Assign each species-specific data frame to the global environment
for (species in unique_species) {
  assign(paste("data_", species, sep = ""), species_data_list[[species]])
}

# Print summary statistics for MeanCount by species
for (i in unique(data_19$Species)) {
  print(i)
  print(summary(data_19$MeanCount[data_19$Species == i]))
  print(sd(data_19$MeanCount[data_19$Species == i]))
}



## Blacktips ---------------------------------------------------------------

#Duplicate DF for binomial version
data_Blacktips_binom<-data_Blacktips

#Transform counts to presence-absence
data_Blacktips_binom$MeanCount[data_Blacktips_binom$MeanCount > 1] <- 1
data_Blacktips_binom$MeanCount[data_Blacktips_binom$MeanCount > 0] <- 1

# Create txt file for r2jags. You could do it inside R, but I like how organized and clear everything looks this way. Specially since we are doing it for multiple species 

sink("Logistic_Blacktips.txt")
cat("
model {

    # PRIORS
    Alpha ~ dnorm(0, 0.001)  # Intercept
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 1)
    Beta_SST ~ dnorm(0, 1)
    Beta_Year ~ dnorm(0, 1)  # New covariate coefficient for year
    Beta_EagleRays ~ dnorm(0, 1)  # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 1)  # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 1)  # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 1)  # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 1)  # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 1)  # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1)  # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 1)  # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 1)  # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 1)  # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 1)  # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 1)  # Covariate for Whitetips counts



    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
        
        # Main process
        Counts[i] ~ dbern(prob[i])
        logit(prob[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i]+
                          offset[i]

        offset[i] <- log(ndives[i])  # Define offset


        # Simulate new data
        Counts.new[i] ~ dbern(prob[i])

        # CHECK MODEL FIT
        presid[i] <- (Counts[i] - prob[i]) / sqrt(prob[i] * (1 - prob[i]))
        presid.new[i] <- (Counts.new[i] - prob[i]) / sqrt(prob[i] * (1 - prob[i]))

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

} # end of model"
    , fill = TRUE)
sink()


#Scale covariates

 <- data_Blacktips_binom %>%
  mutate(across(
    c(Visibility, SST, EagleRays_Count, Galapagos_Count, Hammerheads_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Silky_Count,
      Silvertips_Count, TigerSharks_Count, Turtles_Count, WhaleSharks_Count,
      Whitetips_Count),
    ~ scale(ceiling(.x))[, 1]
  ))


#Check data is ok
head(data_Blacktips_binom)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(Counts = ceiling(data_Blacktips_binom$MeanCount),
                 week = as.numeric(data_Blacktips_binom$week),
                 Year = as.numeric(as.factor(data_Blacktips_binom$Year)),
                 ndives = data_Blacktips_binom$ndives,
                 Visibility = data_Blacktips_binom$Visibility,
                 SST = data_Blacktips_binom$SST,
                 EagleRays_Count = data_Blacktips_binom$EagleRays_Count,
                 Galapagos_Count = data_Blacktips_binom$Galapagos_Count,
                 Hammerheads_Count = data_Blacktips_binom$Hammerheads_Count,
                 MantaRays_Count = data_Blacktips_binom$MantaRays_Count,
                 MarbledRays_Count = data_Blacktips_binom$MarbledRays_Count,
                 MobulaRays_Count = data_Blacktips_binom$MobulaRays_Count,
                 Silky_Count = data_Blacktips_binom$Silky_Count,
                 Silvertips_Count = data_Blacktips_binom$Silvertips_Count,
                 TigerSharks_Count = data_Blacktips_binom$TigerSharks_Count,
                 Turtles_Count = data_Blacktips_binom$Turtles_Count,
                 WhaleSharks_Count = data_Blacktips_binom$WhaleSharks_Count,
                 Whitetips_Count =data_Blacktips_binom$Whitetips_Count,
                 n = nrow(data_Blacktips)
)

# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_EagleRays", "Beta_Galapagos", "Beta_Hammerheads", "Beta_MantaRays",
            "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
            "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )


# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_Blacktips <- jags(win.data, inits, params, "Logistic_Blacktips.txt",
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      working.directory = getwd())
#print model
print(out_Blacktips)

o_nb_Blacktips <- out_Blacktips$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_Blacktips)

# Calculate summaries for each parameter
summary_nb_Blacktips <- lapply(params, function(param) {
  samples <- o_nb_Blacktips[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_Blacktips_df <- do.call(rbind, summary_nb_Blacktips)

write_csv(summary_nb_Blacktips_df, "BlacktipsModelSummaryResults.csv")

# Significant parameters
SummarySig_Blacktips <- summary_nb_Blacktips_df[with(summary_nb_Blacktips_df, Lower * Upper > 0), ]
rownames(SummarySig_Blacktips) <- 1:nrow(SummarySig_Blacktips)
Summaryfilter_Blacktips <- SummarySig_Blacktips[grep("^Beta_", SummarySig_Blacktips$Parameter), ]
rownames(Summaryfilter_Blacktips) <- 1:nrow(Summaryfilter_Blacktips)

write_csv(Summaryfilter_Blacktips, "BlacktipsModelSummarySignificant.csv")


# Pearson Residuals plot
res_samples <- o_nb_Blacktips$presid
res_samples_new <- o_nb_Blacktips$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals Blacktips Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_BlacktipsSeparate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals Blacktips Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Blacktips_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_Blacktips$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_Blacktips_binom)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_Blacktips <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts Blacktips") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_blacktipsSeparate_plot.jpg", plot = observed_predicted_plot_Blacktips, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_Blacktips <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted Blacktips Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_blacktipsSeparate.jpg", plot = join_difference_plot_Blacktips, width = 12, height = 8, dpi = 300)

## Galapagos ---------------------------------------------------------------


# Sink output to a separate file for Galapagos species
sink("Negbin_Galapagos.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2 
    sigma ~ dunif(0,20) 
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()


#Scale covariates
data_Galapagos <- data_Galapagos %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, EagleRays_Count, Hammerheads_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Silky_Count,
      Silvertips_Count, TigerSharks_Count, Turtles_Count, WhaleSharks_Count,
      Whitetips_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_Galapagos)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS

win.data <- list(Counts = ceiling(data_Galapagos$MeanCount),
                 week = as.numeric(data_Galapagos$week),
                 Year = as.numeric(as.factor(data_Galapagos$Year)),
                 ndives = data_Galapagos$ndives,
                 Visibility = data_Galapagos$Visibility,
                 SST = data_Galapagos$SST,
                 EagleRays_Count = data_Galapagos$EagleRays_Count,
                 Blacktips_Count = data_Galapagos$Blacktips_Count,
                 Hammerheads_Count = data_Galapagos$Hammerheads_Count,
                 MantaRays_Count = data_Galapagos$MantaRays_Count,
                 MarbledRays_Count = data_Galapagos$MarbledRays_Count,
                 MobulaRays_Count = data_Galapagos$MobulaRays_Count,
                 Silky_Count = data_Galapagos$Silky_Count,
                 Silvertips_Count = data_Galapagos$Silvertips_Count,
                 TigerSharks_Count = data_Galapagos$TigerSharks_Count,
                 Turtles_Count = data_Galapagos$Turtles_Count,
                 WhaleSharks_Count = data_Galapagos$WhaleSharks_Count,
                 Whitetips_Count = data_Galapagos$Whitetips_Count,
                 n = nrow(data_Galapagos))

# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_EagleRays", "Beta_Blacktips", "Beta_Hammerheads", "Beta_MantaRays",
            "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_Galapagos <- jags(win.data, inits, params, "Negbin_Galapagos.txt",
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      working.directory = getwd())

print(out_Galapagos)

o_nb_Galapagos <- out_Galapagos$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_Galapagos)

# Calculate summaries for each parameter
summary_nb_Galapagos <- lapply(params, function(param) {
  samples <- o_nb_Galapagos[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_Galapagos_df <- do.call(rbind, summary_nb_Galapagos)

write_csv(summary_nb_Galapagos_df, "GalapagosModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_Galapagos <- summary_nb_Galapagos_df[with(summary_nb_Galapagos_df, Lower * Upper > 0), ]
rownames(SummarySig_Galapagos) <- 1:nrow(SummarySig_Galapagos)
Summaryfilter_Galapagos <- SummarySig_Galapagos[grep("^Beta_", SummarySig_Galapagos$Parameter), ]
rownames(Summaryfilter_Galapagos) <- 1:nrow(Summaryfilter_Galapagos)

write_csv(Summaryfilter_Galapagos, "GalapagosModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_Galapagos$presid
res_samples_new <- o_nb_Galapagos$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals Galapagos Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Galapagos_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals Galapagos", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Galapagos_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_Galapagos$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_Galapagos)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_Galapagos <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts Galapagos Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_GalapagosSeparate_plot.jpg", plot = observed_predicted_plot_Galapagos, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_Galapagos <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted Galapagos Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_Galapagos.jpg", plot = join_difference_plot_Galapagos, width = 12, height = 8, dpi = 300)



## Hammerheads ---------------------------------------------------------------


# Sink output to a separate file for Hammerheads species
sink("Negbin_Hammerheads.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2    
    sigma ~ dunif(0,20)  
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    

    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()

#Scale covariates
data_Hammerheads <- data_Hammerheads %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, EagleRays_Count, Galapagos_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Silky_Count,
      Silvertips_Count, TigerSharks_Count, Turtles_Count, WhaleSharks_Count,
      Whitetips_Count),
    ~ scale(ceiling(.x))[, 1]
  ))
head(data_Hammerheads)
# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(Counts = ceiling(data_Hammerheads$MeanCount),
                 week = as.numeric(data_Hammerheads$week),
                 Year = as.numeric(as.factor(data_Hammerheads$Year)),
                 ndives = data_Hammerheads$ndives,
                 Visibility = data_Hammerheads$Visibility,
                 SST = data_Hammerheads$SST,
                 EagleRays_Count = data_Hammerheads$EagleRays_Count,
                 Blacktips_Count = data_Hammerheads$Blacktips_Count,
                 Galapagos_Count = data_Hammerheads$Galapagos_Count,
                 MantaRays_Count = data_Hammerheads$MantaRays_Count,
                 MarbledRays_Count = data_Hammerheads$MarbledRays_Count,
                 MobulaRays_Count = data_Hammerheads$MobulaRays_Count,
                 Silky_Count = data_Hammerheads$Silky_Count,
                 Silvertips_Count = data_Hammerheads$Silvertips_Count,
                 TigerSharks_Count = data_Hammerheads$TigerSharks_Count,
                 Turtles_Count = data_Hammerheads$Turtles_Count,
                 WhaleSharks_Count = data_Hammerheads$WhaleSharks_Count,
                 Whitetips_Count = data_Hammerheads$Whitetips_Count,
                 n = nrow(data_Hammerheads))


# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_EagleRays", "Beta_Blacktips", "Beta_Galapagos", "Beta_MantaRays",
            "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_Hammerheads <- jags(win.data, inits, params, "Negbin_Hammerheads.txt",
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      working.directory = getwd())
print(out_Hammerheads)


o_nb_Hammerheads <- out_Hammerheads$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_Hammerheads)

# Calculate summaries for each parameter
summary_nb_Hammerheads <- lapply(params, function(param) {
  samples <- o_nb_Hammerheads[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_Hammerheads_df <- do.call(rbind, summary_nb_Hammerheads)

write_csv(summary_nb_Hammerheads_df, "HammerheadsModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_Hammerheads <- summary_nb_Hammerheads_df[with(summary_nb_Hammerheads_df, Lower * Upper > 0), ]
rownames(SummarySig_Hammerheads) <- 1:nrow(SummarySig_Hammerheads)
Summaryfilter_Hammerheads <- SummarySig_Hammerheads[grep("^Beta_", SummarySig_Hammerheads$Parameter), ]
rownames(Summaryfilter_Hammerheads) <- 1:nrow(Summaryfilter_Hammerheads)

write_csv(Summaryfilter_Hammerheads, "HammerheadsModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_Hammerheads$presid
res_samples_new <- o_nb_Hammerheads$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals Hammerheads Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals Hammerheads Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_Hammerheads$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_Hammerheads)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_Hammerheads <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts Hammerheads Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_Hammerheads_plot_Separate.jpg", plot = observed_predicted_plot_Hammerheads, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_Hammerheads <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted Hammerheads Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_Hammerheads.jpg", plot = join_difference_plot_Hammerheads, width = 12, height = 8, dpi = 300)


## Silky ---------------------------------------------------------------

# Sink output to a separate file for Silky species
sink("Negbin_Silky.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2  
    sigma ~ dunif(0,20)  
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i]+
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()


#Scale covariates
data_Silky <- data_Silky %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, EagleRays_Count, Galapagos_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Hammerheads_Count,
      Silvertips_Count, TigerSharks_Count, Turtles_Count, WhaleSharks_Count,
      Whitetips_Count),
    ~ scale(ceiling(.x))[, 1]
  ))
head(data_Silky)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(Counts = ceiling(data_Silky$MeanCount),
                 week = as.numeric(data_Silky$week),
                 Year = as.numeric(as.factor(data_Silky$Year)),
                 ndives = data_Silky$ndives,
                 Visibility = data_Silky$Visibility,
                 SST = data_Silky$SST,
                 EagleRays_Count = data_Silky$EagleRays_Count,
                 Blacktips_Count = data_Silky$Blacktips_Count,
                 Galapagos_Count = data_Silky$Galapagos_Count,
                 Hammerheads_Count = data_Silky$Hammerheads_Count,
                 MantaRays_Count = data_Silky$MantaRays_Count,
                 MarbledRays_Count = data_Silky$MarbledRays_Count,
                 MobulaRays_Count = data_Silky$MobulaRays_Count,
                 Silvertips_Count = data_Silky$Silvertips_Count,
                 TigerSharks_Count = data_Silky$TigerSharks_Count,
                 Turtles_Count = data_Silky$Turtles_Count,
                 WhaleSharks_Count = data_Silky$WhaleSharks_Count,
                 Whitetips_Count = data_Silky$Whitetips_Count,
                 n = nrow(data_Silky))


# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_EagleRays", "Beta_Blacktips", "Beta_Galapagos", "Beta_Hammerheads",
            "Beta_MantaRays", "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_Silky <- jags(win.data, inits, params, "Negbin_Silky.txt",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                  working.directory = getwd())
print(out_Silky)

o_nb_Silky <- out_Silky$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_Silky)

# Calculate summaries for each parameter
summary_nb_Silky <- lapply(params, function(param) {
  samples <- o_nb_Silky[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_Silky_df <- do.call(rbind, summary_nb_Silky)

write_csv(summary_nb_Silky_df, "SilkyModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_Silky <- summary_nb_Silky_df[with(summary_nb_Silky_df, Lower * Upper > 0), ]
rownames(SummarySig_Silky) <- 1:nrow(SummarySig_Silky)
Summaryfilter_Silky <- SummarySig_Silky[grep("^Beta_", SummarySig_Silky$Parameter), ]
rownames(Summaryfilter_Silky) <- 1:nrow(Summaryfilter_Silky)

write_csv(Summaryfilter_Silky, "SilkyModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_Silky$presid
res_samples_new <- o_nb_Silky$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals Silky Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Silky_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals Silky Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_silky_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_Silky$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_Silky)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_Silky <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts Silky Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_Silky_plot_Separate.jpg", plot = observed_predicted_plot_Silky, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_Silky <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted Silky Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_Silky_Separate.jpg", plot = join_difference_plot_Silky, width = 12, height = 8, dpi = 300)

## Silvertips ---------------------------------------------------------------

# Sink output to a separate file for Silvertips species
sink("Negbin_Silvertips.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2  
    sigma ~ dunif(0,20)  
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts


    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()

#Scale covariates
data_Silvertips <- data_Silvertips %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, EagleRays_Count, Galapagos_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Hammerheads_Count,
      Silky_Count, TigerSharks_Count, Turtles_Count, WhaleSharks_Count,
      Whitetips_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_Silvertips)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(Counts = ceiling(data_Silvertips$MeanCount),
                 week = as.numeric(data_Silvertips$week),
                 Year = as.numeric(as.factor(data_Silvertips$Year)),
                 ndives = data_Silvertips$ndives,
                 Visibility = data_Silvertips$Visibility,
                 SST = data_Silvertips$SST,
                 EagleRays_Count = data_Silvertips$EagleRays_Count,
                 Blacktips_Count = data_Silvertips$Blacktips_Count,
                 Galapagos_Count = data_Silvertips$Galapagos_Count,
                 Hammerheads_Count = data_Silvertips$Hammerheads_Count,
                 MantaRays_Count = data_Silvertips$MantaRays_Count,
                 MarbledRays_Count = data_Silvertips$MarbledRays_Count,
                 MobulaRays_Count = data_Silvertips$MobulaRays_Count,
                 Silky_Count = data_Silvertips$Silky_Count,
                 TigerSharks_Count = data_Silvertips$TigerSharks_Count,
                 Turtles_Count = data_Silvertips$Turtles_Count,
                 WhaleSharks_Count = data_Silvertips$WhaleSharks_Count,
                 Whitetips_Count = data_Silvertips$Whitetips_Count,
                 n = nrow(data_Silvertips))

# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_EagleRays", "Beta_Blacktips", "Beta_Galapagos", "Beta_Hammerheads",
            "Beta_MantaRays", "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_Silvertips <- jags(win.data, inits, params, "Negbin_Silvertips.txt",
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       working.directory = getwd())

print(out_Silvertips)

o_nb_Silvertips <- out_Silvertips$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_Silvertips)

# Calculate summaries for each parameter
summary_nb_Silvertips <- lapply(params, function(param) {
  samples <- o_nb_Silvertips[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_Silvertips_df <- do.call(rbind, summary_nb_Silvertips)

write_csv(summary_nb_Silvertips_df, "SilvertipsModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_Silvertips <- summary_nb_Silvertips_df[with(summary_nb_Silvertips_df, Lower * Upper > 0), ]
rownames(SummarySig_Silvertips) <- 1:nrow(SummarySig_Silvertips)
Summaryfilter_Silvertips <- SummarySig_Silvertips[grep("^Beta_", SummarySig_Silvertips$Parameter), ]
rownames(Summaryfilter_Silvertips) <- 1:nrow(Summaryfilter_Silvertips)

write_csv(Summaryfilter_Silvertips, "SilvertipsModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_Silvertips$presid
res_samples_new <- o_nb_Silvertips$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals Silvertips Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Silvertips_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals Silvertips Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Silvertips_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_Silvertips$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_Silvertips)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_Silvertips <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts Silvertips Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_Silvertips_plot_Separate.jpg", plot = observed_predicted_plot_Silvertips, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_Silvertips <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted Silvertips Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_Silvertips_Separate.jpg", plot = join_difference_plot_Silvertips, width = 12, height = 8, dpi = 300)

## TigersSharks ---------------------------------------------------------------

# Sink output to a separate file for TigerSharks species
sink("Negbin_TigerSharks.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2   
    sigma ~ dunif(0,20) 
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()

#Scale covariates
data_TigerSharks <- data_TigerSharks %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, EagleRays_Count, Galapagos_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Hammerheads_Count,
      Silky_Count, Silvertips_Count, Turtles_Count, WhaleSharks_Count,
      Whitetips_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_TigerSharks)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(Counts = ceiling(data_TigerSharks$MeanCount),
                 week = as.numeric(data_TigerSharks$week),
                 Year = as.numeric(as.factor(data_TigerSharks$Year)),
                 ndives = data_TigerSharks$ndives,
                 Visibility = data_TigerSharks$Visibility,
                 SST = data_TigerSharks$SST,
                 EagleRays_Count = data_TigerSharks$EagleRays_Count,
                 Blacktips_Count = data_TigerSharks$Blacktips_Count,
                 Galapagos_Count = data_TigerSharks$Galapagos_Count,
                 Hammerheads_Count = data_TigerSharks$Hammerheads_Count,
                 MantaRays_Count = data_TigerSharks$MantaRays_Count,
                 MarbledRays_Count = data_TigerSharks$MarbledRays_Count,
                 MobulaRays_Count = data_TigerSharks$MobulaRays_Count,
                 Silky_Count = data_TigerSharks$Silky_Count,
                 Silvertips_Count = data_TigerSharks$Silvertips_Count,
                 Turtles_Count = data_TigerSharks$Turtles_Count,
                 WhaleSharks_Count = data_TigerSharks$WhaleSharks_Count,
                 Whitetips_Count = data_TigerSharks$Whitetips_Count,
                 n = nrow(data_TigerSharks))


# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_EagleRays", "Beta_Blacktips", "Beta_Galapagos", "Beta_Hammerheads",
            "Beta_MantaRays", "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky",
            "Beta_Silvertips", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_TigerSharks <- jags(win.data, inits, params, "Negbin_TigerSharks.txt",
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        working.directory = getwd())
print(out_TigerSharks)

o_nb_TigerSharks <- out_TigerSharks$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_TigerSharks)

# Calculate summaries for each parameter
summary_nb_TigerSharks <- lapply(params, function(param) {
  samples <- o_nb_TigerSharks[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_TigerSharks_df <- do.call(rbind, summary_nb_TigerSharks)

write_csv(summary_nb_TigerSharks_df, "TigerSharksModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_TigerSharks <- summary_nb_TigerSharks_df[with(summary_nb_TigerSharks_df, Lower * Upper > 0), ]
rownames(SummarySig_TigerSharks) <- 1:nrow(SummarySig_TigerSharks)
Summaryfilter_TigerSharks <- SummarySig_TigerSharks[grep("^Beta_", SummarySig_TigerSharks$Parameter), ]
rownames(Summaryfilter_TigerSharks) <- 1:nrow(Summaryfilter_TigerSharks)

write_csv(Summaryfilter_TigerSharks, "TigerSharksModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_TigerSharks$presid
res_samples_new <- o_nb_TigerSharks$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals TigerSharks Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals TigerSharks Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_TigerSharks$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_TigerSharks)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_TigerSharks <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts TigerSharks Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_TigerSharks_plot_Separate.jpg", plot = observed_predicted_plot_TigerSharks, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_TigerSharks <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted TigerSharks Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_TigerSharks_Separate.jpg", plot = join_difference_plot_TigerSharks, width = 12, height = 8, dpi = 300)


## Whitetips ---------------------------------------------------------------

# Sink output to a separate file for Whitetips species
sink("Negbin_Whitetips.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2   
    sigma ~ dunif(0,20) 
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()

#Scale covariates
data_Whitetips <- data_Whitetips %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, EagleRays_Count, Galapagos_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Hammerheads_Count,
      Silky_Count, Silvertips_Count, Turtles_Count, WhaleSharks_Count,
      TigerSharks_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_Whitetips)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(Counts = ceiling(data_Whitetips$MeanCount),
                 week = as.numeric(data_Whitetips$week),
                 Year = as.numeric(as.factor(data_Whitetips$Year)),
                 ndives = data_Whitetips$ndives,
                 Visibility = data_Whitetips$Visibility,
                 SST = data_Whitetips$SST,
                 EagleRays_Count = data_Whitetips$EagleRays_Count,
                 Blacktips_Count = data_Whitetips$Blacktips_Count,
                 Galapagos_Count = data_Whitetips$Galapagos_Count,
                 Hammerheads_Count = data_Whitetips$Hammerheads_Count,
                 MantaRays_Count = data_Whitetips$MantaRays_Count,
                 MarbledRays_Count = data_Whitetips$MarbledRays_Count,
                 MobulaRays_Count = data_Whitetips$MobulaRays_Count,
                 Silky_Count = data_Whitetips$Silky_Count,
                 Silvertips_Count = data_Whitetips$Silvertips_Count,
                 TigerSharks_Count = data_Whitetips$TigerSharks_Count,
                 Turtles_Count = data_Whitetips$Turtles_Count,
                 WhaleSharks_Count = data_Whitetips$WhaleSharks_Count,
                 n = nrow(data_Whitetips))


# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_EagleRays", "Beta_Blacktips", "Beta_Galapagos", "Beta_Hammerheads",
            "Beta_MantaRays", "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky",
            "Beta_Silvertips", "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_Whitetips <- jags(win.data, inits, params, "Negbin_Whitetips.txt",
                      n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                      working.directory = getwd())
print(out_Whitetips)

o_nb_Whitetips <- out_Whitetips$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_Whitetips)

# Calculate summaries for each parameter
summary_nb_Whitetips <- lapply(params, function(param) {
  samples <- o_nb_Whitetips[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_Whitetips_df <- do.call(rbind, summary_nb_Whitetips)

write_csv(summary_nb_Whitetips_df, "WhitetipsModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_Whitetips <- summary_nb_Whitetips_df[with(summary_nb_Whitetips_df, Lower * Upper > 0), ]
rownames(SummarySig_Whitetips) <- 1:nrow(SummarySig_Whitetips)
Summaryfilter_Whitetips <- SummarySig_Whitetips[grep("^Beta_", SummarySig_Whitetips$Parameter), ]
rownames(Summaryfilter_Whitetips) <- 1:nrow(Summaryfilter_Whitetips)

write_csv(Summaryfilter_Whitetips, "WhitetipsModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_Whitetips$presid
res_samples_new <- o_nb_Whitetips$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals Whitetips Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Whitetips_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals Whitetips Separate Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_Whitetips$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_Whitetips)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_Whitetips <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts Whitetips Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_Whitetips_plot_Separate.jpg", plot = observed_predicted_plot_Whitetips, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_Whitetips <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted Whitetips Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_Whitetips.jpg", plot = join_difference_plot_Whitetips, width = 12, height = 8, dpi = 300)

## WhaleSharks -------------------------------------------------------------

# Sink output to a separate file for WhaleSharks species
sink("Negbin_WhaleSharks.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2   
    sigma ~ dunif(0,20)  
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()

#Scale covariates
data_WhaleSharks <- data_WhaleSharks %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, EagleRays_Count, Galapagos_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Hammerheads_Count,
      Silky_Count, Silvertips_Count, Turtles_Count, Whitetips_Count,
      TigerSharks_Count),
    ~ scale(ceiling(.x))[, 1]
  ))
head(data_WhaleSharks)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(Counts = ceiling(data_WhaleSharks$MeanCount),
                 week = as.numeric(data_WhaleSharks$week),
                 Year = as.numeric(as.factor(data_WhaleSharks$Year)),
                 ndives = data_WhaleSharks$ndives,
                 Visibility = data_WhaleSharks$Visibility,
                 SST = data_WhaleSharks$SST,
                 EagleRays_Count = data_WhaleSharks$EagleRays_Count,
                 Blacktips_Count = data_WhaleSharks$Blacktips_Count,
                 Galapagos_Count = data_WhaleSharks$Galapagos_Count,
                 Hammerheads_Count = data_WhaleSharks$Hammerheads_Count,
                 MantaRays_Count = data_WhaleSharks$MantaRays_Count,
                 MarbledRays_Count = data_WhaleSharks$MarbledRays_Count,
                 MobulaRays_Count = data_WhaleSharks$MobulaRays_Count,
                 Silky_Count = data_WhaleSharks$Silky_Count,
                 Silvertips_Count = data_WhaleSharks$Silvertips_Count,
                 TigerSharks_Count = data_WhaleSharks$TigerSharks_Count,
                 Turtles_Count = data_WhaleSharks$Turtles_Count,
                 Whitetips_Count = data_WhaleSharks$Whitetips_Count,
                 n = nrow(data_WhaleSharks))


# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_EagleRays", "Beta_Blacktips", "Beta_Galapagos", "Beta_Hammerheads",
            "Beta_MantaRays", "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky",
            "Beta_Silvertips", "Beta_TigerSharks", "Beta_Turtles", "Beta_Whitetips",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_WhaleSharks <- jags(win.data, inits, params, "Negbin_WhaleSharks.txt",
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        working.directory = getwd())
print(out_WhaleSharks)

o_nb_WhaleSharks <- out_WhaleSharks$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_WhaleSharks)

# Calculate summaries for each parameter
summary_nb_WhaleSharks <- lapply(params, function(param) {
  samples <- o_nb_WhaleSharks[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_WhaleSharks_df <- do.call(rbind, summary_nb_WhaleSharks)

write_csv(summary_nb_WhaleSharks_df, "WhaleSharksModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_WhaleSharks <- summary_nb_WhaleSharks_df[with(summary_nb_WhaleSharks_df, Lower * Upper > 0), ]
rownames(SummarySig_WhaleSharks) <- 1:nrow(SummarySig_WhaleSharks)
Summaryfilter_WhaleSharks <- SummarySig_WhaleSharks[grep("^Beta_", SummarySig_WhaleSharks$Parameter), ]
rownames(Summaryfilter_WhaleSharks) <- 1:nrow(Summaryfilter_WhaleSharks)

write_csv(Summaryfilter_WhaleSharks, "WhaleSharksModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_WhaleSharks$presid
res_samples_new <- o_nb_WhaleSharks$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals WhaleSharks Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals WhaleSharks Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Whalesharks_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_WhaleSharks$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_WhaleSharks)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_WhaleSharks <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts WhaleSharks Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_WhaleSharks_plot_Separate.jpg", plot = observed_predicted_plot_WhaleSharks, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_WhaleSharks <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted WhaleSharks Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_WhaleSharks.jpg", plot = join_difference_plot_WhaleSharks, width = 12, height = 8, dpi = 300)



## EagleRays ---------------------------------------------------------------
data_EagleRays_binom<-data_EagleRays

# Make sure observations are transformed to Detection / non detection
data_EagleRays_binom$MeanCount[data_EagleRays_binom$MeanCount > 1] <- 1
data_EagleRays_binom$MeanCount[data_EagleRays_binom$MeanCount > 0] <- 1

#remove extreme values crashing this model
data_EagleRays_binom$Silky_Count[data_EagleRays_binom$Silky_Count>20]<-10

sink("Logistic_EagleRays.txt")
cat("
model {

    # PRIORS
    Alpha ~ dnorm(0, 0.001)  # Intercept
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 1)
    Beta_SST ~ dnorm(0, 1)
    Beta_Year ~ dnorm(0, 1)  # New covariate coefficient for year
    Beta_Blacktips ~ dnorm(0, 1)  # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 1)  # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 1)  # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 1)  # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 1)  # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 1)  # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 1)  # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 1)  # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 1)  # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 1)  # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 1)  # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 1)  # Covariate for Whitetips counts



    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
        
        # Main process
        Counts[i] ~ dbern(prob[i])
        logit(prob[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i]+
                          offset[i]

        offset[i] <- log(ndives[i])  # Define offset


        # Simulate new data
        Counts.new[i] ~ dbern(prob[i])

        # CHECK MODEL FIT
        presid[i] <- (Counts[i] - prob[i]) / sqrt(prob[i] * (1 - prob[i]))
        presid.new[i] <- (Counts.new[i] - prob[i]) / sqrt(prob[i] * (1 - prob[i]))

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

} # end of model"
    , fill = TRUE)
sink()

#Scale Covariates
data_EagleRays_binom <- data_EagleRays_binom %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, WhaleSharks_Count, Galapagos_Count,
      MantaRays_Count, MarbledRays_Count, MobulaRays_Count, Hammerheads_Count,
      Silky_Count, Silvertips_Count, Turtles_Count, Whitetips_Count,
      TigerSharks_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_EagleRays_binom)

# Set up initial values
inits <- function() list()


# Prepare data for JAGS
win.data <- list(
  Counts = ceiling(data_EagleRays_binom$MeanCount),
  week = as.numeric(data_EagleRays_binom$week),
  Year = as.numeric(as.factor(data_EagleRays_binom$Year)),
  ndives = data_EagleRays_binom$ndives,
  Visibility = data_EagleRays_binom$Visibility,
  SST = data_EagleRays_binom$SST,
  Blacktips_Count = data_EagleRays_binom$Blacktips_Count,
  Galapagos_Count = data_EagleRays_binom$Galapagos_Count,
  Hammerheads_Count = data_EagleRays_binom$Hammerheads_Count,
  MantaRays_Count = data_EagleRays_binom$MantaRays_Count,
  MarbledRays_Count = data_EagleRays_binom$MarbledRays_Count,
  MobulaRays_Count = data_EagleRays_binom$MobulaRays_Count,
  Silky_Count = data_EagleRays_binom$Silky_Count,
  Silvertips_Count = data_EagleRays_binom$Silvertips_Count,
  TigerSharks_Count = data_EagleRays_binom$TigerSharks_Count,
  Turtles_Count = data_EagleRays_binom$Turtles_Count,
  WhaleSharks_Count = data_EagleRays_binom$WhaleSharks_Count,
  Whitetips_Count = data_EagleRays_binom$Whitetips_Count,
  n = nrow(data_EagleRays_binom)
)

# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_Blacktips", "Beta_Galapagos", "Beta_Hammerheads", "Beta_MantaRays",
            "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
            "D.sum", "D.new.sum", "Bayes.P","Counts.new" ,"presid",
            "presid.new", "week_sin_cos")

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_EagleRays <- jags(win.data, inits, params, "Logistic_EagleRays.txt",
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          working.directory = getwd()) 

print(out_EagleRays)

o_nb_EagleRays <- out_EagleRays$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_EagleRays)

# Calculate summaries for each parameter
summary_nb_EagleRays <- lapply(params, function(param) {
  samples <- o_nb_EagleRays[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_EagleRays_df <- do.call(rbind, summary_nb_EagleRays)

write_csv(summary_nb_EagleRays_df, "EagleRaysModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_EagleRays <- summary_nb_EagleRays_df[with(summary_nb_EagleRays_df, Lower * Upper > 0), ]
rownames(SummarySig_EagleRays) <- 1:nrow(SummarySig_EagleRays)
Summaryfilter_EagleRays <- SummarySig_EagleRays[grep("^Beta_", SummarySig_EagleRays$Parameter), ]
rownames(Summaryfilter_EagleRays) <- 1:nrow(Summaryfilter_EagleRays)

write_csv(Summaryfilter_EagleRays, "EagleRaysModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_EagleRays$presid
res_samples_new <- o_nb_EagleRays$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals EagleRays Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Eaglerays_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals EagleRays Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_EagleRays_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_EagleRays$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_EagleRays_binom)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_EagleRays <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts EagleRays Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_EagleRays_plot_Separate.jpg", plot = observed_predicted_plot_EagleRays, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_EagleRays <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted EagleRays Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_EagleRays_Separate.jpg", plot = join_difference_plot_EagleRays, width = 12, height = 8, dpi = 300)


## MantaRays ----------------------------------------------------------------

data_MantaRays_binom<-data_MantaRays


data_MantaRays_binom$MeanCount[data_MantaRays_binom$MeanCount > 1] <- 1
data_MantaRays_binom$MeanCount[data_MantaRays_binom$MeanCount > 0] <- 1


sink("Logistic_MantaRays.txt")
cat("
model {

    # PRIORS
    Alpha ~ dnorm(0, 0.001)  # Intercept
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 1)
    Beta_SST ~ dnorm(0, 1)
    Beta_Year ~ dnorm(0, 1)  # New covariate coefficient for year
    Beta_Blacktips ~ dnorm(0, 1)  # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 1)  # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 1)  # Covariate for Hammerheads counts
    Beta_EagleRays ~ dnorm(0, 1)  # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 1)  # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 1)  # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 1)  # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 1)  # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 1)  # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 1)  # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 1)  # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 1)  # Covariate for Whitetips counts



    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
        
        # Main process
        Counts[i] ~ dbern(prob[i])
        logit(prob[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i]+
                          offset[i]

        offset[i] <- log(ndives[i])  # Define offset


        # Simulate new data
        Counts.new[i] ~ dbern(prob[i])

        # CHECK MODEL FIT
        presid[i] <- (Counts[i] - prob[i]) / sqrt(prob[i] * (1 - prob[i]))
        presid.new[i] <- (Counts.new[i] - prob[i]) / sqrt(prob[i] * (1 - prob[i]))

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

} # end of model"
    , fill = TRUE)
sink()

#Scale Covariates

data_MantaRays_binom <- data_MantaRays_binom %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, WhaleSharks_Count, Galapagos_Count,
      EagleRays_Count, MarbledRays_Count, MobulaRays_Count, Hammerheads_Count,
      Silky_Count, Silvertips_Count, Turtles_Count, Whitetips_Count,
      TigerSharks_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_MantaRays_binom)

# Set up initial values
inits <- function() list()


# Prepare data for JAGS

# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_Blacktips", "Beta_Galapagos", "Beta_Hammerheads", "Beta_EagleRays",
            "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
            "D.sum", "D.new.sum", "Bayes.P","Counts.new" ,"presid",
            "presid.new", "week_sin_cos")

win.data <- list(
  Counts = ceiling(data_MantaRays_binom$MeanCount),
  week = as.numeric(data_MantaRays_binom$week),
  Year = as.numeric(as.factor(data_MantaRays_binom$Year)),
  ndives = data_MantaRays_binom$ndives,
  Visibility = data_MantaRays_binom$Visibility,
  SST = data_MantaRays_binom$SST,
  EagleRays_Count = data_MantaRays_binom$EagleRays_Count,
  Galapagos_Count = data_MantaRays_binom$Galapagos_Count,
  Hammerheads_Count = data_MantaRays_binom$Hammerheads_Count,
  Blacktips_Count = data_MantaRays_binom$Blacktips_Count,
  MarbledRays_Count = data_MantaRays_binom$MarbledRays_Count,
  MobulaRays_Count = data_MantaRays_binom$MobulaRays_Count,
  Silky_Count = data_MantaRays_binom$Silky_Count,
  Silvertips_Count = data_MantaRays_binom$Silvertips_Count,
  TigerSharks_Count = data_MantaRays_binom$TigerSharks_Count,
  Turtles_Count = data_MantaRays_binom$Turtles_Count,
  WhaleSharks_Count = data_MantaRays_binom$WhaleSharks_Count,
  Whitetips_Count = data_MantaRays_binom$Whitetips_Count,
  n = nrow(data_MantaRays_binom)
)



# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_MantaRays <- jags(win.data, inits, params, "Logistic_MantaRays.txt",
                          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                          working.directory = getwd()) 
print(out_MantaRays)

o_nb_MantaRays <- out_MantaRays$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_MantaRays)

# Calculate summaries for each parameter
summary_nb_MantaRays <- lapply(params, function(param) {
  samples <- o_nb_MantaRays[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_MantaRays_df <- do.call(rbind, summary_nb_MantaRays)

write_csv(summary_nb_MantaRays_df, "MantaRaysModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_MantaRays <- summary_nb_MantaRays_df[with(summary_nb_MantaRays_df, Lower * Upper > 0), ]
rownames(SummarySig_MantaRays) <- 1:nrow(SummarySig_MantaRays)
Summaryfilter_MantaRays <- SummarySig_MantaRays[grep("^Beta_", SummarySig_MantaRays$Parameter), ]
rownames(Summaryfilter_MantaRays) <- 1:nrow(Summaryfilter_MantaRays)

write_csv(Summaryfilter_MantaRays, "MantaRaysModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_MantaRays$presid
res_samples_new <- o_nb_MantaRays$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals MantaRays Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_MantaRays_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals MantaRays Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_mantarays_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_MantaRays$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_MantaRays_binom)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_MantaRays <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts MantaRays Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_MantaRays_plot_Separate.jpg", plot = observed_predicted_plot_MantaRays, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_MantaRays <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted MantaRays Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_MantaRays_Separate.jpg", plot = join_difference_plot_MantaRays, width = 12, height = 8, dpi = 300)



## MarbledRays -------------------------------------------------------------

# Sink output to a separate file for MarbledRays species
sink("Negbin_MarbledRays.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2 
    sigma ~ dunif(0,20)   
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()

#scale Covariates

data_MarbledRays <- data_MarbledRays %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, WhaleSharks_Count, Galapagos_Count,
      EagleRays_Count, MantaRays_Count, MobulaRays_Count, Hammerheads_Count,
      Silky_Count, Silvertips_Count, Turtles_Count, Whitetips_Count,
      TigerSharks_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_MarbledRays)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(
  Counts = ceiling(data_MarbledRays$MeanCount),
  week = as.numeric(data_MarbledRays$week),
  Year = as.numeric(as.factor(data_MarbledRays$Year)),
  ndives = data_MarbledRays$ndives,
  Visibility = data_MarbledRays$Visibility,
  SST = data_MarbledRays$SST,
  Blacktips_Count = data_MarbledRays$Blacktips_Count,
  EagleRays_Count = data_MarbledRays$EagleRays_Count,
  Galapagos_Count = data_MarbledRays$Galapagos_Count,
  Hammerheads_Count = data_MarbledRays$Hammerheads_Count,
  MantaRays_Count = data_MarbledRays$MantaRays_Count,
  MobulaRays_Count = data_MarbledRays$MobulaRays_Count,
  Silky_Count = data_MarbledRays$Silky_Count,
  Silvertips_Count = data_MarbledRays$Silvertips_Count,
  TigerSharks_Count = data_MarbledRays$TigerSharks_Count,
  Turtles_Count = data_MarbledRays$Turtles_Count,
  WhaleSharks_Count = data_MarbledRays$WhaleSharks_Count,
  Whitetips_Count = data_MarbledRays$Whitetips_Count,
  n = nrow(data_MarbledRays)
)

# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_Blacktips", "Beta_EagleRays", "Beta_Galapagos", "Beta_Hammerheads",
            "Beta_MantaRays", "Beta_MobulaRays", "Beta_Silky", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_MarbledRays <- jags(win.data, inits, params, "Negbin_MarbledRays.txt",
                        n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                        working.directory = getwd())

print(out_MarbledRays)

o_nb_MarbledRays <- out_MarbledRays$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_MarbledRays)

# Calculate summaries for each parameter
summary_nb_MarbledRays <- lapply(params, function(param) {
  samples <- o_nb_MarbledRays[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_MarbledRays_df <- do.call(rbind, summary_nb_MarbledRays)

write_csv(summary_nb_MarbledRays_df, "MarbledRaysModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_MarbledRays <- summary_nb_MarbledRays_df[with(summary_nb_MarbledRays_df, Lower * Upper > 0), ]
rownames(SummarySig_MarbledRays) <- 1:nrow(SummarySig_MarbledRays)
Summaryfilter_MarbledRays <- SummarySig_MarbledRays[grep("^Beta_", SummarySig_MarbledRays$Parameter), ]
rownames(Summaryfilter_MarbledRays) <- 1:nrow(Summaryfilter_MarbledRays)

write_csv(Summaryfilter_MarbledRays, "MarbledRaysModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_MarbledRays$presid
res_samples_new <- o_nb_MarbledRays$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals MarbledRays Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Marbledrays_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals MarbledRays Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Marbledrays_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_MarbledRays$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_MarbledRays)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_MarbledRays <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts MarbledRays Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_MarbledRays_plot_Separate.jpg", plot = observed_predicted_plot_MarbledRays, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_MarbledRays <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted MarbledRays Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_MarbledRays_Separate.jpg", plot = join_difference_plot_MarbledRays, width = 12, height = 8, dpi = 300)


## MobulaRays --------------------------------------------------------------
# Sink output to a separate file for MobulaRays species
sink("Negbin_MobulaRays.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2    
    sigma ~ dunif(0,20)      
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_Turtles ~ dnorm(0, 0.1) # Covariate for Turtles counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_Turtles * Turtles_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()

#scale Covariates

data_MobulaRays <- data_MobulaRays %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, WhaleSharks_Count, Galapagos_Count,
      EagleRays_Count, MantaRays_Count, MarbledRays_Count, Hammerheads_Count,
      Silky_Count, Silvertips_Count, Turtles_Count, Whitetips_Count,
      TigerSharks_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_MobulaRays)

# Set up initial values
inits <- function() list()

# Prepare data for JAGS
win.data <- list(
  Counts = ceiling(data_MobulaRays$MeanCount),
  week = as.numeric(data_MobulaRays$week),
  Year = as.numeric(as.factor(data_MobulaRays$Year)),
  ndives = data_MobulaRays$ndives,
  Visibility = data_MobulaRays$Visibility,
  SST = data_MobulaRays$SST,
  Blacktips_Count = data_MobulaRays$Blacktips_Count,
  EagleRays_Count = data_MobulaRays$EagleRays_Count,
  Galapagos_Count = data_MobulaRays$Galapagos_Count,
  Hammerheads_Count = data_MobulaRays$Hammerheads_Count,
  MarbledRays_Count = data_MobulaRays$MarbledRays_Count,
  MantaRays_Count = data_MobulaRays$MantaRays_Count,
  Silky_Count = data_MobulaRays$Silky_Count,
  Silvertips_Count = data_MobulaRays$Silvertips_Count,
  TigerSharks_Count = data_MobulaRays$TigerSharks_Count,
  Turtles_Count = data_MobulaRays$Turtles_Count,
  WhaleSharks_Count = data_MobulaRays$WhaleSharks_Count,
  Whitetips_Count = data_MobulaRays$Whitetips_Count,
  n = nrow(data_MobulaRays)
)

# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_Blacktips", "Beta_EagleRays", "Beta_Galapagos", "Beta_Hammerheads",
            "Beta_MarbledRays", "Beta_MantaRays","Beta_Silky", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_Turtles", "Beta_WhaleSharks", "Beta_Whitetips",
             "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_MobulaRays <- jags(win.data, inits, params, "Negbin_MobulaRays.txt",
                       n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                       working.directory = getwd())

print(out_MobulaRays)

o_nb_MobulaRays <- out_MobulaRays$BUGSoutput$sims.list

# Extract main results
params <- names(o_nb_MobulaRays)

# Calculate summaries for each parameter
summary_nb_MobulaRays <- lapply(params, function(param) {
  samples <- o_nb_MobulaRays[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_MobulaRays_df <- do.call(rbind, summary_nb_MobulaRays)

write_csv(summary_nb_MobulaRays_df, "MobulaRaysModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_MobulaRays <- summary_nb_MobulaRays_df[with(summary_nb_MobulaRays_df, Lower * Upper > 0), ]
rownames(SummarySig_MobulaRays) <- 1:nrow(SummarySig_MobulaRays)
Summaryfilter_MobulaRays <- SummarySig_MobulaRays[grep("^Beta_", SummarySig_MobulaRays$Parameter), ]
rownames(Summaryfilter_MobulaRays) <- 1:nrow(Summaryfilter_MobulaRays)

write_csv(Summaryfilter_MobulaRays, "MobulaRaysModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_MobulaRays$presid
res_samples_new <- o_nb_MobulaRays$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals MobulaRays Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals MobulaRays Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_MobulaRays$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_MobulaRays)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_MobulaRays <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts MobulaRays Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_MobulaRays_plot_Separate.jpg", plot = observed_predicted_plot_MobulaRays, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_MobulaRays <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted MobulaRays Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_MobulaRays_Separate.jpg", plot = join_difference_plot_MobulaRays, width = 12, height = 8, dpi = 300)


## Turtles -------------------------------------------------------------

# Sink output to a separate file for Turtles species
sink("Negbin_Turtles.txt")
cat("
    model {

    # PRIORS
    tau <- 1/sigma^2   
    sigma ~ dunif(0,20)  
    Alpha ~ dnorm(0, 0.001) #intercept

    # Priors for covariates slope coefficients
    Beta_week ~ dunif(-2, 2)
    Beta_Visibility ~ dnorm(0, 0.1)
    Beta_SST ~ dnorm(0, 0.1)
    Beta_Year ~ dnorm(0, 0.1) # New covariate coefficient for year
    
    Beta_Blacktips ~ dnorm(0, 0.1) # Covariate for Blacktips counts
    Beta_EagleRays ~ dnorm(0, 0.1) # Covariate for EagleRays counts
    Beta_Galapagos ~ dnorm(0, 0.1) # Covariate for Galapagos counts
    Beta_Hammerheads ~ dnorm(0, 0.1) # Covariate for Hammerheads counts
    Beta_MantaRays ~ dnorm(0, 0.1) # Covariate for MantaRays counts
    Beta_MarbledRays ~ dnorm(0, 0.1) # Covariate for MarbledRays counts
    Beta_MobulaRays ~ dnorm(0, 0.1) # Covariate for MobulaRays counts
    Beta_Silky ~ dnorm(0, 0.1) # Covariate for Silky counts
    Beta_Silvertips ~ dnorm(0, 0.1) # Covariate for Silvertips counts
    Beta_TigerSharks ~ dnorm(0, 0.1) # Covariate for TigerSharks counts
    Beta_WhaleSharks ~ dnorm(0, 0.1) # Covariate for WhaleSharks counts
    Beta_Whitetips ~ dnorm(0, 0.1) # Covariate for Whitetips counts

     
     

    # Define cyclic transformation function for the week variable
    week_cyclic <- 2 * 3.1416 * week / 52  # Cyclic transformation for week

    # Define combined sine and cosine transformation
    week_sin_cos <- sin(week_cyclic) + cos(week_cyclic)

    # LIKELIHOOD
    for (i in 1:n) {
         Counts[i] ~ dpois(lambda[i])
        log(lambda[i]) <- Alpha + Beta_week * week_sin_cos[i] +
                          Beta_Visibility * Visibility[i] +
                          Beta_SST * SST[i] +
                          Beta_Year * Year[i] +
                          Beta_Blacktips * Blacktips_Count[i] +
                          Beta_EagleRays * EagleRays_Count[i] +
                          Beta_Galapagos * Galapagos_Count[i] +
                          Beta_Hammerheads * Hammerheads_Count[i] +
                          Beta_MantaRays * MantaRays_Count[i] +
                          Beta_MarbledRays * MarbledRays_Count[i] +
                          Beta_MobulaRays * MobulaRays_Count[i] +
                          Beta_Silky * Silky_Count[i] +
                          Beta_Silvertips * Silvertips_Count[i] +
                          Beta_TigerSharks * TigerSharks_Count[i] +
                          Beta_WhaleSharks * WhaleSharks_Count[i] +
                          Beta_Whitetips * Whitetips_Count[i] +
                          eps[i]+
                          offset[i]
                          
                          
        eps[i] ~ dnorm(0,tau)


        offset[i] <- log(ndives[i])  # Define offset

         

        # CHECK MODEL FIT
        presid[i]<- (Counts[i] - lambda[i])/sqrt(lambda[i])
        Counts.new[i] ~ dpois(lambda[i])
        presid.new[i] <- (Counts.new[i] - lambda[i])/sqrt(lambda[i])  

        D[i] <- presid[i]^2
        D.new[i] <- presid.new[i]^2
    }

    # CHECK MODEL FIT - BAYESIAN P-VALUE
    D.sum <- sum(D[])
    D.new.sum <- sum(D.new[])
    Bayes.P <- step(D.new.sum / D.sum - 1)

    } # end of model"
    , fill = TRUE)
sink()

# Set up initial values
inits <- function() list()


#scale Covariates

data_Turtles <- data_Turtles %>%
  mutate(across(
    c(Visibility, SST, Blacktips_Count, WhaleSharks_Count, Galapagos_Count,
      EagleRays_Count, MantaRays_Count, MarbledRays_Count, Hammerheads_Count,
      Silky_Count, Silvertips_Count, MobulaRays_Count, Whitetips_Count,
      TigerSharks_Count),
    ~ scale(ceiling(.x))[, 1]
  ))

head(data_Turtles)


# Prepare data for JAGS
win.data <- list(
  Counts = ceiling(data_Turtles$MeanCount),
  week = as.numeric(data_Turtles$week),
  Year = as.numeric(as.factor(data_Turtles$Year)),
  ndives = data_Turtles$ndives,
  Visibility = data_Turtles$Visibility,
  SST = data_Turtles$SST,
  Blacktips_Count = data_Turtles$Blacktips_Count,
  EagleRays_Count = data_Turtles$EagleRays_Count,
  Galapagos_Count = data_Turtles$Galapagos_Count,
  Hammerheads_Count = data_Turtles$Hammerheads_Count,
  MantaRays_Count = data_Turtles$MantaRays_Count,
  MarbledRays_Count = data_Turtles$MarbledRays_Count,
  MobulaRays_Count = data_Turtles$MobulaRays_Count,
  Silky_Count = data_Turtles$Silky_Count,
  Silvertips_Count = data_Turtles$Silvertips_Count,
  TigerSharks_Count = data_Turtles$TigerSharks_Count,
  WhaleSharks_Count = data_Turtles$WhaleSharks_Count,
  Whitetips_Count = data_Turtles$Whitetips_Count,
  n = nrow(data_Turtles)
)

# Parameters monitored
params <- c("Alpha", "Beta_week", "Beta_Visibility", "Beta_SST", "Beta_Year",
            "Beta_Blacktips", "Beta_EagleRays", "Beta_Galapagos", "Beta_Hammerheads","Beta_MantaRays",
            "Beta_MarbledRays", "Beta_MobulaRays", "Beta_Silky", "Beta_Silvertips",
            "Beta_TigerSharks", "Beta_WhaleSharks", "Beta_Whitetips",
            "D.sum", "D.new.sum", "Counts.new", "Bayes.P", "presid",
            "presid.new", "week_sin_cos"  )

# MCMC settings
ni <- 10000
nt <- 1
nb <- 1000
nc <- 3

# Run the model
out_Turtles <- jags(win.data, inits, params, "Negbin_Turtles.txt",
                    n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                    working.directory = getwd())
o_nb_Turtles <- out_Turtles$BUGSoutput$sims.list

print(out_Turtles)
# Extract main results
params <- names(o_nb_Turtles)

# Calculate summaries for each parameter
summary_nb_Turtles <- lapply(params, function(param) {
  samples <- o_nb_Turtles[[param]]
  mean_val <- mean(samples)
  sd_value <- sd(samples)
  hpdi <- HPDinterval(as.mcmc(samples), prob = 0.95)
  return(data.frame(Parameter = param, Mean = mean_val, SD = sd_value, Lower = hpdi[1], Upper = hpdi[2]))
})

# Combine into a data frame
summary_nb_Turtles_df <- do.call(rbind, summary_nb_Turtles)

write_csv(summary_nb_Turtles_df, "TurtlesModelSummaryResultsSeparate.csv")

# Significant parameters
SummarySig_Turtles <- summary_nb_Turtles_df[with(summary_nb_Turtles_df, Lower * Upper > 0), ]
rownames(SummarySig_Turtles) <- 1:nrow(SummarySig_Turtles)
Summaryfilter_Turtles <- SummarySig_Turtles[grep("^Beta_", SummarySig_Turtles$Parameter), ]
rownames(Summaryfilter_Turtles) <- 1:nrow(Summaryfilter_Turtles)

write_csv(Summaryfilter_Turtles, "TurtlesModelSummarySignificantSeparate.csv")


# Pearson Residuals plot
res_samples <- o_nb_Turtles$presid
res_samples_new <- o_nb_Turtles$presid.new

# Separate for observed and simulated residuals
res_means <- colMeans(res_samples)
res_means_05 <- apply(res_samples, 2, quantile, probs = 0.05)
res_means_95 <- apply(res_samples, 2, quantile, probs = 0.95)

res_new_means <- colMeans(res_samples_new)
res_new_means_05 <- apply(res_samples_new, 2, quantile, probs = 0.05)
res_new_means_95 <- apply(res_samples_new, 2, quantile, probs = 0.95)

# Create the histogram for observed residuals
observed_residuals_plot <- ggplot(data.frame(res_means), aes(x = res_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Observed Residuals Turtles Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("observed_residuals_histogram_Separate.jpg", plot = observed_residuals_plot, width = 8, height = 6, dpi = 300)

# Create the histogram for simulated residuals
simulated_residuals_plot <- ggplot(data.frame(res_new_means), aes(x = res_new_means)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Simulated Residuals Turtles Separate", x = "Residuals", y = "Frequency") +
  theme_minimal()

# Save the plot
ggsave("simulated_residuals_histogram_Separate.jpg", plot = simulated_residuals_plot, width = 8, height = 6, dpi = 300)

# Observed vs Predicted Data
sim_samples <- o_nb_Turtles$Counts.new
column_means <- colMeans(sim_samples)
percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)

plot_data <- data.frame(
  x = 1:ncol(sim_samples),
  mean = column_means,
  lower = percentile_05,
  upper = percentile_95
)

plot_data <- cbind(plot_data, data_Turtles)
plot_data <- plot_data %>%
  mutate(week = ifelse(week < 10, paste0("0", week), as.character(week))) %>%
  mutate(Date = as.Date(refdate)) %>%
  mutate(YearWeek = paste(year(Date), week, sep = "-")) %>%
  arrange(refdate) %>%
  mutate(MeanCount = ceiling(MeanCount))

# Create the plot for observed vs predicted counts
observed_predicted_plot_Turtles <- ggplot(plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean"), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count"), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts Turtles Separate") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

# Save the plot
ggsave("observed_predicted_Turtles_plot_Separate.jpg", plot = observed_predicted_plot_Turtles, width = 12, height = 8, dpi = 300)

# Difference plot
join_difference_plot_Turtles <- ggplot(plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Observed vs Predicted Turtles Separate") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.spacing.x = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10))

ggsave("difference_plot_Turtles_Separate.jpg", plot = join_difference_plot_Turtles, width = 12, height = 8, dpi = 300)

# Save the entire global environment to an .RData file
save.image("PoissonScaledPremise.RData")





# Plots and figures -------------------------------------------------------


# Define the model list with the appropriate reference data frames for each model
model_list <- list(
  "Blacktips" = list(model = o_nb_Blacktips, data = data_Blacktips_binom),
  "Galapagos" = list(model = o_nb_Galapagos, data = data_Galapagos),
  "Hammerheads" = list(model = o_nb_Hammerheads, data = data_Hammerheads),
  "Silky" = list(model = o_nb_Silky, data = data_Silky),
  "Silvertips" = list(model = o_nb_Silvertips, data = data_Silvertips),
  "TigerSharks" = list(model = o_nb_TigerSharks, data = data_TigerSharks),
  "Whitetips" = list(model = o_nb_Whitetips, data = data_Whitetips),
  "WhaleSharks" = list(model = o_nb_WhaleSharks, data = data_WhaleSharks),
  "EagleRays" = list(model = o_nb_EagleRays, data = data_EagleRays_binom),
  "MantaRays" = list(model = o_nb_MantaRays, data = data_MantaRays_binom),
  "MarbledRays" = list(model = o_nb_MarbledRays, data = data_MarbledRays),
  "MobulaRays" = list(model = o_nb_MobulaRays, data = data_MobulaRays),
  "Turtles" = list(model = o_nb_Turtles, data = data_Turtles)
)

# Initialize an empty list to store plot data for all species
all_plot_data <- list()

# Loop through each species to create plot data
for (species in names(model_list)) {
  
  # Extract model output and corresponding data
  model_output <- model_list[[species]][["model"]][["Counts.new"]]
  species_data <- model_list[[species]][["data"]]
  
  # Calculate column means and percentiles
  sim_samples <- model_output
  head(sim_samples)
  class(sim_samples)
  column_means <- colMeans(sim_samples)
  percentile_05 <- apply(sim_samples, 2, quantile, probs = 0.05)
  percentile_95 <- apply(sim_samples, 2, quantile, probs = 0.95)
  
  # Prepare plot data
  plot_data <- data.frame(
    refdate = species_data$refdate,
    MeanCount = species_data$MeanCount,
    mean = column_means,
    lower = percentile_05,
    upper = percentile_95,
    Species = species
  )
  
  # Add plot data to the list
  all_plot_data[[species]] <- plot_data
}

# Combine all plot data into a single data frame
combined_plot_data <- do.call(rbind, all_plot_data)

# Plot the observed vs predicted counts for all species
observed_predicted_join_plot <- ggplot(combined_plot_data, aes(x = refdate, y = mean, ymin = lower, ymax = upper, group = Species)) +  
  geom_point(aes(y = MeanCount), size = 3, shape = 16, alpha = 0.5) +
  labs(x = "Weeks", y = "Predicted counts") +
  geom_ribbon(alpha = 0.3, aes(fill = "95th percentile")) +
  geom_line(aes(y = mean, color = "Predicted Mean", group = Species), size = 1) +
  geom_line(aes(y = MeanCount, color = "Observed Count", group = Species), size = 1, linetype = "dashed") +
  labs(title = "Observed and Predicted Counts") +
  scale_color_manual(values = c("Observed Count" = "black", "Predicted Mean" = "firebrick")) +
  scale_fill_manual(values = "blue", guide = "legend", name = "Confidence Interval") +
  theme_minimal(base_size = 15) +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank(), 
        legend.spacing.x = unit(0.5, "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10)) +
  facet_wrap(~ Species, scales = "free_y", nrow = 3)

# Save the combined observed vs predicted plot
ggsave("Poiss_predicted_join_plot.jpg", plot = observed_predicted_join_plot, width = 12, height = 8, dpi = 300)

# Plot the difference (Predicted - Observed counts) for all species
join_difference_plot <- ggplot(combined_plot_data, aes(x = refdate, y = ceiling(MeanCount) - mean, group = Species)) +
  geom_errorbar(aes(ymin = ceiling(MeanCount) - mean - sd(ceiling(MeanCount) - mean), ymax = ceiling(MeanCount) - mean + sd(ceiling(MeanCount) - mean)), width = 0.2, alpha = 0.5) +
  geom_point(size = 3, shape = 16, alpha = 0.5, color = "firebrick") +
  labs(x = "Weeks", y = "Difference (Predicted - Observed counts)") +
  theme_minimal(base_size = 20) +
  labs(title = "Difference Obs vs Pred") +
  theme(legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank(), 
        legend.spacing.x = unit(0.5, "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=10)) +
  facet_wrap(~ Species, scales = "free_y", nrow = 3)

# Save the combined difference plot
ggsave("Poiss_difference_plot.jpg", plot = join_difference_plot, width = 12, height = 8, dpi = 300)





## Descriptives ------------------------------------------------------------


data_year<-data_19%>%group_by(Year,Species)%>%
  summarise(SD=sd(MeanCount),
            MeanCount=mean(MeanCount),
            LogCount=log(MeanCount),
            LogSD=log(SD),
            SampleSize=n())
data_year$SD[is.na(data_year$SD)]<-0


box1<-ggplot(data_year) +
  geom_boxplot(aes(x=factor(Year), y=LogCount))+
  theme_classic()+
  ggtitle("Yearly Elasmobranch Counts")+
  labs(y = "Counts")+
  labs(x = "Year")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~Species)

ggsave("BoxDescriptive.jpg", plot = box1, width = 12, height = 8, dpi = 300)

density1<-ggplot(data_year, aes(x=LogCount, fill = Species)) +
  geom_density(alpha = 0.5) +
  labs(x = "Log Count", y = "Density") +
  facet_wrap(~Species) +
  theme_minimal()+
  theme(legend.position = "n",
        # strip.background = element_blank(),
        strip.text = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("DensityDescriptive.jpg", plot = density1, width = 12, height = 8, dpi = 300)

point1<-ggplot(data_year, aes(x=factor(Year), y=LogCount, color= Species)) +
  geom_point(aes(y = LogCount), size = 3, alpha = 0.5) +
  labs(x = "Years", y = "Elasmobranch mean counts", color = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5))+
  theme(legend.position = "n") +
  facet_wrap(~ Species, scales = "free_y")

ggsave("PointDescriptive.jpg", plot = point1, width = 12, height = 8, dpi = 300)




## Heatmaps ----------------------------------------------------------------

#create trimed vector of covariates

covariates <- data_19 %>%
  pivot_wider(
    names_from = Species, # Species will become the column names
    values_from = MeanCount # MeanCount will populate the new columns
  )

covariates_trim<-covariates%>%dplyr::select(c("Blacktips" ,  "EagleRays", "Galapagos", "Hammerheads", "MantaRays", "MarbledRays", "MobulaRays", "Silky", "Silvertips", "TigerSharks", "Turtles", "WhaleSharks", "Whitetips"))


## Correlation Heatmap -------------------------------------------------------------------


#Reset data arrangements

data_19 <- read.csv("data.csv")[,-1]

data_19$refdate <- as.Date(data_19$refdate, format="%Y-%m-%d")

options(max.print=500)

data_19 <- data_19 %>%
  mutate(
    week = ifelse(week < 10, paste0("0", week), as.character(week)),
    Date = as.Date(refdate),
    YearWeek = paste(year(Date), week, sep = "-")
  ) %>%
  arrange(refdate)


data_19$week<-as.numeric(data_19$week)
#make week cyclic

data_19$week <- 2 * 3.1416 * data_19$week / 52  # Cyclic transformation for week

# Define combined sine and cosine transformation
data_19$week <- sin(data_19$week) + cos(data_19$week)

# Get the unique species
unique_species <- unique(data_19$Species)


# Compute the correlation matrix
cor_matrix <- cor(covariates_trim, use = "complete.obs")

# Convert the correlation matrix into a long format
cor_melt <- melt(cor_matrix)

# Define the mako color palette with the desired colors
col_palette <- c("white", "aliceblue", "azure",mako(5, direction = -1))  # Adjusting colors

# Add row and column indices for filtering upper and lower triangles
cor_melt <- cor_melt %>%
  mutate(
    Var1 = factor(Var1, levels = unique(Var1)),  # Ensure correct factor levels
    Var2 = factor(Var2, levels = unique(Var2)),
    row = as.numeric(Var1),
    col = as.numeric(Var2)
  )


####Colors top Text Bottom----
heatmap_plot <- ggplot() +
  # Lower triangle: Color + text
  geom_tile(data = cor_melt %>% filter(row <= col), aes(Var1, Var2, fill = value), color = "white") +
  geom_text(data = cor_melt %>% filter(row >= col), aes(Var1, Var2, label = round(value, 2)), size = 3, color = "black") +
  # Upper triangle: Blank tiles only
  geom_tile(data = cor_melt %>% filter(row < col), aes(Var1, Var2), fill = NA, color = "white") +
  # Color scale with the mako palette
  scale_fill_gradientn(
    colors = col_palette,         # Use the custom color palette
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),  # Rescale correlation values
    name = "Correlation",         # Legend title
    limits = c(-1, 1),            # Set limits for the color scale
    guide = "colourbar"           # Ensure it shows as a continuous color bar
  ) +
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Keep x-axis labels vertical
    axis.title = element_blank(),  # Remove axis titles
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.border = element_blank(),  # Remove panel borders
    axis.ticks = element_blank()  # Remove axis ticks
  ) +
  coord_fixed()  # Keep the aspect ratio square

# Display the heatmap
print(heatmap_plot)

ggsave("correlation_heatmap.jpg", plot = heatmap_plot, width = 8, height = 6)

## ModelHeatmap ------------------------------------------------------------

#Read the CSV file containing matrix-arranged mean beta values from model.

model_betas <- read.csv("model_betas.csv", header = TRUE, row.names = 1)

#Read CSV file containing numerical values for signifficant betas only
signif_betas <- read.csv("Signifficant.csv", header = TRUE, row.names = 1)

# Convert both to matrices
model_betas <- as.matrix(model_betas)
signif_betas <- as.matrix(signif_betas)

# Melt both matrices
cor_melt <- melt(model_betas)
signif_melt <- melt(signif_betas)

# Ensure correct factor levels for rows and columns
cor_melt <- cor_melt %>%
  mutate(
    Var1 = factor(Var1, levels = unique(Var1)),
    Var2 = factor(Var2, levels = unique(Var2))
  )

# Add significant values as a new column
cor_melt <- cor_melt %>%
  left_join(signif_melt, by = c("Var1", "Var2"), suffix = c("", "_signif"))

# Define the mako color palette with desired colors
col_palette <- c("white", "aliceblue", "azure", mako(5, direction = -1))  # Adjusting colors

# Create the heatmap
heatmap_plot <- ggplot(data = cor_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +  # Add white grid lines
  geom_text(
    data = filter(cor_melt, !is.na(value_signif)),  # Only show text for significant values
    aes(label = round(value, 2)),
    size = 3, 
    color = "black"
  ) +
  scale_fill_gradientn(
    colors = col_palette,  # Use the custom color palette
    values = scales::rescale(c(-1, -0.5, 0, 0.5, 1)),  # Rescale correlation values
    na.value = "black",
    name = "Effect",  # Legend title
    limits = c(-1, 1),  # Set limits for the color scale
    guide = "colourbar"  # Ensure it shows as a continuous color bar
  ) +
  labs(title = "",
       x = "Response variable",
       y = "Predictor variable") +
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Keep x-axis labels vertical
    axis.title.x = element_text(angle = 0, vjust = 0.5, hjust = 0,size=20),  # Remove axis titles
    axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0,size=20), 
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.border = element_blank(),  # Remove panel borders
    axis.ticks = element_blank()  # Remove axis ticks
  ) +
  coord_fixed()+  # Keep the aspect ratio square
  coord_flip()

# Print the heatmap
print(heatmap_plot)
ggsave("betas_heatmap_Signif.jpg", plot = heatmap_plot, width = 8, height = 6)

