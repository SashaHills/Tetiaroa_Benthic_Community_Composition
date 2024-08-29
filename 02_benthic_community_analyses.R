#### DESCRIPTION ####
# This script generates a Bray-Curtis uses marine 
# benthic community data collected in Tetiaroa Atoll, French Polynesia.
# It performs Principal Components Analysis (PCA) to explore patterns in 
# community composition and uses distance-based Redundancy Analysis (dbRDA)
# to assess relationships between benthic communities and biogeophysical 
# predictors. The data used in these analyses were pre-processed from 
# CoralNet annotations and spatially co-occurring biogeophysical data (see
# 00_covariate_preparation.py and 01_data_preprocessing.R). Data were collected
# by Rosalie Wright (University of Oxford) and Casey Benkwitt (Lancaster 
# University) from November 6–14, 2021.
# This script also contains the code for plotting the data and the purpose of each plot is explained throughout.

#### REQUIRED PACKAGES ####
# install packages (first run only)
install.packages(c("tidyverse", "dplyr", "readxl", "ggplot2",
                   "here", "usdm", "easypackages", "PNWColors", 
                   "corrplot", "Cairo", "vegan", "stats", "cowplot", "pheatmap", "ggrepel")) 

library(easypackages)
libraries("tidyverse", "dplyr", "readxl", "ggplot2",
          "here", "usdm", "easypackages", "PNWColors", 
          "corrplot", "Cairo", "vegan", "stats", "cowplot", "pheatmap", "ggrepel")
# multiple packages have the functions select() and filter(), tell R to 
# prioritize these functions from the dplyr package
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
set.seed(999) # for reproducibility 

#### BENTHIC DATA PREP ####
# check working directory
getwd()
# set working directory and point the here package to where we are
setwd("C:/tetiaroa_benthic_analysis/tetiaroa_2021/")
here::i_am(".here")
here::here() # double check that here is looking in the right place.
# here will allow us to use relative file paths moving forward.

#### LOAD BENTHIC COMMUNITY DATA ####
community = read.csv(here("Benthic_Analysis", "Data", 
                          "Benthic_Community_Composition_RQ1.csv"))

# convert depth from negative values to positive values for ease of 
# interpretation later.
community = community %>%
  mutate(Depth = (-1 * Depth))

#### CLEVELAND DOT PLOTS  ####
# checking for outliers in the explanatory variables
# depth
a = ggplot(community,
           aes(x = Depth, y = Motu)) +
  geom_point() +
  labs(x = "Depth (m)", y = "Motu") +
  theme_bw() +
  theme(panel.grid = element_blank()) 
# slope
b = ggplot(community,
           aes(x = Slope, y = Motu)) +
  geom_point() +
  labs(x = "Slope (degrees)", y = "Motu") +
  theme_bw() +
  theme(panel.grid = element_blank()) 
# offshore distance 
c = ggplot(community,
           aes(x = Offshore_distance, y = Motu)) +
  geom_point() +
  labs(x = "Distance offshore (m)", y = "Motu") +
  theme_bw() +
  theme(panel.grid = element_blank()) 
# N15
d = ggplot(community,
           aes(x = N15, y = Motu)) +
  geom_point() +
  labs(x = "Mean macroalgal δ15N within 200 m (‰)", y = "Motu") +
  theme_bw() +
  theme(panel.grid = element_blank()) 
# intercardinal bearing
e = ggplot(community,
           aes(x = Intercardinal_bearing, y = Motu)) +
  geom_point() +
  labs(x = "Intercardinal bearing", y = "Motu") +
  theme_bw() +
  theme(panel.grid = element_blank())

# arrange the plots into a multi-figure panel with two rows
cleveland_dotplots = plot_grid(
  plot_grid(a, b, c, ncol = 3), # first row with three plots
  plot_grid(d, e, NULL, ncol = 3), # second row with two plots and one empty space
  nrow = 2)
plot(cleveland_dotplots)

# save the multipanel figure to a file
ggsave(here("Benthic_Analysis", "Figures", "Cleveland_Dotplots.png"),
       cleveland_dotplots, width = 10, height = 8, units = "in", dpi = 450,
       bg = "white")

# remove the intermediate plots
rm(list = c("a", "b", "c", "d", "e"))

# exploratory histograms of some benthic classes 
par(mar = c(2,2,2,1))
hist(community$Acropora_arborescent)
hist(community$Favites_encrusting)
hist(community$Montipora_encrusting)
hist(community$Porites_rus)

#### DATA TRANSFORMATION ####
# apply a log (x + 1) transformation to the benthic cover columns to stabilize 
# variance and make the data more suitable for later multivariate analyses. for 
# ordination of our benthic communities, all species are measured in the same 
# units (percent cover), so the data do not need to be standardized.
community = community %>%
  mutate(across(Acropora_arborescent:Rubble, ~ log(. + 1)))

# update site names for consistency and clarity
community$Site = ifelse(community$Site == "north", "North", community$Site)
community$Site = ifelse(community$Site == "south", "South", community$Site)
community$Site = ifelse(community$Site == "site1", "Southwest", community$Site)
community$Site = ifelse(community$Site == "site2", "Southeast", community$Site)

# add a sample column by combining the motu + site information (where present)
community = community %>%
  mutate(Sample = ifelse(is.na(Site), Motu, paste(Motu, Site, sep = "_"))) %>%
  relocate(Sample, .after = Time)

# save names so we can use them as labels later
cover_labels = colnames(community[,10:26]) #species/category names
species_labels = colnames(community[,10:17]) # specific species names
site_labels = community[, 4] # site names
bearing_labels = community[, 31] # bearing names
slope_labels = community[, 28] # slope 
covariate_labels = colnames(community[,27:31]) # covariates (env variables)
print(covariate_labels)

# are there any NAs in the community data? there shouldn't be, but let's 
# double check to prevent errors later
colSums(is.na(community[1:12, 10:26])) # no problems

# isolate only the percent cover data
cover = community[,10:26]
rownames(cover) = site_labels # assign row labels
colnames(cover) = cover_labels # assign column labels

#### Principal Components Analysis (PCA) ####
# use prcomp() to conduct a PCA on the COVER data, all 
# species use the same scale and units so set scale = FALSE
cover_pca = prcomp(cover, scale = FALSE) 

# check the model summary results
summary(cover_pca)
# PC1  explains ~ 38% of the total variance, PC2 explains an additional
# ~ 21% of the total variance, such that the cumulative proportion explained
# by PC1 and PC2 is ~ 60%

# loadings indicate the strength and direction of the relationship between
# each original variable (species) and the principal components
cover_pca$rotation  # how each species contributes to each PC

# scores of each observation (sample) on each principal component provide a
# transformed representation of the data in terms of the principal components
cover_pca$x # how each sample (or observation) is represented by each PC


# make quick plots using base R, we can improve on these later
par(mar = c(4,4,4,4))

# arrows are plotted to show the directionality and angle of the descriptors 
# in the ordination, where:
# descriptors at 180 degrees of each other are negatively correlated;
# descriptors at 90 degrees of each other have zero correlation;
# descriptors at 0 degrees of each other are positively correlated.

# type 2 scaling (default): distances among objects are not approximations of Euclidean 
# distances; angles between descriptor (species) vectors reflect their
# correlations.
biplot(cover_pca, scaling = 2) 

# type 1 scaling: attempts to preserve the Euclidean distance (in multidimensional
# space) among objects (sites): the angles among descriptor (species) vector are
# not meaningful.
biplot(cover_pca, scaling = 1) 


# the broken-stick model retains components that explain more variance than 
# would be expected by randomly dividing the variance into p parts.
head(bstick(cover_pca))
screeplot(cover_pca, bstick = TRUE, type = "lines", legend = NULL) # continue despite legend error
legend("topright", legend = c("Ordination", "Broken Stick"), 
       col = c("black", "red"), lty = 1, cex = 0.5, horiz = FALSE)

##### plotting PCA in ggplot #####

install.packages("ggplot2")
install.packages("ggfortify")
install.packages("ggrepel")
library(ggplot2)
library(ggfortify)
library(ggrepel)
autoplot(cover_pca)
PC1 <- cover_pca$x[,1] # these are the two most important Principal Components taken from the PCA
PC2 <- cover_pca$x[,2] # using the $x column plots how each sample (or observation) is represented by each PC, this is effectively what is visualised in a PCA

# prepare PCA scores for sites
pca_data <- as.data.frame(cover_pca$x) # make PCA a dataframe
pca_data$Site <- site_labels # see code lines 127+, add sites to dataframe
pca_data$Slope <- slope_labels # add slope to dataframe
pca_data$Intercardinal_bearing <- bearing_labels # add intercardinal bearing to dataframe
pca_data$Intercardinal_bearing <- as.factor(pca_data$Intercardinal_bearing)
pca_data$Intercardinal_bearing <- factor(pca_data$Intercardinal_bearing, levels = c("ENE", "SSE", "NNE", "WSW", "WNW"))


# prepare PCA loadings for species
pca_loadings <- as.data.frame(cover_pca$rotation)
pca_loadings$Variable <- rownames(pca_loadings)

# scale loadings for better visibility 
pca_loadings$PC1 <- pca_loadings$PC1 * 2
pca_loadings$PC2 <- pca_loadings$PC2 * 2

print(pca_loadings)

# Create PCA biplot of just species loadings (contributions to PC1 and PC2 i.e. differences between sites)
pca_plot_species_loadings <- ggplot() +
  geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "red") +
  geom_text(data = pca_loadings, aes(x = PC1, y = PC2, label = Variable), color = "red", vjust = -0.5, hjust = 1.0, size = 2.5) +
  labs(title = "PCA Biplot: Species Contributions to Benthic Community Differences") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1)
print(pca_plot_species_loadings)

# PCA plot trial: sites coloured by slope, shapes based on intercardinal bearing

pca_plot1 <- ggplot(pca_data, aes(x = PC1, y = PC2, shape = Intercardinal_bearing, colour = Slope, label = Site)) +
  geom_point(size = 4) + 
  scale_shape_manual(values = c(19, 20, 15, 18, 17)) +
  scale_colour_gradientn(colours = c("darkblue","lightblue","green","yellow","red","darkred")) +
  xlab("PC1") + 
  ylab("PC2") +
  guides(colour=guide_colourbar(title = "Slope", barwidth = 0.9, barheight = 10)) +
  geom_text_repel(size = 2.0, color = "black", fontface = "bold",
                  nudge_y = 0.015,
                  nudge_x = 0.005) +
  #stat_ellipse(data = pca_data_2, aes(x = PC1, y = PC2), level = 0.95, color = "black", linetype = "dashed") +
  theme(panel.background = element_rect(fill = "white", color = NA),  # Set the background to white
        plot.background = element_rect(fill = "white", color = NA), # Set the plot area background to white
        panel.grid.major = element_line(color = "lightgrey"),  # Add major grid lines
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # Add a black border around the panel
        axis.line = element_line(color = "black"))  # Add axis lines
print(pca_plot1)

# PCA plot trial: sites coloured by intercardinal bearing

pca_plot2 <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = Intercardinal_bearing, label = Site)) +
  geom_point(size = 3) + 
  labs(x = "PC1 (38.3%)", 
       y = "PC2 (21.1%)",
       title = "PCA Plot",
       colour = "Intercardinal Bearing") +
  scale_colour_manual(values = c("lightblue", "blue", "darkgreen", "#99CC00", "#00FF00")) +
  geom_text_repel(size = 2.0, color = "black", fontface = "bold",
                  nudge_y = 0.015,
                  nudge_x = 0.005) +
  # stat_ellipse(data = pca_data_2, aes(x = PC1, y = PC2), level = 0.5, color = "black", linetype = "dashed") +
  theme(panel.background = element_rect(fill = "white", color = NA),  # Set the background to white
        plot.background = element_rect(fill = "white", color = NA), # Set the plot area background to white
        panel.grid.major = element_line(color = "lightgrey"),  # Add major grid lines
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # Add a black border around the panel
        axis.line = element_line(color = "black"))  # Add axis lines
print(pca_plot2)

##### adding environmental variables to the PCA #####
# using envfit in vegan package

# isolate the environmental columns from Depth to Intercardinal Bearing to create covariate dataframe named environment
environment = community[,27:31] 
environment$Intercardinal_bearing = as.factor(environment$Intercardinal_bearing) # convert intercardinal bearing to factor values

envfit_PCA <- envfit(cover_pca, environment, permutations = 999)
envfit_PCA_vectors<- as.data.frame(envfit_PCA$vectors$arrows)
envfit_PCA_vectors$Variable <- rownames(envfit_PCA_vectors)
envfit_PCA_vectors$Variable[envfit_PCA_vectors$Variable == "Offshore_distance"] <- "Offshore distance"  # Rename environmental variables for aesthestic

# now combine with previous PCA plots to customise overlay of environmental variables 

# basic plot of sites and covariates
env_pca_plot1 <- ggplot(cover_pca$x, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = environment$Site)) +  # Change to your site grouping variable
  geom_segment(data = as.data.frame(envfit_PCA$vectors$arrows),
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")),
               color = "red", linewidth = 1.2) +
  geom_text(data = as.data.frame(envfit_PCA$vectors$arrows),
            aes(x = PC1, y = PC2, label = rownames(envfit_PCA$vectors$arrows)), 
            hjust = 0.5, vjust = 0.5, color = "red") +
  geom_segment(data = as.data.frame(envfit_PCA$factors$centroids),
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")),
               color = "red", linewidth = 1.2) +
  geom_text(data = as.data.frame(envfit_PCA$factors$centroids),
            aes(x = PC1, y = PC2, label = rownames(envfit_PCA$factors$centroids)), 
            hjust = 0.5, vjust = 0.5, color = "red") +
  labs(title = "PCA Biplot with Environmental Variables",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

print(env_pca_plot1)

# now we want to change value names so we get neater plots
# rename rownames for bearings (factors$centroids)
rownames(envfit_PCA$factors$centroids)[rownames(envfit_PCA$factors$centroids) == "Intercardinal_bearingENE"] <- "Bearing ENE"
rownames(envfit_PCA$factors$centroids)[rownames(envfit_PCA$factors$centroids) == "Intercardinal_bearingSSE"] <- "Bearing SSE"
rownames(envfit_PCA$factors$centroids)[rownames(envfit_PCA$factors$centroids) == "Intercardinal_bearingNNE"] <- "Bearing NNE"
rownames(envfit_PCA$factors$centroids)[rownames(envfit_PCA$factors$centroids) == "Intercardinal_bearingWSW"] <- "Bearing WSW"
rownames(envfit_PCA$factors$centroids)[rownames(envfit_PCA$factors$centroids) == "Intercardinal_bearingWNW"] <- "Bearing WNW"

# rename rownames for covariates (vectors$arrows)
rownames(envfit_PCA$vectors$arrows)[rownames(envfit_PCA$vectors$arrows) == "Offshore_distance"] <- "Offshore distance"

# Rename Site data so it presents better

pca_data <- pca_data %>% 
  mutate(Site = recode(Site, 
                       "Aie_North" = "'Ā'ie (N)",
                       "Aie_South" = "'Ā'ie (S)",
                       "Rei_Southwest" = "Reiono (SW)",
                       "Rei_Southeast" = "Reiono (SE)",
                       "Rim_Southwest" = "Rimatu'u (SW)",
                       "Rim_Southeast" = "Rimatu'u (SE)",
                       "Aur" = "Ahuroa",
                       "Hir" = "Hīra'a'ānae",
                       "Hon" = "Honuea",
                       "One" = "Onetahi",
                       "Oro" = "Horoāterā",
                       "Tau" = "Tauini"))

# plot with sites labels, sites coloured by Intercardinal bearing and covariate vectors (arrows) displayed, all adjusted for aesthestic
env_pca_plot2 <- ggplot() +
  geom_segment(data = as.data.frame(envfit_PCA$vectors$arrows), # all covariates are vectors except intercardinal bearing so we do vectors and factors separately
               aes(x = 0, y = 0, xend = PC1, yend = PC2), # geom_segment plots arrows
               arrow = arrow(type = "open", length = unit(0.2, "cm")),
               color = "#B39DDB", linewidth = 0.8, size = 1) +
  geom_text_repel(data = as.data.frame(envfit_PCA$vectors$arrows), # use text repel to avoid overlapping vector labels
                  aes(x = PC1, y = PC2, label = rownames(envfit_PCA$vectors$arrows)), 
                  size = 5, color = "#B39DDB", nudge_x = 0.05, nudge_y = -0.02) +
  geom_segment(data = as.data.frame(envfit_PCA$factors$centroids),
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(type = "open", length = unit(0.2, "cm")),
               color = "#B39DDB", linewidth = 0.8, size = 1) +
  geom_text(data = as.data.frame(envfit_PCA$factors$centroids), aes(x = PC1, y = PC2, label = rownames(envfit_PCA$factors$centroids)), 
            size = 5, color = "#B39DDB",
            nudge_x = ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing WNW", 0.2, # here we use nudge if else to madjust position of specific labels 
                             ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing WSW", 0.23,
                                    ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing ENE", 0.01,
                                           ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing SSE", -0.12, 0)))),  # Example nudge for a specific label
            nudge_y = ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing WNW", -0.03, 
                             ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing WSW", 0.05,
                                    ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing NNE", 0.08,
                                           ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing ENE", -0.07,
                                                  ifelse(rownames(envfit_PCA$factors$centroids) == "Bearing SSE", 0.07, 0)))))) +
  geom_point(data = pca_data, aes(x = PC1, y = PC2, colour = Intercardinal_bearing), size = 2.5) + # geom_point plots points (i.e. sites)
  geom_text(data = pca_data, aes(x = PC1, y = PC2, label = Site), 
            size = 5, color = "black", fontface = "bold",
            nudge_x = ifelse(pca_data$Site == "Horoāterā", -0.07,
                             ifelse(pca_data$Site == "Rimatuꞌu (SW)", -0.04,
                                    ifelse(pca_data$Site == "'Ā'ie (N)", -0.02,
                                           ifelse(pca_data$Site == "Onetahi", -0.07,
                                                  ifelse(pca_data$Site == "Reiono (SW)", -0.05, 0))))),
            nudge_y = ifelse(pca_data$Site == "Reiono (SW)", 0.1,
                             ifelse(pca_data$Site == "Reiono (SE)", 0.1,
                                    ifelse(pca_data$Site == "Rimatuꞌu (SW)", -0.08,
                                           ifelse(pca_data$Site == "Rimatuꞌu (SE)", 0.1, 
                                                  ifelse(pca_data$Site == "'Ā'ie (N)", -0.1, 
                                                         ifelse(pca_data$Site == "'Ā'ie (S)", 0.1,
                                                                ifelse(pca_data$Site == "Ahuroa", 0.1,
                                                                       ifelse(pca_data$Site == "Hīra'a'ānae", 0.1,
                                                                              ifelse(pca_data$Site == "Honuea", 0.1,
                                                                                     ifelse(pca_data$Site == "Tauini", 0.1,
                                                                                            ifelse(pca_data$Site == "Horoāterā", -0.06, 
                                                                                                   ifelse(pca_data$Site == "Onetahi", -0.07, 0))))))))))))) +
  labs(x = "PC1 (38.86%)", 
       y = "PC2 (21.42%)",
       title = "",
       colour = "Intercardinal Bearing",
       size = 16) +
  scale_colour_manual(values = c("#99CC00", "#43A047", "#64B5F6", "#3366FF", "navy")) +
  
  # stat_ellipse(data = pca_data_2, aes(x = PC1, y = PC2), level = 0.5, color = "black", linetype = "dashed") +
  theme(panel.background = element_rect(fill = "white", color = NA),  # Set the background to white
        plot.background = element_rect(fill = "white", color = NA), # Set the plot area background to white
        panel.grid.major = element_blank(),  # Add major grid lines
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # Add a black border around the panel
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = "top",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 17))  
print(env_pca_plot2)
ggsave(here("Benthic_Analysis", "Figures", "env_pca_plot2.png"),
       env_pca_plot2, width = 10, height = 10, units = "in", dpi = 450,
       bg = "white")

#### Analysis of Species Contributions to PCA ####

# The following analyses and plots visualise the differences in site community composition based on cover species/categories

# Identify species with highest loadings for PC1 and PC2
top_species_pc1 <- pca_loadings[order(abs(pca_loadings$PC1), decreasing = TRUE), ]
top_species_pc2 <- pca_loadings[order(abs(pca_loadings$PC2), decreasing = TRUE), ]
# Print top species
print("Top species for PC1:")
print(top_species_pc1)
# output:PC1          
# Rubble                      1.20212764
# Turf                        0.90833016
# Pavement                    0.78976688 
# Sand                       -0.75985195 
# Sand_rubble                 0.37233785 
# Montipora_encrusting       -0.36123374 (only included top five in manuscript)

print("Top species for PC2:")
print(top_species_pc2)
# output:           PC2  
# Macroalgae                   1.505193325  
# Porites_rus                 -0.657822097  
# Rubble                       0.625328203 
# Sand_rubble                 -0.542631917 
# Turf                        -0.395431849


# Perform ANOVA for top species identified in PC1 and PC2 loadings - CAN DO THIS IF NEEDED FOR MANUSCRIPT (IF SIGNIFICANCE OF EACH INDIVIDUAL SPECIES CONTRIBUTION IS REQUIRED)
# Example for one species: 'Acropora_arborescent'
# anova_pc1 <- aov(PC1 ~ Acropora_arborescent, data = pca_loadings)
# summary(anova_pc1)

# anova_pc2 <- aov(PC2 ~ Acropora_arborescent, data = pca_loadings)
# summary(anova_pc2)

# OR RUN THIS CODE IF NEEDED

# Function to run ANOVA for each species
# run_anova <- function(species, pc, pca_loadings) {
#  formula <- as.formula(paste(pc, "~", species))
#  aov_result <- aov(formula, data = pca_loadings)
#  return(summary(aov_result))
# }

# Run ANOVA for top species for PC1
# for(species in top_species_pc1$Species) {
#  cat("ANOVA for", species, "on PC1:\n")
#  print(run_anova(species, "PC1", data))
}

# Run ANOVA for top species for PC2
#for(species in top_species_pc2$Species) {
#  cat("ANOVA for", species, "on PC2:\n")
#  print(run_anova(species, "PC2", data))
#}

##### Species-Site Plots #####
# Here we plot and visualise the species/category abundance at each site using COVER dataframe

#recode species names so they present better
pca_loadings <- pca_loadings %>%
  mutate(Variable = recode(Variable, 
                           "Acropora_arborescent" = "Acropora arborescent",
                           "Acropora_corymbose" = "Acropora corymbose",
                           "Dead_coral" = "Dead Coral",
                           "Favites_encrusting" = "Favites encrusting",
                           "Live_benthos" = "Live benthos",
                           "Montipora_encrusting" = "Montipora encrusting",
                           "Pavona_cactus" = "Pavona cactus",
                           "Pavona_varians" = "Pavona varians",
                           "Porites_massive_submassive" = "Porites massive/submassive",
                           "Porites_rus" = "Porites rus",
                           "Sand_rubble" = "Sand/rubble"))

# Plot of species loadings in PCA space
pca_plot_species_loadings <- ggplot() +
  geom_segment(data = pca_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "navy") +
  labs(title = "PCA Biplot Displaying Species Contributions to Benthic Community Differences", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text_repel(data = pca_loadings, aes(x = PC1, y = PC2, label = Variable), size = 2.0, color = "black", fontface = "bold", nudge_x = -0.01, nudge_y = -0.02 ) +
  coord_fixed(ratio = 1)
print(pca_plot_species_loadings)

###### Heatmap of species abundance at each site ######

devtools::install_github("jokergoo/ComplexHeatmap")
install.packages("circlize")
install.packages("RColorBrewer")
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# original percent cover data (%)
original_data = read.csv(here("Benthic_Analysis", "Data", "Benthic_Community_Composition_RQ1.csv"))
original_data = original_data %>%
  select(Acropora_arborescent:Rubble)

cover_matrix <- as.matrix(original_data) # turn cover data into matrix

# Rename columns using colnames()
colnames(cover_matrix) <- c("Acropora arborescent", "Acropora corymbose", 
                            "Favites encrusting", "Montipora encrusting", 
                            "Pavona cactus", "Pavona varians", 
                            "Porites rus", "Porites massive/submassive",
                            "Sand", "Dead coral", "Pavement", "Sand rubble",
                            "CCA", "Live benthos", "Turf", "Macroalgae", "Rubble")

# Assign new row names using rownames()
rownames(cover_matrix) <- c("'Ă'ie (N)", "'Ă'ie (S)", "Ahuroa", "Hīra'a'ānae", "Honuea",
                            "Onetahi", "Horoāterā", "Reiono (SW)", "Reiono (SE)", 
                            "Rimatu'u (SW)", "Rimatu'u (SE)", "Tauini")

# Create the heatmap with ComplexHeatmap
heatmap = Heatmap(cover_matrix, 
                  col = rev(brewer.pal(n = 10, name = "Spectral")), 
                  cluster_rows = TRUE, 
                  cluster_columns = TRUE,  
                  row_names_side = "left",
                  row_dend_side = "right", 
                  column_names_side = "top", 
                  column_dend_side = "bottom",
                  rect_gp = gpar(col = "black", lwd = 0.75),
                  column_names_gp = gpar(fontsize = 8, fontface = "italic"),
                  row_names_gp = gpar(fontsize = 8),
                  heatmap_legend_param = list(title = "Percent cover (%)",
                                              labels_gp = gpar(fontsize = 6),   # Change legend text size
                                              title_gp = gpar(fontsize = 8,
                                                              fontface = "bold")))
plot(heatmap)

# Export the heatmap as a png
png(filename = "Percent_Cover_Heatmap.png", width = 1800, height = 1200, res = 300, bg = "white")
ComplexHeatmap::draw(heatmap)
dev.off()  # Close the device


###### Facet Grid Plot ######

#create site name data to add to end dataframe

cover_long <- reshape2::melt(cover) # melt the cover data to prepare for facetwrap
colnames(cover_long) <- c("Species", "Abundance") # only Species and Abundance were transferred, label these as the columns
cover_long$Site <- factor(rep(site_labels, times = 17)) # add a site column using site_labels vector created in line 130 , which repeats 18 times (corresponsing with the number of species/categories we have)
# this ensures that in the dataset we have the abundance for each species/category abundance at each site, in separate columns

# Create the facet grid plot, separating each site (showing species abundane broken down by site)
ggplot(cover_long, aes(x = Species, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Species Composition at Each Site") +
  ylab("Abundance") +
  xlab("Species")

# there are a lot of colours going on that are not necessary if corresponding with x axis labels (species)
# therefore we are going to colour them based on the categories they fall into in CoralNet...

# Create a vector that associates each species to its category
cover_categories <- c(
  "Acropora_arborescent" = "Hard Coral",
  "Acropora_corymbose" = "Hard Coral",
  "Favites_encrusting" = "Hard Coral",
  "Montipora_encrusting" = "Hard Coral",
  "Pavona_cactus" = "Hard Coral",
  "Pavona_varians" = "Hard Coral",
  "Porites_rus" = "Hard Coral",
  "Porites_massive_submassive" = "Hard Coral",
  "Sand" = "Soft Substrate",
  "Dead_coral" = "Hard Substrate",
  "Pavement" = "Hard Substrate", 
  "Sand_rubble" = "Other",
  "CCA" = "Algae",
  "Live_benthos" = "Other Invertebrates",
  "Turf" = "Algae",
  "Macroalgae" = "Algae",
  "Rubble" = "Other")

# Create a new data frame to map species names to categories
species_categories <- data.frame(
  Category = cover_categories,
  stringsAsFactors = FALSE
)
# make the rownames a column called Species
species_categories$Species <- rownames(species_categories) # make a columns for species (rather than species being rownames)

# create a vector that can be used to produce the order we want in the legend
# legend items must be ordered as factors in the dataframe i.e. must be done before plotting
desired_order <- c("Hard Coral", "Hard Substrate", "Soft Substrate", "Other Invertebrates", "Algae", "Other")

# Merge the category mapping with the cover data
cover_long <- merge(cover_long, species_categories, by = "Species")
cover_long$Category <- factor(cover_long$Category, levels = desired_order) # change the category order using vector above so the legend appears how we want

# create a colour palette for category
category_colours <- c("Hard Coral" = "#D1C4E9", "Algae" = "#C8E6C9", "Hard Substrate" = "#FFCDD2", "Soft Substrate" = "#FFE082", "Other Invertebrates" = "#FFF59D", "Other" = "#CCCCCC")

# rename site and species so they present better
# dplyr approach
cover_long <- cover_long %>%
  mutate(Species = recode(Species, 
                          "Acropora_arborescent" = "Acropora arborescent",
                          "Acropora_corymbose" = "Acropora corymbose",
                          "Dead_coral" = "Dead Coral",
                          "Favites_encrusting" = "Favites encrusting",
                          "Live_benthos" = "Live benthos",
                          "Montipora_encrusting" = "Montipora encrusting",
                          "Pavona_cactus" = "Pavona cactus",
                          "Pavona_varians" = "Pavona varians",
                          "Porites_massive_submassive" = "Porites massive/submassive",
                          "Porites_rus" = "Porites rus",
                          "Sand_rubble" = "Sand/rubble"))

cover_long <- cover_long %>%
  mutate(Site = recode(Site, 
                       "Aie_north" = "Aie (N)",
                       "Aie_South" = "Aie (S)",
                       "Rei_Southwest" = "Rei (SW)",
                       "Rei_Southeast" = "Rei (SE)",
                       "Rim_Southwest" = "Rim (SW)",
                       "Rim_Southeast" = "Rim (SE)"))

# Create the facet grid plot with bars coloured by Category (as used in CoralNet)
# facet wrap grid with desired legend order and colours 
barplot_species_sites <- ggplot(cover_long, aes(x = Species, y = Abundance, fill = Category)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.25, colour = "black")) +
  scale_fill_manual(values = c(category_colours)) +
  ggtitle("Species Composition at Each Site by Category") +
  ylab("Abundance") +
  xlab("Species") 
print(barplot_species_sites)

#### Bray-curtis Dissimilarity and Plots - NOT USED FOR ANALYSIS, JUST FOR INTEREST/COMPLETENESS ####
# Cover data is converted to bray-curtis automatically in dbRDA using the distance = bray order

# convert the benthic community data to a Bray-Curtis dissimilarity matrix
bray_curtis_dist = vegdist(cover, method = "bray")
print(bray_curtis_dist) # view the matrix

##### Dendogram visualisation of Bray Curtis #####

# Perform hierarchical clustering
hc <- hclust(bray_curtis_dist, method = "average")

# Plot dendrogram
bray_curtis_dendogram <- plot(hc, main = "Hierarchical Clustering Dendrogram (Bray-Curtis)", xlab = "", sub = "", cex = 0.9)
print(bray_curtis_dendogram)

#### Distance-based Redundancy Analysis (dbRDA) ####

# Perform dbRDA with Bray-Curtis distance
# using environment data to include the covariates
dbrda_model <- capscale(cover ~ Depth + Slope + Offshore_distance + N15 + Intercardinal_bearing, data = environment, distance = "bray", na.action = na.fail, 
                        sqrt.dist = TRUE)

# Extract species scores
species_scores <- scores(dbrda_model, display = "species")

# View species scores
print(species_scores)

# basic plot
plot(dbrda_model)

# Extract other scores
site_scores <- scores(dbrda_model, display = "sites") # extracts site scores
constraint_scores <- scores(dbrda_model, display = "bp") # extracts biplot scores for the constraints (predictor variables)

# Convert to data frames for ggplot
species_scores_df <- as.data.frame(species_scores)
species_scores_df$species <- rownames(species_scores_df)

site_scores_df <- as.data.frame(site_scores)
site_scores_df$site <- rownames(site_scores_df)

constraint_scores_df <- as.data.frame(constraint_scores)
constraint_scores_df$constraint <- rownames(constraint_scores_df)

# Plotting with ggplot2
# Initial plot of sites and species and covariate contributions to dbRDA1 and dbRDA2
dbrda_plot <- ggplot() +
  geom_point(data = site_scores_df, aes(x = CAP1, y = CAP2), color = "blue", size = 3) +
  geom_text(data = site_scores_df, aes(x = CAP1, y = CAP2, label = site), vjust = -1, hjust = 1, color = "blue") +
  geom_point(data = species_scores_df, aes(x = CAP1, y = CAP2), color = "red", size = 3) +
  geom_text(data = species_scores_df, aes(x = CAP1, y = CAP2, label = species), vjust = -1, hjust = 1, color = "red") +
  geom_segment(data = constraint_scores_df, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.2, "cm")), color = "darkgreen") +
  geom_text(data = constraint_scores_df, aes(x = CAP1, y = CAP2, label = constraint), vjust = -0.5, color = "darkgreen") +
  theme_minimal() +
  labs(x = "dbRDA1", y = "dbRDA2", title = "Test dbRDA Plot")
print(dbrda_plot)
# this plot is very busy but shows all information: arrows for constraints and points for species however could switch
# can also show just the influence of constraints filtered by backward selection if we want

##### dbRDA Significance Testing #####

anova(dbrda_model) # overall test of the significance of the analysis
anova(dbrda_model, by = "axis") # test axes for significance, which axis/vector is significant
anova(dbrda_model, by = "terms") # test for significant environmental variables, THIS IS OF MOST INTEREST TO US
# OUTPUT was depth (p < 0.01), slope (p < 0.01), offshore distance (p <0.001), intercardinal bearing (p < 0.1)

# this output provides the sum of squares and p value but not the r^2 so we want to pull this out 

##### Calculating R^2 for dbRDA model #####

#using sum of squares (SS) that was calculated above
total_SS <- sum(dbrda_model$tot.chi)
explained_SS <- sum(dbrda_model$CCA$tot.chi)
residual_SS <- sum(dbrda_model$CA$tot.chi)

# Calculate R^2
R2 <- explained_SS / total_SS

# Print results
cat("Total Sum of Squares (TSS):", total_SS, "\n")

cat("Explained Sum of Squares (ESS):", explained_SS, "\n")

cat("Residual Sum of Squares (RSS):", residual_SS, "\n")

cat("R^2:", R2, "\n")
# R^2: 0.8359056 

##### Calculating the R^2 for each env variable #####

# Perform an ANOVA-like permutation test for each variable 
anova_dbrda_margin <- anova(dbrda_model, by = "margin", permutations = 999)
print(anova_dbrda_margin)

# Extract total sum of squares
total_SS_margin <- sum(dbrda_model$tot.chi)

# Extract sum of squares for each environmental variable from the anova_margin results
ss_per_variable <- anova_dbrda_margin[, "SumOfSqs"]

# Calculate R^2 for each variable
R2_per_variable <- ss_per_variable / total_SS_margin

# Print R^2 for each variable
print(R2_per_variable)

# Output:
# Depth: 0.07526694, Slope: 0.08944483, Offshore_distance: 0.07866974, N15: 0.05685790, Intercardinal_bearing: 0.31562040, Residual: 0.16409438

#### Backward selection ####

# perform backward selection using ordistep() to identify which biophysical
# variables are most important in explaining the variation in the distance matrix.
selected_model = ordistep(dbrda_model, scope = formula(dbrda_model), 
                          direction = "backward")
summary(selected_model)
# slope and intercardinal bearing were pulled out as most important (both p <0.01)

# extract the formula of the selected model
selected_formula = formula(selected_model)
print(selected_formula)

# extract the final retained variables
final_vars = attr(terms(selected_model), "term.labels")
print(final_vars)

# plot the final model for visualization
plot(selected_model)

##### Backward Selection Significance Testing #####

# get the ANOVA table for the final selected model
anova_table = anova(selected_model)
print(anova_table)
# model is signficiant (p < 0.001)

# this output provides the sum of squares and p value but not the R^2 so we want to pull this out 

##### Calculating R^2 for Backward Selection Model #####

total_SS_backward_selection_model <- sum(selected_model$tot.chi)
explained_SS_backward_selection_model <- sum(selected_model$CCA$tot.chi)
residual_SS_backward_selection_model <- sum(selected_model$CA$tot.chi)

# Calculate R^2
R2_backward_selection_model <- explained_SS_backward_selection_model / total_SS_backward_selection_model

# Print results
cat("Total Sum of Squares (TSS):", total_SS_backward_selection_model, "\n")

cat("Explained Sum of Squares (ESS):", explained_SS_backward_selection_model, "\n")

cat("Residual Sum of Squares (RSS):", residual_SS_backward_selection_model, "\n")

cat("R^2:", R2_backward_selection_model, "\n")
# R^2 for backward selection model = 0.654344

##### Calculating r^2 for each env variable #####

# Perform an ANOVA-like permutation test for each variable (similar to above)
anova_margin_backward_selection <- anova(selected_model, by = "margin", permutations = 999)
print(anova_margin_backward_selection)
# Significance: Slope (p< 0.05), intercardinal bearing (p < 0.001)

# Extract total sum of squares
selection_total_SS <- sum(selected_model$tot.chi)

# Extract sum of squares for each environmental variable from the anova_margin results
selection_ss_per_variable <- anova_margin_backward_selection[, "SumOfSqs"]

# Calculate R^2 for each variable
selection_R2_per_variable <- selection_ss_per_variable / selection_total_SS

# Print R^2 for each variable
print(selection_R2_per_variable)

# Output:
# Slope: 0.1231472, Intercardinal bearing: 0.5244532, Residual: 0.3456560
