##############################
## Create distribution maps ##
##############################

## Set working directory
#setwd("Path_to_the_working_directory")


## Load libraries
library(ggplot2)
theme_set(theme_bw())
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(tmap)
library(spData)


## Load data
world <- ne_countries(scale = "medium", returnclass = "sf")

sampling_points <- st_read("Lautoconus ventricosus.Map.kml")

# Define the colors to be used
clade <- vector(mode="character", length=109)
clade[grep("Blue", sampling_points$description)] <- "#2000de"
clade[grep("Brown", sampling_points$description)] <- "#a46301"
clade[grep("Green", sampling_points$description)] <- "#09af02"
clade[grep("Orange", sampling_points$description)] <- "#fba400"
clade[grep("Red", sampling_points$description)] <- "#ff0600"
clade[grep("Violet", sampling_points$description)] <- "#b902cc"
sampling_points$Clade <- clade
  
sites <- st_as_sf(sampling_points, coords = c("X", "Y"), crs = 4326, agr = "constant")

## Make subsets of the data
sites_blue <- subset(sites, sites$Clade == "#2000de")
sites_brown <- subset(sites, sites$Clade == "#a46301")
sites_green <- subset(sites, sites$Clade == "#09af02")
sites_orange <- subset(sites, sites$Clade == "#fba400")
sites_red <- subset(sites, sites$Clade == "#ff0600")
sites_violet <- subset(sites, sites$Clade == "#b902cc")

# Note: these names were defined in an earlier version of the manuscript where
# the three species were spatially isolated. That is not true anymore, but I kept
# the names
Western_species <- subset(sites, sites$Clade == "#b902cc")
Central_species <- subset(sites, sites$Clade == "#09af02")
Eastern_species <- rbind(sites_blue, sites_orange, sites_red)
Eastern_species$Clade <- "#c9c913"

All_species <- rbind(Western_species, Eastern_species, Central_species)


## Create map
# Create a basic map
ggplot(data = world) + geom_sf(color = "grey77 ", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue'))

# Focus on the Mediterranean Sea
ggplot(data = world) + geom_sf(color = "grey77 ", fill = "#faf2e3") + 
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE) + 
  theme(panel.background = element_rect(fill = 'aliceblue'))

# Add sampling points
ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = sites, size = 2.5, shape = 23, fill = "#657b83") +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

# Add each clade color sequentially, so you can choose the order of the colors
ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = sites_violet, size = 2.5, shape = 23, fill = sites_violet$Clade) +
  geom_sf(data = sites_blue, size = 2.5, shape = 23, fill = sites_blue$Clade) +
  geom_sf(data = sites_orange, size = 2.5, shape = 23, fill = sites_orange$Clade) +
  geom_sf(data = sites_green, size = 2.5, shape = 23, fill = sites_green$Clade) +
  geom_sf(data = sites_red, size = 2.5, shape = 23, fill = sites_red$Clade) +
  geom_sf(data = sites_brown, size = 2.5, shape = 23, fill = sites_brown$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

# Make one map per clade
ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = sites_red, size = 2.5, shape = 23, fill = sites_red$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = sites_green, size = 2.5, shape = 23, fill = sites_green$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = sites_orange, size = 2.5, shape = 23, fill = sites_orange$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = sites_blue, size = 2.5, shape = 23, fill = sites_blue$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = sites_violet, size = 2.5, shape = 23, fill = sites_violet$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = sites_brown, size = 2.5, shape = 23, fill = sites_brown$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

# Make one map per species
ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = Western_species, size = 2.5, shape = 23, fill = Western_species$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = Central_species, size = 2.5, shape = 23, fill = Central_species$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)

ggplot(data = world) + geom_sf(color = "grey77", fill = "#faf2e3") + 
  theme(panel.background = element_rect(fill = 'aliceblue')) + 
  geom_sf(data = Eastern_species, size = 2.5, shape = 23, fill = Eastern_species$Clade) +
  coord_sf(xlim = c(36.50, -9.90), ylim = c(30.20, 46.00), expand = FALSE)


# To save the last map
#ggsave("File_name.pdf")

