##############################
## Skmer: Distance analysis ##
##############################

## See: http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2
## See: https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html


## Set working directory
#setwd("Path_to_the_working_directory")


## Load the libraries
library(ggplot2)
library(ggcorrplot)
library(spaa)
library(ape)
library(extrafont)


## Load the data
data <- read.table("skmer_distances.txt", header = TRUE, sep = "\t", 
                   dec = ".", row.names = 1)
corr_matrix <- as.matrix(round(data,3))

# Make sure it is in the right format
ggcorrplot(corr_matrix)


## Plot the matrix
# Plot a simple matrix
corr_plot <- ggcorrplot(corr_matrix)
corr_plot

# Add the correlation coefficient
corr_plot <- ggcorrplot(corr_matrix, lab = TRUE)
corr_plot

# Plot only the upper triangle
corr_plot <- ggcorrplot(corr_matrix, lab = TRUE, type = "upper")
corr_plot

# Change the number of digits in the label
corr_plot <- ggcorrplot(corr_matrix, lab = TRUE, type = "upper", digits = 3, 
                        lab_size = 4)
corr_plot

# Change the scale
corr_plot <- corr_plot + scale_fill_gradient2(limit = c(0,0.04), 
                                 low = "white", high =  "red", mid = "orange", 
                                 midpoint = 0.02, name = "Distance")

corr_plot

# Change the font of the legend
corr_plot + theme(text=element_text(size=12,  family = "Arial"))

