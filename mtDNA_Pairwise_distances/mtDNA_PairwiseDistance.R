###############################
## mtDNA: Pairwise Distances ##
###############################

## Set working directory
#setwd("Path_to_the_working_directory")


## Load libraries
library(ggplot2)


## Load data
distances <- read.table("mtDNA_PairwiseDistance.txt", header = TRUE, sep = "\t", dec = ",")


## Plot data
# Basic plot to see how it looks. First sort the categories in the x axis (default = alphabetical)
distances$Clade <- factor(distances$Clade,levels = c("Blue", "Brown", "Green", "Orange", "Red", "Violet", "Outgroups"))

d_plot <- ggplot(data = distances, aes(x = Clade, y = Distance, fill = Parameter)) + 
                 geom_bar(stat="identity", color="black", position=position_dodge())
d_plot

# Change the color palette
d_plot <- d_plot + scale_fill_grey()
d_plot

# Change the colors of x labels to match the mitochondrial clades
d_plot <- d_plot + theme(axis.text.x = element_text(size = 15, face = "bold", 
                                       colour= c("#2000de", "#a46301", "#09af02", "#fba400", "#ff0600", "#b902cc",  "black")))
d_plot

# Remove title of x axis and legend
d_plot <- d_plot + theme(axis.title.x = element_blank(), legend.title = element_blank())
d_plot

# Change the font of y labels and title
d_plot <- d_plot + theme(axis.text.y = element_text(size = 12, face = "bold"), 
                         axis.title.y = element_text(size = 15, face = "bold"))
d_plot

# Change legend position and size
d_plot <- d_plot + theme(legend.position="top", legend.text = element_text(size=15))
d_plot <- d_plot + theme(legend.margin = margin(-10, -10, -10, -10), legend.box.margin=margin(5, 5, 5, 5))
d_plot

# Add axis lines. Note: the x axis is not added as an axis but as a horizontal line
# crossing y at 0
d_plot <- d_plot + theme(axis.line.y = element_line(size = 0.5, colour = "black", linetype=1)) + 
          geom_hline(yintercept = 0, size = 0.5)
d_plot

# change scale in y axis
d_plot <- d_plot + scale_y_continuous(breaks = pretty(distances$Distance, n = 20), 
                                      limits = c(0, 0.095))
d_plot

