################################
## MCMCtreeR: Plot chronogram ##
################################

## Set working directory
#setwd("")


## Load libraries
library(MCMCtreeR, quietly = TRUE, warn.conflicts = FALSE)


################################################################################
################################################################################

## Mitochondrial chronogram

## Load data
mtDNA_tree <- readMCMCtree(inputPhy = "Chronogram_mtDNA.tre", from.file = TRUE)
mcmc_mtDNA <- read.table(file = "XX.mcmc", head=TRUE)


## Plot tree
# Write the node ages to a file
write.table(mtDNA_tree$nodeAges, file = "Node_ages_mtDNA.txt", sep = "\t")

# Just plot the tree
MCMC.tree.plot(phy = mtDNA_tree, cex.tips = 1, time.correction = 1, 
               scale.res = c("Epoch", "Age"), plot.type = "phylogram", 
               cex.age = 0.6, cex.labels = 0.6, relative.height = 0.08, 
               col.tree = "black", label.offset = 4, node.method = "none", 
               no.margin = TRUE, edge.width = 4)

# Add the uncertainty to each node
MCMC.tree.plot(phy = mtDNA_tree, cex.tips = 1, time.correction = 1, 
               scale.res = c("Epoch", "Age"), plot.type = "phylogram", 
               cex.age = 1, cex.labels = 1, relative.height = 0.08, 
               col.tree = "black", label.offset = 0.5, node.method = "bar", 
               lwd.bar = 5, col.age = "red", no.margin = TRUE, edge.width = 4)

# Add the full distribution of the posterior to each node
MCMC.tree.plot(phy = mtDNA_tree, cex.tips = 1, time.correction = 1, MCMC.chain = mcmc_mtDNA, 
               scale.res = c("Epoch", "Age"), plot.type = "distributions", 
               cex.age = 1, cex.labels = 1, relative.height = 0.08, 
               col.tree = "black", label.offset = 0.5, density.col = "#00000050", 
               density.border.col = "#00000080", no.margin = TRUE, edge.width = 4)

# Color the tips according to the species:
# blue:   #2000de
# brown:  #a46301
# green:  #09af02
# orange: #fba400
# red:    #ff0600
# violet: #b902cc

tip_colors <- c("#2000de", "#2000de", "#2000de", "#2000de", "#2000de", "#2000de", 
                "#b902cc", "#b902cc", "#b902cc", "#b902cc", "#b902cc", "#ff0600", 
                "#ff0600", "#ff0600", "#fba400", "#fba400", "#fba400", "#a46301", 
                "#a46301", "#09af02", "#09af02", "#09af02", "black", "black", 
                "black", "black", "black", "black", "black", "black", "black", "black")

MCMC.tree.plot(phy = mtDNA_tree, cex.tips = 1, time.correction = 1, MCMC.chain = mcmc_mtDNA, 
               scale.res = c("Epoch", "Age"), plot.type = "distributions", 
               cex.age = 1, cex.labels = 1, relative.height = 0.08, 
               label.offset = 0.5, density.col = "#00000050", 
               density.border.col = "#00000080", no.margin = TRUE, edge.width = 4,
               tip.color = tip_colors, font = 2)

# Color the branches as per the mitochondrial clades
# blue:   #2000de
# brown:  #a46301
# green:  #09af02
# orange: #fba400
# red:    #ff0600
# violet: #b902cc

branch_colors = c("black", "black", "black", "black", "black", "black", "black", 
                  "#2000de", "#2000de", "#2000de", "#2000de", "#2000de", "#2000de", 
                  "#2000de", "#2000de", "#2000de", "#2000de", "#2000de", "#b902cc", 
                  "#b902cc", "#b902cc", "#b902cc", "#b902cc", "#b902cc", "#b902cc", 
                  "#b902cc", "#b902cc", "#ff0600", "#ff0600", "#ff0600", "#ff0600", 
                  "#ff0600", "#fba400", "#fba400", "#fba400", "#fba400", "#fba400", 
                  "#a46301", "#a46301", "#a46301", "#09af02", "#09af02", "#09af02", 
                  "#09af02", "#09af02", "black", "black", "black", "black", "black",
                  "black", "black", "black", "black", "black", "black", "black", 
                  "black", "black", "black", "black", "black")

MCMC.tree.plot(phy = mtDNA_tree, cex.tips = 1, time.correction = 1, MCMC.chain = mcmc_mtDNA, 
               scale.res = c("Epoch", "Age"), plot.type = "distributions", 
               cex.age = 1, cex.labels = 1, relative.height = 0.08, 
               label.offset = 0.5, density.col = "#00000050", 
               density.border.col = "#00000080", no.margin = TRUE, edge.width = 4,
               font = 2, col.tree = branch_colors)

################################################################################
################################################################################

## Nuclear chronogram

## Load data
mtDNA_tree <- readMCMCtree(inputPhy = "Chronogram_nuclear.tre", from.file = TRUE)
mcmc_nuclear <- read.table(file = "XX.mcmc", head=TRUE)


## Plot tree
# Write the node ages to a file
write.table(my_tree$nodeAges, file = "Node_ages_nuclear.txt", sep = "\t")

# Just plot the tree
MCMC.tree.plot(phy = my_tree, cex.tips = 1, time.correction = 1, 
               scale.res = c("Epoch", "Age"), plot.type = "phylogram", 
               cex.age = 0.6, cex.labels = 0.6, relative.height = 0.08, 
               col.tree = "black", label.offset = 4, node.method = "none", 
               no.margin = TRUE, edge.width = 4)

# Add the uncertainty to each node
MCMC.tree.plot(phy = my_tree, cex.tips = 1, time.correction = 1, 
               scale.res = c("Epoch", "Age"), plot.type = "phylogram", 
               cex.age = 1, cex.labels = 1, relative.height = 0.08, 
               col.tree = "black", label.offset = 0.5, node.method = "bar", 
               lwd.bar = 5, col.age = "red", no.margin = TRUE, edge.width = 4)

# Add the full distribution of the posterior to each node
MCMC.tree.plot(phy = my_tree, cex.tips = 1, time.correction = 1, MCMC.chain = mcmc_nuclear, 
               scale.res = c("Epoch", "Age"), plot.type = "distributions", 
               cex.age = 1, cex.labels = 1, relative.height = 0.08, 
               col.tree = "black", label.offset = 0.5, density.col = "#00000050", 
               density.border.col = "#00000080", no.margin = TRUE, edge.width = 4)

# Color the tips according to the species:
# green:  #09af02
# violet: #b902cc
# yellow: #ecfe36

tip_colors <- c("#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", 
                "#ecfe36", "#ecfe36", "#ecfe36", "#b902cc", "#b902cc", "#b902cc", 
                "#b902cc", "#b902cc", "#09af02", "#09af02", "black", "black")

MCMC.tree.plot(phy = my_tree, cex.tips = 1, time.correction = 1, MCMC.chain = mcmc_nuclear, 
               scale.res = c("Epoch", "Age"), plot.type = "distributions", 
               cex.age = 1, cex.labels = 1, relative.height = 0.08, 
               label.offset = 0.5, density.col = "#00000050", 
               density.border.col = "#00000080", no.margin = TRUE, edge.width = 4,
               tip.color = tip_colors, font = 2)

# Color the branches as per the nuclear clades:
# green:  #09af02
# violet: #b902cc
# yellow: #c9c913

branch_colors = c("black", "black", "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", 
                  "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", 
                  "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", "#ecfe36", 
                  "#ecfe36", "#b902cc", "#b902cc", "#b902cc", "#b902cc", "#b902cc", 
                  "#b902cc", "#b902cc", "#b902cc", "#b902cc", "#09af02", "#09af02", 
                  "#09af02", "black", "black", "black")

MCMC.tree.plot(phy = my_tree, cex.tips = 1, time.correction = 1, MCMC.chain = mcmc_nuclear, 
               scale.res = c("Epoch", "Age"), plot.type = "distributions", 
               cex.age = 1, cex.labels = 1, relative.height = 0.08, 
               label.offset = 0.5, density.col = "#00000050", 
               density.border.col = "#00000080", no.margin = TRUE, edge.width = 4,
               font = 2, col.tree = branch_colors)


