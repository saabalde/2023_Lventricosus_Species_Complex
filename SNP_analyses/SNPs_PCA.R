############################################
## Plot the PCA generated using Plink 1.9 ##
############################################

## See: https://speciationgenomics.github.io/


## Set working directory
#setwd("Path_to_the_working_directory")


## Load the library
library(tidyverse)
library(extrafont)


## Load the input data
# Read the eigenvectors of the PCA results
pca <- read_table2("./Lventricosus.snps.1-PCA.eigenvec", col_names = FALSE)

# Read the eigenvalues
eigenval <- scan("./Lventricosus.snps.1-PCA.eigenval")


## Format the PCA data
# Remove the duplicated column with the individual IDs
pca <- pca[,-1]

# Change the column names from something a bit more informative
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Add a column to the PCA with the clade information (we will use ths to color
# the individuals accordingly)
clade <- rep(NA, length(pca$ind))
clade[grep("CD02", pca$ind)] <- "purple"
clade[grep("CH01", pca$ind)] <- "orange"
clade[grep("CR17", pca$ind)] <- "blue"
clade[grep("CR18", pca$ind)] <- "red"
clade[grep("CR29", pca$ind)] <- "blue"
clade[grep("CR30", pca$ind)] <- "red"
clade[grep("CV0064", pca$ind)] <- "outgroup"
clade[grep("CV0385", pca$ind)] <- "outgroup"
clade[grep("CV1501", pca$ind)] <- "purple"
clade[grep("MN04", pca$ind)] <- "purple"
clade[grep("OK0273", pca$ind)] <- "blue"
clade[grep("SZN034", pca$ind)] <- "brown"
clade[grep("SZN035", pca$ind)] <- "orange"
clade[grep("SZN047", pca$ind)] <- "green"
clade[grep("SZN049", pca$ind)] <- "orange"
clade[grep("SZN054", pca$ind)] <- "brown"
clade[grep("SZN064", pca$ind)] <- "green"

pca <- as.tibble(data.frame(pca, clade))


## Plot the eigenvalues to see how much of the variance they explain
pve <- data.frame(PC = 1:17, pve = eigenval/sum(eigenval)*100)
eigenvalues_plot <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
eigenvalues_plot + ylab("Percentage variance explained") + theme_light()

# We can easily see that PC1 and PC2 together explain more than 30% of the total
# variance. Yet, we can see the exact values:
cumsum(pve$pve)


## Plot the PCA, using for that PC1 and PC2
# plot pca
pca_plot <- ggplot(pca, aes(PC1, PC2, col = clade)) + 
            geom_point(size = 3)
pca_plot <- pca_plot + scale_colour_manual(values = c("#2000de", "#a46301", 
                                                      "#09af02", "#fba400", 
                                                      "black", "#b902cc",
                                                      "#ff0600"))
pca_plot <- pca_plot + coord_equal() + theme_light()
pca_plot + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
           ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + 
           theme(text=element_text(size=14,  family = "Arial"))

################################################################################
################################################################################

###############################################
## Sub sampling: no green clade or outgroups ##
###############################################

## Load the input data
# Read the eigenvectors of the PCA results
pca_subsampling <- read_table2("./Lventricosus.snps.3-PCA_subset.eigenvec", col_names = FALSE)

# Read the eigenvalues
eigenval_subsampling <- scan("./Lventricosus.snps.3-PCA_subset.eigenval")


## Format the PCA data
# Remove the duplicated column with the individual IDs
pca_subsampling <- pca_subsampling[,-1]

# Change the column names from something a bit more informative
names(pca_subsampling)[1] <- "ind"
names(pca_subsampling)[2:ncol(pca_subsampling)] <- paste0("PC", 1:(ncol(pca_subsampling)-1))

# Add a column to the PCA with the clade information (we will use ths to color
# the individuals accordingly)
clade <- rep(NA, length(pca_subsampling$ind))
clade[grep("CD02", pca_subsampling$ind)] <- "purple"
clade[grep("CH01", pca_subsampling$ind)] <- "orange"
clade[grep("CR17", pca_subsampling$ind)] <- "blue"
clade[grep("CR18", pca_subsampling$ind)] <- "red"
clade[grep("CR29", pca_subsampling$ind)] <- "blue"
clade[grep("CR30", pca_subsampling$ind)] <- "red"
clade[grep("CV1501", pca_subsampling$ind)] <- "purple"
clade[grep("MN04", pca_subsampling$ind)] <- "purple"
clade[grep("OK0273", pca_subsampling$ind)] <- "blue"
clade[grep("SZN034", pca_subsampling$ind)] <- "brown"
clade[grep("SZN035", pca_subsampling$ind)] <- "orange"
clade[grep("SZN049", pca_subsampling$ind)] <- "orange"
clade[grep("SZN054", pca_subsampling$ind)] <- "brown"

pca_subsampling <- as.tibble(data.frame(pca_subsampling, clade))


## Plot the eigenvalues to see how much of the variance they explain
pve_subsampling <- data.frame(PC = 1:13, pve = eigenval_subsampling/sum(eigenval_subsampling)*100)
eigenvalues_subsampling_plot <- ggplot(pve_subsampling, aes(PC, pve)) + 
                                geom_bar(stat = "identity")
eigenvalues_subsampling_plot + ylab("Percentage variance explained") + theme_light()

# We can easily see that PC1 and PC2 together explain almost 25% of the total
# variance. Yet, we can see the exact values:
cumsum(pve_subsampling$pve)


## Plot the PCA, using for that PC1 and PC2
# plot pca
pca_subsampling_plot <- ggplot(pca_subsampling, aes(PC1, PC2, col = clade)) + 
  geom_point(size = 3)
pca_subsampling_plot <- pca_subsampling_plot + 
                        scale_colour_manual(values = c("#2000de", "#a46301",
                                                      "#fba400", "#b902cc",
                                                      "#ff0600"))
pca_subsampling_plot <- pca_subsampling_plot + coord_equal() + theme_light()
pca_subsampling_plot + xlab(paste0("PC1 (", signif(pve_subsampling$pve[1], 3), "%)")) + 
                       ylab(paste0("PC2 (", signif(pve_subsampling$pve[2], 3), "%)"))

# With labels
pca_subsampling_plot + xlab(paste0("PC1 (", signif(pve_subsampling$pve[1], 3), "%)")) + 
                       ylab(paste0("PC2 (", signif(pve_subsampling$pve[2], 3), "%)")) + 
                       geom_text(aes(label = ind), hjust = 0.5,  vjust = -1, family = "Arial", size = 4) + 
                       theme(text=element_text(size=14,  family = "Arial"))
