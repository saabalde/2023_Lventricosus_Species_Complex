############################################
## MCMCTreeR: Prepare MCMCTree input file ##
############################################

## Set the working directory if you need to
#sedwd(")


## Load libraries
library(MCMCtreeR, quietly = TRUE, warn.conflicts = FALSE)
library(ggtree)


################################################################################
################################################################################

## Calibrations:
## Calibration 1, Africonus - Lautoconus;     node 19: 20.0 - 28 my
## Calibration 2, L ventricosus;              node 20: 10.0 - 20 my 
## Calibration 3, A. maioensis - A. crotchii; node 35:  0.5 -  2 my
## Calibration 4, Root (mtDNA);               node 33: 30.0 - 40 my

################################################################################
################################################################################

## Mitochondrial data

## Input data
# Load the tree
mtDNA_tree <- read.tree("")

# Visualize it
ggtree(mtDNA_tree, size = 1)  + geom_treescale() + geom_tiplab() + 
  geom_text(aes(label=node), hjust=-0.3)


## Estimate the parameters for each of the nodes
# Calibration 1: Skew-Normal distribution
Africonus_Lautoconus <- tipDes(mtDNA_tree, c(35))
calibration1_tree <- estimateSkewNormal(minAge = 20, maxAge = 28, monoGroups = Africonus_Lautoconus,
                                        addMode = 6, phy = mtDNA_tree, writeMCMCtree = FALSE, 
                                        minProb = 0.01, maxProb = 0, estimateScale = FALSE, 
                                        scale = 1.85, shape = -10, plot = TRUE)
calibration1_tree

# Calibration 2: Uniform distribution
# Define the node
Lventricosus = tipDes(mtDNA_tree, c(36))

# Estimate the node priors
calibration2_tree <- estimateBound(minAge = 10, maxAge = 20, minProb = 0, rightTail = 0.1,
                                   phy = calibration1_tree$apePhy, monoGroups = Lventricosus,
                                   plot = TRUE, writeMCMCtree = FALSE)
calibration2_tree

# Calibration 3: Cauchy distribution (Normal but with longer tails)
# Define the node
Africonus <- c("0385_A_crotchii", "0064_A_maioensis")

# Estimate the node priors
calibration3_tree <- estimateBound(minAge = 0.5, maxAge = 2, minProb = 0.05, rightTail = 0.05,
                                    phy = calibration2_tree$apePhy, monoGroups = Africonus, 
                                    plot = TRUE, writeMCMCtree = FALSE)
calibration3_tree


# Calibration 4: Uniform distribution at the root
calibration4_tree <- estimateBound(minAge = 30, maxAge = 40, minProb = 0.1, rightTail = 0.1,
                                   phy = calibration3_tree$apePhy, monoGroups = mtDNA_tree$tip.label, 
                                   plot = TRUE, writeMCMCtree = TRUE, MCMCtreeName = "Lventricosus_MCMCtree_mtDNA.tre")
calibration4_tree


################################################################################
################################################################################

## Nuclear data

## Input data
# Load the tree
nuclear_tree <- read.tree("")

# Visualize it
ggtree(nuclear_tree, size = 1)  + geom_treescale() + geom_tiplab() + 
  geom_text(aes(label=node), hjust=-0.3)


## Estimate the parameters for each of the nodes
# Calibration 1: Skew-Normal distribution
calibration1_tree <- estimateSkewNormal(minAge = 20, maxAge = 28, monoGroups = nuclear_tree$tip.label,
                                        addMode = 6, phy = nuclear_tree, writeMCMCtree = FALSE, 
                                        minProb = 0.01, maxProb = 0, estimateScale = FALSE, 
                                        scale = 1.85, shape = -10, plot = TRUE)
calibration1_tree

# Calibration 2: Uniform distribution
# Define the node
Lventricosus = c("CR29", "OK0273", "SZN035", "SZN034", "SZN049", "CR17", "CH01", 
                 "CR18", "CR30", "MN04", "CV1501", "genome", "CD02", "SZN054", 
                 "SZN064", "SZN047")

# Estimate the node priors
calibration2_tree <- estimateBound(minAge = 10, maxAge = 20, minProb = 0, rightTail = 0.1,
                                   phy = calibration1_tree$apePhy, monoGroups = Lventricosus,
                                   plot = TRUE, writeMCMCtree = FALSE)
calibration2_tree

# Calibration 3: Cauchy distribution (Normal but with longer tails)
# Define the node
Africonus <- c("CV0385", "CV0064")

# Estimate the node priors
calibration3_tree <- estimateBound(minAge = 0.5, maxAge = 2, minProb = 0.05, rightTail = 0.05,
                                   phy = calibration2_tree$apePhy, monoGroups = Africonus, 
                                   plot = TRUE,
                                   writeMCMCtree = TRUE, MCMCtreeName = "Lventricosus_MCMCtree_nuclear.tre")
calibration3_tree
