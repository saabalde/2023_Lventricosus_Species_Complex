# Seed for the random number generator (-1 selects a random seed)
seed =  12345

# input and output files
seqfile  = ./Lventricosus.Supermatrix1.phy           # sequence file (per-locus alignments)
Imapfile = ./Lventricosus.Imap.txt                   # assignments of samples to species
outfile  = Lventricosus.A10.a0.out    # output log file
mcmcfile = Lventricosus.A10.a0.mcmc   # file to log mcmc samples

# Selection of analysis type by setting options "speciesdelimitation" and "speciestree"
# A00 - estimation of parameters on fixed phylogeny
# A10 - species delimitation using a guide tree
# A01 - species tree inference
# A11 - joint species tree inference and species delimitation

# enable species delimitation (two available algorithms):
# speciesdelimitation = 0            # species delimitation disabled (default)
speciesdelimitation = 1 0 2        # species delimitation algorithm 0 finetune (e)
# speciesdelimitation = 1 1 2 0.5    # species delimitation algorithm 1 finetune (a m)

# enable species tree inference
# speciestree = 1

# specification of: # of species, whitespace-separated list of species
# followed by max number of sequences for each species at a locus,
# and a starting (fixed for A00 and A10) species tree
species&tree = 5  green blue orange red purple
                  2     3    4      2   5
                  ((green, ((blue, orange), red)), purple);

usedata = 1      #  0: do no use data (prior); 1: use sequence data
nloci = 50       # number of data sets (alignments) in seqfile

cleandata = 0    # remove sites with ambiguity data (1:yes, 0:no)?

# species model prior (four potential priors: 0,1,2,3)
# Method A01 (species tree inference) uses speciesmodelprior = 1
# Method A10 (species delimitation with a fixed guide tree) uses either speciesmodelprior=0,1
# Method A11 (joint species delimitation and species tree inference) uses either speciesmodelprior=0,1,2,3
speciesmodelprior = 1         * 0: uniform labeled histories; 1:uniform rooted trees

# theta and gamma priors
thetaprior = 3 0.004 e       # Inverse-Gamma(a, b) for theta
tauprior   = 3 0.004         # Inverse-Gamma(a, b) for root tau
#thetaprior = gamma 2 1000     # Gamma(a,b) for theta
#tauprior   = gamma 2 1000     # Gamma(a,b) for root tau

# auto-tune step-length parameters during burnin (1: yes, 0: no)
# Potentially, followed by a colon and starting values (otherwise defaults are used)
finetune = 1
# finetune = 1: .01 .0001 .005 .0005 .2 .01 .01 .01  # GBtj, GBspr, theta, tau, mix, locusrate, ...

# binary flags on what to log
print = 1 0 0 0    # MCMC samples, locusrate, heredity scalars, Gene trees

# MCMC chain information
# Total chain length is: burnin + sampfreq*nsample
# First burnin samples are discarded, then we log every sampfreq-th sample
burnin = 10000       # discard first 8000 steps
sampfreq = 10        # log sample every 2nd step (after burnin)
nsample = 1000000    # number of samples to log in mcmcfile

# threads
threads = 4          # threads = threads
# threads = 2 1 1        # threads = threads starting_core step