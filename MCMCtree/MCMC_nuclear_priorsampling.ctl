          seed = 12345
       seqfile = ./Lventricosus_MCMCtree_Input.phy
      treefile = ./Lventricosus_MCMCtree_Input.tre
      mcmcfile = Lventricosus_priorsampling_chain1.mcmc
       outfile = Lventricosus_priorsampling_chain1.out

         ndata = 50
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 0    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates

         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5  * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0    * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = 2 4200 1  * gamma prior for overall rates for genes
  sigma2_gamma = 1 10 1  * gamma prior for sigma^2     (for clock=2 or 3)

         print = 1
        burnin = 1000000
      sampfreq = 100
       nsample = 100000
