#  ***Ghat.R:***
Ghat is an R function to estimate G-hat from allele frequency and effect size data. Ghat.R can be loaded into R by typing: `source("https://raw.githubusercontent.com/timbeissinger/ComplexSelection/master/Function/Ghat.R")`. There are two versions of Ghat that this will load. The first, **Ghat_func()** first calculates effect sizes using RR-BLUP and then tests for selection. The second, **Ghat_func_effectsKnown()** requires externally-estimated effects. Both functions require a vector of allele frequency change as input. The **Ghat_func()** funcion also requires input of phenotypes and genotypes so that effects can be estimated, while **Ghat_func_effectsKnown()** simply expects a vector of SNP effects. All inputs are listed below:

* geno:   Matrix of genotype data. Individuals in rows, genotypes (-1, 0, 1) in columns. (required when using Ghat_func)
* phen:   Matrix of phenotype data. Individuals in rows, phenotype in single column. (required when using Ghat_func)
* effects: Vector of effect sizes. (required when using Ghat_func_effectsKnown)
* change: Change in allele frequency. Must record the change at the allele corresponding to "1" in the genotype file.
* method: "trim", "scale", or "vanilla". How should LD be accounted for?
* perms:  Number of permutations to run
* plot:   "Ghat", "Cor", or "Both", Should a plot of the Ghat or correlation test be returned?
* blockSize: How large should blocks for trimming be? Only required if method = "trim"
* num_eff: What is the effective number of independent markers?