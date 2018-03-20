# ComplexSelection
This repository houses code and data to identify selection on traits controlled by many loci. The corresponding paper, "a simple test identifies selection on complex traits", is available at [Genetics](http://www.genetics.org/content/early/2018/03/14/genetics.118.300857.long). In this repository, we have divided files into subdirectories as outlined below. For analyses performed using R-markdown, we have deposited both raw files and html outputs.

## *Function:* R function to estimate G-hat from allele frequency and effect size data.
* ***Ghat.R:*** R function to estimate G-hat from allele frequency and effect size data.



## *Maize:* Code to analyze maize data
* ***Frequency_change_rrBlup_reverseInlay.Rmd:*** Code to calculate allele frequencies, estimate effects, and run analysis on maize data.

* ***Frequency_change_rrBlup_reverseInlay.html:*** HTML-readable code.

* ***ld_calc_maize_03.r:*** Code to calculate LD decay in maize data.


## *Chickens:* Code to analyze chicken data
* ***Frequency_change_rrBlup_reverse_inset_EW1.Rmd:*** Code to calculate allele frequencies, estimate effects, and run analysis on chicken data.

* ***Frequency_change_rrBlup_reverse_inset_EW1.html:*** HTML-readable code.

* ***Frequency_change_rrBlup_outliers_removed_EW1.Rmd:*** Code to calculate allele frequencies, estimate effects, and run analysis on chicken data, with outliers removed.

* ***Frequency_change_rrBlup_outliers_removed_EW1.html:*** HTML-readable code.

* ***ld_calc_chicken.r:*** Code to calculate LD decay in chicken data.

## *Simulations:* Code to run and analyze simulation code
* ***qtl_10.prm:*** 10-QTL simulation parameter file
* ***run_10.sh:*** 10-QTL simulation submission script for OpenLava cluster
* ***Analyze10.Rmd:*** Code to analyze 10-QTL simulation
* ***Analyze10.html:*** HTML-readable code
* ***Map10.Rmd:*** Code to do seleciton mapping in 10-QTL simulation
* ***Map10.html:*** HTML-readable code
* ***Analyze10_L44.Rmd:*** Code to tune LD threshold for determining M_eff. Code determines mean LD decay over a specific distance (Specified in the 44th line of the QMSim LD output file for the 10-QTL simulation)and calculate the false positive rate. This example script was repeatedly run for lines 44-58 and used to create raw data for generating the M_eff supplemental figure in the manuscript.


* ***qtl_50.prm:*** 50-QTL simulation parameter file
* ***run_50.sh:*** 50-QTL simulation submission script for OpenLava cluster
* ***Analyze50.Rmd:*** Code to analyze 50-QTL simulation
* ***Analyze50.html:*** HTML-readable code
* ***Map50.Rmd:*** Code to do seleciton mapping in 50-QTL simulation
* ***Map50.html:*** HTML-readable code
* ***Analyze50_L44.Rmd:*** Code to tune LD threshold for determining M_eff. Code determines mean LD decay over a specific distance (Specified in the 44th line of the QMSim LD output file for the 50-QTL simulation)and calculate the false positive rate. This example script was repeatedly run for lines 44-58 and used to create raw data for generating the M_eff supplemental figure in the manuscript.

* ***qtl_100.prm:*** 100-QTL simulation parameter file
* ***run_100.sh:*** 100-QTL simulation submission script for OpenLava cluster
* ***Analyze100.Rmd:*** Code to analyze 100-QTL simulation
* ***Analyze100.html:*** HTML-readable code
* ***Map100.Rmd:*** Code to do selection mapping in 100-QTL simulation
* ***Map100.html:*** HTML-readable code
* ***Analyze100_L44.Rmd:*** Code to tune LD threshold for determining M_eff. Code determines mean LD decay over a specific distance (Specified in the 44th line of the QMSim LD output file for the 100-QTL simulation)and calculate the false positive rate. This example script was repeatedly run for lines 44-58 and used to create raw data for generating the M_eff supplemental figure in the manuscript.

* ***qtl_1000.prm:*** 1000-QTL simulation parameter file
* ***run_1000.sh:*** 1000-QTL simulation submission script for OpenLava cluster
* ***Analyze1000.Rmd:*** Code to analyze 1000-QTL simulation
* ***Analyze1000.html:*** HTML-readable code
* ***Map1000.Rmd:*** Code to do seleciton mapping in 1000-QTL simulation
* ***Map1000.html:*** HTML-readable code
* ***Analyze1000_L44.Rmd:*** Code to tune LD threshold for determining M_eff. Code determines mean LD decay over a specific distance (Specified in the 44th line of the QMSim LD output file for the 1000-QTL simulation)and calculate the false positive rate. This example script was repeatedly run for lines 44-58 and used to create raw data for generating the M_eff supplemental figure in the manuscript.

* ***qtl_10000.prm:*** 10000-QTL simulation parameter file
* ***run_10000.sh:*** 10000-QTL simulation submission script for OpenLava cluster
* ***Analyze10000.Rmd:*** Code to analyze 10000-QTL simulation
* ***Analyze10000.html:*** HTML-readable code
* ***Map10000.Rmd:*** Code to do seleciton mapping in 10000-QTL simulation
* ***Map10000.html:*** HTML-readable code

* ***LD_Thresh_Fig.R:*** R code to compile LD threshold information for generating M_eff supplemental figure and determining appropriate M_eff.
* ***PowerFigure.R:*** Code to compile power information comparing G-hat to Fst mapping and generate Figure 1 in the manuscript.

* ***intermittent.prm*** Intermittent selection parameter file.
* ***on20off5.Rmd*** Code to analyze simulation involving 20 gens with selection followed by 5 without.
* ***on20off5.html*** HTML-readable code
* ***on20off20.Rmd*** Code to analyze simulation involving 20 gens with selection followed by 20 without.
* ***on20off20.html*** HTML-readable code
* ***on20off50.Rmd*** Code to analyze simulation involving 20 gens with selection followed by 50 without.
* ***on20off50.html*** HTML-readable code
* ***on20off100.Rmd*** Code to analyze simulation involving 20 gens with selection followed by 100 without.
* ***on20off100.html*** HTML-readable code
* ***on1off5.Rmd*** Code to analyze simulation involving 1 gens with selection followed by 5 without.
* ***on1off5.html*** HTML-readable code
* ***on1off20.Rmd*** Code to analyze simulation involving 1 gens with selection followed by 20 without.
* ***on1off20.html*** HTML-readable code

* ***intensity.prm*** Parameter file for variable selection intensity.
* ***intensity01.Rmd*** Code to analyze simulation involving a selected proportion of 1 percent.
* ***intensity01.html*** HTML-readable code
* ***intensity05.Rmd*** Code to analyze simulation involving a selected proportion of 5 percent.
* ***intensity05.html*** HTML-readable code
* ***intensity5.Rmd*** Code to analyze simulation involving a selected proportion of 50 percent.
* ***intensity5.html*** HTML-readable code

* ***sampleSize.RData*** Code to analyze effect of sample size.
* ***sampleSize.html*** HTML-readable code.

* ***gens.prm*** Parameter file for variable numbers of generations
* ***Gens1.Rmd*** Code to analyze 1 generation of selection
* ***Gens1.html*** HTML-readable code
* ***Gens10.Rmd*** Code to analyze 10 generation of selection
* ***Gens10.html*** HTML-readable code
* ***Gens50.Rmd*** Code to analyze 50 generation of selection
* ***Gens50.html*** HTML-readable code
* ***Gens100.Rmd*** Code to analyze 100 generation of selection
* ***Gens100.html*** HTML-readable code

* ***phenotypeGen.prm*** Parameter file for varying phenotyping generation.
* ***Phen0.Rmd*** Code to analyze phenotyping in generation 0
* ***Phen0.html*** HTML-readable code
* ***Phen10.Rmd*** Code to analyze phenotyping in generation 10
* ***Phen10.html*** HTML-readable code

## Data
* Maize data are published in:
  + Lorenz, A. J., Beissinger, T.M., Rodrigues, R., de Leon, N. 2015. [Selection for silage yield and composition did not affect genomic diversity within the Wisconsin Quality Synthetic maize population.](http://www.g3journal.org/content/early/2015/02/02/g3.114.015263.abstract) Genes Genomes Genetics. DOI: 10.1534/g3.114.015263.

* Chicken data are available from [Figshare](https://figshare.com/articles/Brown_White_AllPhenotypes_FrequencyChange_Effects_RDS/5899267).
