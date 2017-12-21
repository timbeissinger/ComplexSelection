# ComplexSelection
This repository houses code and data to identify selection on traits controlled by many loci. The corresponding manuscript, "a simple test identifies selection on complex traits in breeding and experimentally evolved populations", is available at [insert bioarxiv link here]. In this repository, we have divided files into subdirectories as outlined below. The majority of analysis code is in the form of R-markdown scripts, for which both raw files and html output have been deposited.

## *Function:* R function to estimate G-hat from allele frequency and effect size data.
***Ghat.R:*** R function to estimate G-hat from allele frequency and effect size data.



## *Maize:* Code to analyze maize data
***Frequency_change_rrBlup_reverseInlay.Rmd:*** Code to calculate allele frequencies, estimate effects, and run analysis on maize data.

***Frequency_change_rrBlup_reverseInlay.html:*** HTML-readable code.

***ld_calc_maize_03.r:*** Code to calculate LD decay in maize data.


## *Chickens:* Code to analyze chicken data
***Frequency_change_rrBlup_reverse_inset_EW1.Rmd:*** Code to calculate allele frequencies, estimate effects, and run analysis on chicken data.

***Frequency_change_rrBlup_reverse_inset_EW1.html:*** HTML-readable code.

***Frequency_change_rrBlup_outliers_removed_EW1.Rmd:*** Code to calculate allele frequencies, estimate effects, and run analysis on chicken data, with outliers removed.

***Frequency_change_rrBlup_outliers_removed_EW1.html:*** HTML-readable code.

***ld_calc_chicken.r:*** Code to calculate LD decay in chicken data.

## *Simulations:* Code to run and analyze simulation code
***qtl_10.prm:*** 10-QTL simulation parameter file
***run_10.sh:*** 10-QTL simulation submission script for OpenLava cluster
***Analyze10.Rmd:*** Code to analyze 10-QTL simulation
***Analyze10.html:*** HTML-readable code
***Map10.Rmd:*** Code to do seleciton mapping in 10-QTL simulation
***Map10.html:*** HTML-readable code

***qtl_50.prm:*** 50-QTL simulation parameter file
***run_50.sh:*** 50-QTL simulation submission script for OpenLava cluster
***Analyze50.Rmd:*** Code to analyze 50-QTL simulation
***Analyze50.html:*** HTML-readable code
***Map50.Rmd:*** Code to do seleciton mapping in 50-QTL simulation
***Map50.html:*** HTML-readable code

***qtl_100.prm:*** 100-QTL simulation parameter file
***run_100.sh:*** 100-QTL simulation submission script for OpenLava cluster
***Analyze100.Rmd:*** Code to analyze 100-QTL simulation
***Analyze100.html:*** HTML-readable code
***Map100.Rmd:*** Code to do seleciton mapping in 100-QTL simulation
***Map100.html:*** HTML-readable code

***qtl_1000.prm:*** 1000-QTL simulation parameter file
***run_1000.sh:*** 1000-QTL simulation submission script for OpenLava cluster
***Analyze1000.Rmd:*** Code to analyze 1000-QTL simulation
***Analyze1000.html:*** HTML-readable code
***Map1000.Rmd:*** Code to do seleciton mapping in 1000-QTL simulation
***Map1000.html:*** HTML-readable code

## Data
