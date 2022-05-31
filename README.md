# priorsplitteR  
R package for estimating prior densities above and below thresholds, using f-modelling. Includes BFDR estimation. Paper: https://doi.org/10.1101/2021.02.05.429962  

#Installation (In R):  
library(devtools)  
install_github("djmcrouch/priorsplitteR")

#Run (in R):  
#Beta and s should be vectors of the same length. Seed is a scalar.  
BFDR.priorsplitteR(beta=beta,s=s,seed=seed,outFileStem="outFileStem")
