# priorsplitteR  
R package for estimating prior densities above and below thresholds, using f-modelling. Includes priorityFDR estimation. Paper: https://doi.org/10.1101/2021.02.05.429962  

#Installation (In R):  
library(devtools)  
install_github("djmcrouch/priorsplitteR")

#Simulated example data is loaded with the priorsplitteR package. This contains effect estimates (beta), and their standard errors (s) required for priorityFDR estimation. 
#It also contains z=beta/s, which gives chi-square P values for each variable using pchisq(z^2,df=1,lower.tail=FALSE). 
#True effects are provided, both unscaled (true_beta), and on the z scale, true_beta/s (true_z), if users wish to compare these with the estimated quantities.

library(priorsplitteR)
colnames(simulatedData) 
#"beta"  "s" "z" "true_beta" "true_z"


#Pseudo-random numbers are used during model fitting and estimation. To ensure results are fully reproducible, we require users to specify a seed value, rather than leaving R to derive it from the system clock.
seed<-12345

#Run priorityFDR estimation. Beta and s need to be vectors of the same length. A filename stem needs to be specified to outFileStem. 
priorityFDR.out<-priorityFDR.priorsplitteR(

beta=simulatedData[,'beta'],

s=simulatedData[,'s'],

seed=seed,

outFileStem="yourFileNameStem",

blockSize=2000
                 )

#The function first estimates mixtures of alt and null effects. It then decomposes the alt distribution into separate mixtures for variables with positive true effects, negative true effects, and effects with indistinguishable signs (effects too small and/or standard errors too large). It prints graphs for each of these estimated distributions alongside histograms of the data, e.g. the following alt-null mixture model for the simulated example data:
![singleSim_poissonMixtureFit_altNull_mixingProp](https://github.com/djmcrouch/priorsplitteR/assets/56267642/9a8240d3-512b-4fe6-86ff-3d49c37ff8cf)

#The following is the decomposition of the alt mixture, the green distribution in the above, into mixtures of variables with negative (green), positive (turquoise) and those with unknown effects (blue):
![singleSim_poissonMixtureFit_falseSign_mixingPropScaled](https://github.com/djmcrouch/priorsplitteR/assets/56267642/47cc5ad5-a954-4a4b-bafd-5699821f866f)

#After model estimation, the function produces estimates of effect priorities for each variable, which requires some random re-sampling. This phase may be time consuming but can be parallelised using the threads setting. By default, the number of threads is set to 1, for use on a personal computer without excessively slowing down other applications, but if multiple threads are available, e.g. on a computing cluster, this should be increased. In total, the priorityFDR function takes between around 30 minutes to 2 hours to run using 4 threads, depending on the number of variables and complexity of the effect distributions. The model estimation phase is optimally parallelised over 4 threads, with additional threads not providing an advantage. However, the following phase, which estimates effect priorities for each variable, always benefits from the addition of more threads.

#The argument blockSize determines how many variables are processed at once during the second phase. On a computing cluster, we recommended not specifying this argument so it can be set by the package. On a personal computer, it should be set sufficiently low so as to not cause memory errors in R, e.g. around 500-3000.

#The main output of interest in the returned object is results, giving local and tail-area estimates of priorityFDRs, effect priorities and FDRs:
names(priorityFDR.out$results). 

#"results"      "ANEfit"       "NNPEDfit"     "NNPEfit"      "scaleFactors"      "pi"           "outputLog" 

#The other elements are parameter estimates and other values produced during each phase of the modelling/estimation process, plus an output log. The value in priorityFDR.out$pi will be of interest to some users, containing the proportion of variables estimated to have null effects (indistinguishable from zero), plus the proportion of alt (non-null) true effects estimated to be positive and negative. Subtracting the latter two quantities from 1 gives the proportion of alt variables that have true effects that are too uncertain to be assigned a sign, despite being from the alt distribution.

#The volcanoPlot function can be used to create a quick graphing of effect sizes versus significance, labelled by priorityFDRs and FDRs, using the estimates for each variable in priorityFDR.out$results. This is the main way we visualise how significance and effect size are traded off during priorityFDR-based variable selection.
volcanoPlot(

beta=simulatedData[,'beta'],

s=simulatedData[,'s'],

label="yourFileName",

priorityFDR.out=priorityFDR.out,

)

#Open the volcano plot file:
![yourFileName_volcanoPlot](https://github.com/djmcrouch/priorsplitteR/assets/56267642/882bccfe-6e09-4b0c-9b65-8cd5bbd2dc38)


#By default, the red points are labelled using the tailpriorityFDR but this can be changed to any other variable in names(priorityFDR.out$results) using the priorityFDRmeasure argument. Likewise, the FDRmeasure argument can be used to colour points in blue - by default this is the tailFDR. The blue points are intended to indicate variables that are selected by significance alone, whereas the red points are those selected using both effect size and significance via one of the priorityFDR estimates (either tail area or local). Thresholds for the blue and red points are provided by arguments FDRthresh and priorityFDRthresh respectively, which are 0.01 by default. 










