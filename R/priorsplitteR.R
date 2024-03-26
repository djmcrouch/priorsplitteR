library(compiler)
findDelta<-function(x,betaRand2,seRand2,deltaRange){ 
	numRows=nrow(x)
	if(is.null(nrow(x)))numRows=1
	if(numRows>1){
		index=rep(1:nrow(x),rep(length(betaRand2),nrow(x)))
		delta=(x[index,1]-rep(betaRand2,nrow(x)))/sqrt(x[index,2]^2+rep(seRand2^2,nrow(x)))
		}else{
			delta=(x[1]-betaRand2)/sqrt(x[2]^2+seRand2^2)
			}
	bool=!between(delta,deltaRange[1],deltaRange[2]) #Squash outliers
	delta[bool]=deltaRange[1*(sign(delta[bool])>0)+1]
	return(delta)
	}

make_poly<-function(x,y){ #Function to make polynomials
	y^x
		}


make_lambdaDeriv2m<-function(x,y){ #y should be a matrix, even if it just has one row
  numCols=nrow(y)
  matrix(rep(x,numCols),nrow=numCols,byrow=TRUE)*y
}

findDelta<-cmpfun(findDelta)
make_poly<-cmpfun(make_poly)
make_lambdaDeriv2m<-cmpfun(make_lambdaDeriv2m)


#' priorityFDR.priorsplitteR() function by Daniel Crouch
#' Estimates priorityFDRs from a set of effect size estimates and standard errors. Effect size estimates are assumed to have normal sampling error.
#' Returns local and tail area priorityFDRs, effect priorities and FDRs. 
#' See Crouch et al. 2022, https://doi.org/10.1101/2021.02.05.429962, for details.
#' @param beta Vector of effect size estimates with normal sampling error. Must be supplied 
#' @param s Vector of standard error estimates, assumed to be close to their true values. Must be supplied
#' @param seed Seed for random number generation. Must be supplied
#' @param excludeFromFit Vector of integers specifying the variables to exclude from distribution modelling
#' @param simulationTruth Vector of true effect sizes if known from a simulation, for plotting
#' @param outFileStem File name stem for output files. Must be supplied
#' @param breakSize Width of bin, on the Z-score scale, in histogram for Poisson regression distribution modelling. Default 0.5. 
#' @param degree_altNull Degree of polynomial model for ANE model, default 10
#' @param degree_doublePolynomial Degree of polynomial-like model for NUPE and NUPE-diff models. Default is to use half the number of parameters used to fit the likelihood for Tweedie estimation (ceiling value if this is not round)
#' @param nSim Number of simulated z-scores to resample for NUPE-diff models, default 10 times the number of variables
#' @param diffSampleSize Number of 'k' variables to sample for estimating pairwise posterior probabilities for SNP k exceeding focal SNP j in effect size, default 1000
#' @param diffSampleSizeHighAcc Number of additional 'k' variables to sample for estimating pairwise posterior probabilities for SNP k exceeding focal SNP j in effect size, when estimated local effectPriority_j<highAccThresh, default 4000 giving 4000+500=50000 total samples
#' @param highAccThresh Threshold for re-sampling variables to improve estimates of priorityfdr_j based on whether it is low enough to be potentially interesting, default 0.02
#' @param blockSize Number of 'j' variables to vectorise priorityfdr_j for, default 1000. When more memory is available, blockSize can be increased with potential improvements in speed
#' @param EMtol Vector of length 2 containing tolerances for likelihood convergence in EM algorithm (ANE and NUPE/NUPE-diff respectively), default c(1e-8,1e-3)
#' @param nStop Number of consecutive EM iterations that must satisfy EM likelihood convergence tolerance before convergence is accepted, default 5
#' @param maxEMit Vector of length 2 containing maximum number of EM iterations before convergence is accepted, for ANE and NUPE/NUPE-diff respectively, default c(100000,10000). Estimates obtained when the maximum number of iterations are reached should be treated with caution
#' @param ml.reltol Maximum likelihood tolerance, passed to constrOptim, default 1e-300
#' @param outlierSD Variables with absolute z-scores greater than this threshold are removed from distribution modelling. For producing priorityFDR estimates, they are squashed to have absolute values no greater than the maximum value used for modelling. Defaults to either 12, or the 99.9 percent quantile of absolute Z-scores (where n is the number of variables), if the latter is larger, with a default maximum of 20 imposed. Modelling the likelihoods when there are values beyond 20 becomes difficult 
#' @param diffOutlierSD Pairwise effect estimate differences with absolute z-scores greater than this threshold are removed from distribution modelling. For producing priorityFDR estimates, they are squashed to have absolute values no greater than the maximum value used for modelling. Default NULL, which sets the threshold to either outlierSD or the 99.9 percent quantile of absolute deltas, whichever is larger.
#' @param truncatedNull Range within which observed estimates are known a priori to exist, a vector of length 2. For example, if it is known that estimates can only be positive, specify c(0,Inf). Default NULL for range truncation.
#' @param threads Number of threads for parallelising Bayesian computations for pairwise variable effect size comparisons, default 1. We recommended higher settings when multiple threads are available
#' @param modelOnly If TRUE, only model estimation is performed, returning estimated parameters. PriorityFDRs, FDRs and effect priorities are not computed for each variant
#' @param performOneTail If true, one-tail inclusive priorityFDRs are computed, estimating the probability that a variable's true effect, rather than effect size, is exceeded by a randomly chosen variant. This requires additional computational time but the results can be used to estimate the prior distribution of effects using priorEst.priorsplitteR

#' @export

priorityFDR.priorsplitteR<-function(beta,se,seed,excludeFromFit=NULL,simulationTruth=NULL,outFileStem,breakSize=NULL,degree_altNull=12,degree_doublePolynomial=NULL,degree_all=12,nSim=10000000,diffSampleSize=1000,diffSampleSizeHighAcc=9000,highAccThresh=0.05,blockSize=NULL,EMtol=c(5e-4,5e-4),nStop=3,maxEMit=c(100000,1000),ml.reltol=1e-300,outlierSD=NULL,diffOutlierSD=NULL,truncatedNull=NULL,threads=1,modelOnly=FALSE,performOneTail=FALSE){
library(distr);library(parallel);library(dplyr);library(alabama);library(compiler);library(foreach)
EMtol[2]=EMtol[2]*nSim/100000
if(is.character(try(.Random.seed,silent=T)[1])){
	 runif(1) #Let R start a seed if it doesn't exist yet
	 oldSeed=.Random.seed 
		}else{ oldSeed=.Random.seed } #Store the current seed to return at the end of the function, as the user may be relying on this outside of the function. If there is no current seed, let R start one to return at the end of the function, in a way that doesn't depend on the user-specified seed variable
on.exit( .Random.seed<<-oldSeed )
set.seed(seed)
outputLog=paste("outFileStem=",outFileStem,"::breakSize=",breakSize,"::nSim=",nSim,"::diffSampleSize=",diffSampleSize,"::degree_altNull=",degree_altNull,"::degree_doublePolynomial=",degree_doublePolynomial,"::blockSize=",blockSize,"::EMtol=",EMtol,"::nStop=",nStop,"::ml.reltol=",ml.reltol,"::outlierSD=",outlierSD,"::threads=",threads)
if((nchar(outFileStem)>0)&(substr(outFileStem,nchar(outFileStem),nchar(outFileStem))!="_"))outFileStem=paste0(outFileStem,"_")
z=beta/se ; nSNP=length(beta) 
if(sum(is.na(z))>0){
	message="Error: NAs in z scores (where z=beta/s). Check for NAs or Inf values in beta or s and remove or replace as appropriate"
	print(message);return(message)
	}
if(sd(z)<1){
	message="Error: variance of z-scores less than 1 indicating poorly calibrated p-values. PriorsplitteR analysis not possible is not possible under these circumstances. Possible solutions: a) Did you use 1-tail tests? If effects go in the unexpected direction this can cause deflated Z-scores - switch to 2-tail tests. b) Try removing statistics for binary variables with very small numbers of observations in either causes or controls (e.g. rare genetic variants)."
	print(message);return(message)
	}
if(is.null(outlierSD))outlierSD=min(c(max(c(12,quantile(abs(z),0.9999))),20))
fullData=data.frame(beta=beta,se=se,z=z) #Store full data before removing outliers

if(is.null(breakSize)){
	#breakSize=3.49*sd(fullData$z)*(nrow(fullData)^(-1/3)) #This is Scott's rule for choosing histogram bin widths
	breakSize=2*IQR(fullData$z)*(nrow(fullData)^(-1/3)) #This is Scott's rule for choosing histogram bin widths
	nBin=2*outlierSD/breakSize
	if(nBin>100){  #Limit to 100 within the 'non-outlier' region for computational tractability
		breakSize=2*outlierSD/100
		}}
increment=0
#repeat{
	z=fullData[,'z'];beta=fullData[,'beta'];se=fullData[,'se'];nSNP=nrow(fullData)
	indexNonOutliers=(1:nSNP)[(abs(z)<=outlierSD)]
	indexNonOutliers_nonExclude=indexNonOutliers[!indexNonOutliers%in%excludeFromFit]
	outputLog=c(outputLog,paste0("There are ",nSNP-length(indexNonOutliers)," out of ",nSNP," variables exceeding the outlier Z-score threshold of ",outlierSD," or near-zero null densities. Removing these and proceeding with ",length(indexNonOutliers)," variables."))
	z=z[indexNonOutliers_nonExclude];beta=beta[indexNonOutliers_nonExclude];se=se[indexNonOutliers_nonExclude];nSNP=length(indexNonOutliers_nonExclude);simulationTruth=simulationTruth[indexNonOutliers_nonExclude]
	by=round(breakSize,digits=2)
	if(sum(z<0)>0){negBins=seq(-by/2,min(z)-by,by=-by)}else{negBins=NULL}
	if(sum(z>0)>0){posBins=seq(by/2,max(z)+by,by=by)}else{posBins=NULL}
	if((sum(z<0)<=0)&(sum(z>0)>0))posBins=seq(min(z)-by/2,max(z)+by,by=by)
	if((sum(z<0)>0)&(sum(z>0)<=0))negBins=seq(max(z)+by/2,max(z)+by,by=by)
	breakPoints=sort(c(negBins,posBins));breakPoints_store=breakPoints
	histogram<-hist(z,breaks=breakPoints,plot=FALSE);histogram_store=histogram #Call histogram but don't plot
	midPoints<-histogram$mids
	#if(sum((histogram$counts/sum(histogram$counts))>0.2)>0){breaks=breaks+1}else{break} #Increase number of breaks, or break the loop
		#}
#Fit model for tweedie estimation
zTweed<-(fullData$beta/fullData$se)
by=round(breakSize,digits=2)
if(sum(zTweed<0)>0){negBins=seq(-by/2,min(zTweed)-by,by=-by)}else{negBins=NULL}
if(sum(zTweed>0)>0){posBins=seq(by/2,max(zTweed)+by,by=by)}else{posBins=NULL}
if((sum(zTweed<0)<=0)&(sum(zTweed>0)>0))posBins=seq(min(zTweed)-by/2,max(zTweed)+by,by=by)
if((sum(zTweed<0)>0)&(sum(zTweed>0)<=0))negBins=seq(max(zTweed)+by/2,max(zTweed)+by,by=by)
h<-hist(zTweed,breaks=sort(c(negBins,posBins)),plot=FALSE)
AIC=rep(Inf,length=degree_all-1)
options(warn=-1)
for(d in 2:degree_all){
	polyReg=try(glm(y~poly(x,d,raw=TRUE),family=poisson(link="log"),data=data.frame(y=h$counts,x=h$mids)),silent=T)
	if(!is.character(polyReg[[1]])){
	  AIC[d-1]=2*(d+1)-2*logLik(polyReg)
	}}
d=(2:degree_all)[order(AIC,decreasing=FALSE)[1]]
polyReg=glm(y~poly(x,d,raw=TRUE),family=poisson(link="log"),data=data.frame(y=h$counts,x=h$mids))
options(warn=0)
coeffs=coef(polyReg)[-1];coeffs[is.na(coeffs)]=0  #Missing coefficients can't be estimated.
derivReg=(((poly(zTweed,degree=d,raw=TRUE)/matrix(rep(zTweed,d),ncol=d,byrow=FALSE))%*%(coeffs*(1:d)))[,1])
derivReg[(sign(derivReg)==sign(fullData$beta))&(((zTweed)<sort(h$mids,decreasing=FALSE)[2])|((zTweed)>sort(h$mids,decreasing=TRUE)[2]))]=0
post=(zTweed+derivReg)*fullData$se #Tweedie estimate of posterior expectation
zNoise=rnorm(length(zTweed),0,1)
derivRegNoise=(((poly(zNoise,degree=d,raw=TRUE)/matrix(rep(zNoise,d),ncol=d,byrow=FALSE))%*%(coeffs*(1:d)))[,1])
derivRegNoise[(sign(derivRegNoise)==sign(zNoise))&(((zNoise)<sort(h$mids,decreasing=FALSE)[2])|((zNoise)>sort(h$mids,decreasing=TRUE)[2]))]=0
postNoise=(zNoise+derivRegNoise)
postse=sd(postNoise)
png(paste0(outFileStem,"likelihood_for_Tweedie.png"),width=3000,height=3000,res=600)
plot(h,main=NULL)
x=seq(min(h$mids),max(h$mids),by=diff(range(h$mids))/1000)
options(warn=-1)
points(x=x,y=predict(polyReg,type="response",newdata=data.frame(x=x)),type="l",col=2)
options(warn=0)
dev.off()
if(is.null(degree_doublePolynomial))degree_doublePolynomial=ceiling(d/2)
histogram<-hist(z,breaks=breakPoints,plot=FALSE);histogram_store=histogram #Call histogram but don't plot
midPoints<-histogram$mids
counts<-histogram$counts
width=midPoints[2]-midPoints[1]
#null=pnorm(breakPoints)[-1]-pnorm(breakPoints)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
#bool=(midPoints>0)
#null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],lower.tail=FALSE)-pnorm(breakPoints[-1][bool],lower.tail=FALSE) 
null=dnorm(midPoints)*width
if(!is.null(truncatedNull))null=null/diff(pnorm(truncatedNull))
lNull=log(null) #Helps avoid numerical issues in function below
#Use EM algorithm to fit two group mixture density
likelihood.poissonRegMixture_altNull<-function(x){
	s=x[1];b=x[-1]
 	p=0.5*exp(s)*(s<0)+(1-0.5*exp(-s))*(s>=0)
	poly=lambdaDeriv1%*%b
	lambda=exp(poly)[,1]
	l1=counts*(1-m)*(log(1-p)+poly)-(1-p)*lambda*nSNP
	l2=counts*m*(log(p))-null*p*nSNP
	l1[(l1==(-Inf))|is.na(l1)]=-1e300
	l2[(l2==(-Inf))|is.na(l2)]=-1e300
	l=sum(l1+l2)
	return(-l)
	}
grad.poissonRegMixture_altNull<-function(x){
	s=x[1];b=x[-1]
 	p=0.5*exp(s)*(s<0)+(1-0.5*exp(-s))*(s>=0)
	lambda=exp(lambdaDeriv1%*%b)[,1]
	deriv=c(sum(lambda*nSNP-counts*(1-m)/(1-p)+counts*m/p-null*nSNP),(counts*(1-m))%*%lambdaDeriv1-(lambda*nSNP*(1-p))%*%lambdaDeriv1)
	deriv[1]=0.5*deriv[1]*(exp(s)*(s<0)+exp(-s)*(s>=0))
	return(-deriv)
	}
degreeMax<-degree_altNull
c<-1
AICstore<-vector(length=degreeMax)
parStore=NULL
for(degree in 1:degreeMax){ #Loop over degrees
	m=rep(0.5,length(counts))   #Distribution membership variables. 
	i=1;breakCounter=0;par=c(0,log(dnorm(0))-1e-5,rep(0,degree))
	lambdaDeriv1=cbind(1,matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree))
	scaleMeans=rep(0,ncol(lambdaDeriv1));scaleSDs=apply(X=lambdaDeriv1,FUN=sd,MARGIN=2)
	scaleList=list(ANE=scaleSDs)
	lambdaDeriv1[,-1]=scale(lambdaDeriv1[,-1],center=FALSE,scale=scaleSDs[-1])
	colnames(lambdaDeriv1)=NULL
 repeat{ #EM loop
	print(paste0("Alt-Null EM iteration=",i,"::likelihood polynomial degree = ",degree))
	ui=matrix(0,nrow=1+nrow(lambdaDeriv1),ncol=degree+2);ci=rep(-log(dnorm(0))+1e-10,1+nrow(lambdaDeriv1)) #All linear combinations Less than log null minus small constant
	ui[1,2]=-1
	ui[-1,2:(degree+2)]=-lambdaDeriv1
repeat{ 
      #Very occassionally constrOptim will fail at certain initial values of par. This loop catches errors and restarts from a new initialisation point after adding some random noise
mixFit=try(constrOptim(theta=par,f=likelihood.poissonRegMixture_altNull,grad=grad.poissonRegMixture_altNull,method="BFGS",control=list(maxit=1000,reltol=ml.reltol),ui=ui,ci=ci),silent=TRUE)
      if(!is.character(mixFit[1])){
	if(likelihood.poissonRegMixture_altNull(par)<mixFit$value){
		mixFit$value=likelihood.poissonRegMixture_altNull(par)
		mixFit$par=par
		}
	break
		}
      par=par+rnorm(length(par),sd=0.01)
      par[2][par[2]>(lNull[midPoints==0]-1e-10)]=lNull[midPoints==0]-1e-10 
      print("Error in likelihood maximisation. Adding random noise to starting parameters and retrying")
    }
	x=par
	s=x[1];b=x[-1]
 	p=0.5*exp(s)*(s<0)+(1-0.5*exp(-s))*(s>=0)
	poly=lambdaDeriv1%*%b
	lambda=exp(poly)[,1]
	l1=counts*(1-m)*(log(1-p)+poly)-(1-p)*lambda*nSNP
	l2=counts*m*(log(p))-null*p*nSNP
	l1[(l1==(-Inf))|is.na(l1)]=-1e300
	l2[(l2==(-Inf))|is.na(l2)]=-1e300
	Qjj=sum(l1+l2) #Don't need to worry about constants in the likelihoods - see tibshirani, hastie and friedman 'elements of statistical learning' EM algorithm section p276/277, equation 8.46 and following short paragraphs of text. The goal is find a value s.t. Q(\theta^{j+1},\theta^j) > Q(\theta^{j},\theta^j), and the denominator is the same in both these terms.
	x=mixFit$par #Find core expected likelihood given new parameters
	s=x[1];b=x[-1]
 	p=0.5*exp(s)*(s<0)+(1-0.5*exp(-s))*(s>=0)
	poly=lambdaDeriv1%*%b
	lambda=exp(poly)[,1]
	l1=counts*(1-m)*(log(1-p)+poly)-(1-p)*lambda*nSNP
	l2=counts*m*(log(p))-null*p*nSNP
	l1[(l1==(-Inf))|is.na(l1)]=-1e300
	l2[(l2==(-Inf))|is.na(l2)]=-1e300
	QjjPlus1=sum(l1+l2) 
	m<-null*p/(null*p+lambda*(1-p)) #Update missing variables
	m[(null*p+lambda*(1-p))==0]=0  #Fix numerical problems with zero denominators	
	par=mixFit$par
	if(((QjjPlus1-Qjj)<0)){message="Warning: M step has failed to improve on previous optimum"}else{message=""}
	print(paste0("Core expected likelihood update = ",QjjPlus1,"-",Qjj,"=",QjjPlus1-Qjj," ",message))
	if(((QjjPlus1-Qjj)<EMtol[1])&((QjjPlus1-Qjj)>=0)){breakCounter=breakCounter+1}else{breakCounter=0}
	if(breakCounter==nStop)break
	if(i>=maxEMit[1]){
		message=paste0("Warning: Maximum iterations (",maxEMit[1],") reached")
		outputLog=c(outputLog,message)
		print(message)
		break
			}
	i=i+1
		} #End of repeat loops
	     #Compute AIC using full likelihood
	     lFull=counts*(log((1-p)*lambda*nSNP+p*null*nSNP))-(1-p)*lambda*nSNP-null*p*nSNP-lfactorial(counts)
	     AIC=2*(degree+2)-2*(sum(lFull))
	     AICstore[c]=AIC;parStore=c(parStore,list(par))
	     c=c+1
	} #End of loop over degree
degree=(1:degreeMax)[order(AICstore,decreasing=FALSE)[1]];par=parStore[[order(AICstore,decreasing=FALSE)[1]]]
degree_doublePolynomial=max(c(2,ceiling(degree/2)))
lambdaDeriv1=cbind(1,matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree))
scaleMeans=rep(0,ncol(lambdaDeriv1));scaleSDs=apply(X=lambdaDeriv1,FUN=sd,MARGIN=2)
scaleList=list(ANE=scaleSDs)
lambdaDeriv1[,-1]=scale(lambdaDeriv1[,-1],center=FALSE,scale=scaleSDs[-1])
colnames(lambdaDeriv1)=NULL
predict.poissonRegMixture_altNull<-function(x){
	s=x[1];b=x[-1]
 	p=0.5*exp(s)*(s<0)+(1-0.5*exp(-s))*(s>=0)
	lambda=exp(lambdaDeriv1%*%b)
	lambdaMix=lambda*(1-p)+null*p
	return(nSNP*cbind(lambdaMix,lambda,null,lambda*(1-p),null*p))
		}
s=par[1];b=par[-1]
pNull=0.5*exp(s)*(s<0)+(1-0.5*exp(-s))*(s>=0)
fit<-predict.poissonRegMixture_altNull(par);fit[is.na(fit)]=0
#pNull=1-(1-pNull)*sum(lambda) #Adjust pNull to correct integrations slightly exceeding 1 
#Check fit - add more x 'midpoints' points to give a continuous appearance
x=seq(min(midPoints),max(midPoints),by=diff(range(z))/10000)
lambdaDeriv1=matrix(sapply(X=1:degree,y=x,FUN=make_poly),ncol=degree)
lambdaDeriv1=scale(lambdaDeriv1,center=scaleMeans[-1],scale=scaleSDs[-1])
lambdaDeriv1=cbind(1,lambdaDeriv1)
null=dnorm(x)*width
if(!is.null(truncatedNull))null=null/diff(pnorm(truncatedNull))
fit<-predict.poissonRegMixture_altNull(par);fit[is.na(fit)]=0
print("Printing graph of fitted alt/null mixture distribution")
#Put error bars on this fit
png(paste0(outFileStem,"poissonMixtureFit_altNull.png"),width=3000,height=3000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
plot(histogram,main=NULL,ylim=c(0,min(c(max(fit[,1:3]),nSNP)))*1.05,ylab="Count")
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,2],col=3,type="l",lwd=2)
points(x,fit[,3],col=4,type="l",lwd=2)
dev.off()
png(paste0(outFileStem,"poissonMixtureFit_altNull_mixingProp.png"),width=3000,height=3000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
plot(histogram,main=NULL,ylim=c(0,min(c(max(fit[,c(1,4,5)]),nSNP)))*1.05,ylab="Count")
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,4],col=3,type="l",lwd=2)
points(x,fit[,5],col=4,type="l",lwd=2)
dev.off()
#Now simulate from the alternative distribution
par_ANE=par
altNull_Density<-function(x){
	lambdaDeriv1=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
	lambdaDeriv1=scale(lambdaDeriv1,center=scaleMeans[-1],scale=scaleSDs[-1])
	lambdaDeriv1=cbind(1,matrix(lambdaDeriv1,nrow=nrow(lambdaDeriv1)))
	lambda=exp(lambdaDeriv1%*%b)/width
	lambda[is.na(lambda)]=0 #Infinities can occasionally cause problems
	return(as.vector(lambda))
	}
dist <-AbscontDistribution(d=altNull_Density,low1=min(midPoints),up1=max(midPoints),withStand=FALSE) #signature for a dist with pdf ~ twoGroupDensity
rdist <- r(dist) # function to create random variates 
fdr=dnorm(fullData$z)*pNull/(altNull_Density(fullData$z)*(1-pNull)+dnorm(fullData$z)*pNull)
fullData=data.frame(fullData,fdr=fdr,altDensity=altNull_Density(fullData[,'z']))
fullData[-indexNonOutliers,'fdr']=0   #Values outside range of the model usually need replacing with fdr=0
print("Sampling from the estimated distribution function")
fdr=fullData[indexNonOutliers_nonExclude,'fdr'];altDensity=fullData[indexNonOutliers_nonExclude,'altDensity']
#Following few lines from 'https://stackoverflow.com/questions/23570952/simulate-from-an-arbitrary-continuous-probability-distribution'
#Simulate differences by drawing nSim further variables.
zSim=rdist(nSim) 
#zSim=z[sample(1:length(z),size=nSim,replace=T,prob=1-fdr)]
png(paste0(outFileStem,"poissonMixtureFit_altNull_altSim.png"),width=3000,height=6000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),mfrow=c(2,1),cex.axis=0.7,cex.lab=0.7)
x<-seq(min(midPoints),max(midPoints),by=0.01)
y<-unlist(sapply(X=x,FUN=altNull_Density))
denom=mean(y)*diff(range(x))
y=y/denom
yNull=dnorm(x)
hist(zSim,freq=FALSE,breaks=breakPoints,col=0,ylim=c(0,max(c(y,yNull))),main=NULL,xlab="",xlim=range(zSim),ylab="Count")
lines(x,y,type='l',col="green")
lines(x,yNull,type='l',col="blue")
hist(zSim,freq=FALSE,breaks=100,col=0,ylim=c(0,max(c(y,yNull))),main=NULL,xlab="",xlim=range(zSim),ylab="Count")
lines(x,y,type='l',col="green")
lines(x,yNull,type='l',col="blue")
dev.off()
#Perform NNPE estimation for finding the false sign probabilities. This is the same as the NNPED model above but using SNP effects rather than simulated effect differences
#Define histogram
degree<-degree_doublePolynomial
sDelta=1
histogram<-hist(zSim,breaks=breakPoints,plot=FALSE);histogram_store_z=histogram #Call histogram but don't plot
counts<-histogram$counts
lambdaDeriv2=matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree)
scaleMeans=rep(0,ncol(lambdaDeriv2));scaleSDs=apply(lambdaDeriv2,FUN=sd,MAR=2)
scaleList=append(scaleList,list(NUPE=scaleSDs))
lambdaDeriv2=scale(lambdaDeriv2,center=FALSE,scale=scaleSDs) #Don't center first term
colnames(lambdaDeriv2)=NULL
lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)))
lambdaDeriv2m=lambdaDeriv2m/ij
ij=matrix(ij[1,],nrow=degree+1,byrow=FALSE) #For gradient function
width=midPoints[2]-midPoints[1]
nullFactor=dnorm(0,mean=0,sd=sDelta)*width
null=pnorm(breakPoints,mean=0,sd=sDelta)[-1]-pnorm(breakPoints,mean=0,sd=sDelta)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
#bool=(midPoints>0)
#null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],mean=0,sd=sDelta,lower.tail=FALSE)-pnorm(breakPoints[-1][bool],mean=0,sd=sDelta,lower.tail=FALSE) 
null=dnorm(midPoints,mean=0,sd=sDelta)*width
EM<-function(degree,counts,lambdaDeriv2,lambdaDeriv2m,ij,nullFactor,null,midPoints,width,printMessage=NULL){
  #Define functions for ML optimisation of NNPED and NNPE models
  likelihood.poissonRegMixture_twoPolynomial<-function(x){
    s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    p3=1-p1-p2;p3[p3<0]=0
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
    l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*n
    l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*n
    l3=counts*(1-m1-m2)*log(p3)-p3*null*n
    l1[(l1==(-Inf))|is.na(l1)]=-1e300
    l2[(l2==(-Inf))|is.na(l2)]=-1e300
    l3[(l3==(-Inf))|is.na(l3)]=-1e300
    l=sum(l1+l2+l3)
    return(-l)
  }
  grad.poissonRegMixture_twoPolynomial<-function(x){
    s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    p3=1-p1-p2;p3[p3<0]=0
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
    cbindlambdaDeriv2=cbind(1,lambdaDeriv2)
    coeffDerivU=2*cbindlambdaDeriv2*cbindlambdaDeriv2%*%(matrix(rep(c(c1a,b1),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    coeffDerivL=2*cbindlambdaDeriv2*cbindlambdaDeriv2%*%(matrix(rep(c(c2a,b2),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    deriv=c(sum(counts*m1/p1-counts*(1-m1-m2)/p3-lambdaL*n+null*n),sum(counts*m2/p2-counts*(1-m1-m2)/p3-lambdaU*n+null*n),((counts*m2-p2*lambdaU*n)*lambdaDeriv2[,1])%*%coeffDerivU[,-1],-((counts*m1-p1*lambdaL*n)*lambdaDeriv2[,1])%*%coeffDerivL[,-1],sum(counts*m2-p2*lambdaU*n),sum(counts*m1-p1*lambdaL*n),((counts*m2-p2*lambdaU*n)*lambdaDeriv2[,1])%*%coeffDerivU[,1],-((counts*m1-p1*lambdaL*n)*lambdaDeriv2[,1])%*%coeffDerivL[,1])
    deriv[1:2]=deriv[1:2]*c(0.5*(exp(s1)*(s1<0)+exp(-s1)*(s1>=0)),0.5*(exp(s2)*(s2<0)+exp(-s2)*(s2>=0))) 
    return(-deriv)
  }
hin<-function(x){
	s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
   	 polyUa=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2
   	 polyLa=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2
	c(-s1-s2,-c1,-c2,c1a,c2a,-polyUa,-polyLa)
	}	
hin.jac<-function(x){
    b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1a=x[2*degree+5];c2a=x[2*degree+6] 
    cbindlambdaDeriv2=cbind(1,lambdaDeriv2)
    coeffDerivU=2*cbindlambdaDeriv2*cbindlambdaDeriv2%*%(matrix(rep(c(c1a,b1),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    coeffDerivL=2*cbindlambdaDeriv2*cbindlambdaDeriv2%*%(matrix(rep(c(c2a,b2),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    jac=matrix(0,nrow=5+length(counts)*2,ncol=degree*2+6) 
    jac[1,1:2]=-1
    diag(jac[2:3,2*degree+3:4])=-1
    diag(jac[4:5,2*degree+5:6])=1
    upperIndex=6:(5+length(counts));lowerIndex=(length(counts)+6):(2*length(counts)+5)
    jac[upperIndex,2+1:degree]=-lambdaDeriv2[,1]*coeffDerivU[,-1]
    jac[lowerIndex,2+(degree+1):(2*degree)]=lambdaDeriv2[,1]*coeffDerivL[,-1]
    jac[upperIndex,2*degree+3]=-1;jac[lowerIndex,2*degree+4]=-1
    jac[upperIndex,2*degree+5]=-lambdaDeriv2[,1]*coeffDerivU[,1]
    jac[lowerIndex,2*degree+6]=lambdaDeriv2[,1]*coeffDerivL[,1]
    return(jac)
}
heq<-function(x){ #Ensures likelihoods sum to 1
    s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] 
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    p3=1-p1-p2;p3[p3<0]=0
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
    c(sum(lambdaU),sum(lambdaL))	
	}	
  likelihood.poissonRegMixture_twoPolynomial<-cmpfun(likelihood.poissonRegMixture_twoPolynomial)
  grad.poissonRegMixture_twoPolynomial<-cmpfun(grad.poissonRegMixture_twoPolynomial)
  hin<-cmpfun(hin);hin.jac<-cmpfun(hin.jac)
  lNull=log(null) #Helps avoid numerical issues in function below
  n<-sum(counts)
  m1<-rep(0.33,length(counts))   #Lower tail responsibilities 
  m2<-rep(0.33,length(counts))   #Upper tail responsibilities 
  i=1;breakCounter=0;par=c(-0.1,-0.1,rep(0,degree),rep(0,degree),-1e-5,-1e-5,1e-5,1e-5)
  repeat{
    print(paste0("Two-polynomial EM iteration (",printMessage,") =",i))
    par_tweak=par #Modify par to avoid boundary constraint issues
    par_tweak[2*degree+3:4][par[2*degree+3:4]>(-1e-2)]=-1e-2#Initialise these two parameters comfortably within the boundary region
    par_tweak[2*degree+5:6][par[2*degree+5:6]<(1e-2)]=1e-2#Initialise these two parameters comfortably within the boundary region
    repeat{ 
      #Very occassionally constrOptim will fail at certain initial values of par. This loop catches errors and restarts from a new initialisation point after adding some random noise
mixFit=try(auglag(par=par_tweak,fn=likelihood.poissonRegMixture_twoPolynomial,gr=grad.poissonRegMixture_twoPolynomial,hin=hin,hin.jac=hin.jac,heq=heq,control.optim=list(maxit=100000000,reltol=ml.reltol),control.outer=list(trace=FALSE,itmax=2,ilack.max=1,mu0=10)),silent=TRUE)
    mixFit_store=mixFit
      if(!is.character(mixFit[1])){  #Use constrOptim.nl if auglag gives a worse result than the optimisation point - this happens fairly frequently
if(likelihood.poissonRegMixture_twoPolynomial(par_tweak)<mixFit$value)mixFit=try(constrOptim.nl(par=par_tweak,fn=likelihood.poissonRegMixture_twoPolynomial,gr=grad.poissonRegMixture_twoPolynomial,hin=hin,hin.jac=hin.jac,heq=heq,control.optim=list(maxit=1000,reltol=ml.reltol),control.outer=list(trace=FALSE,itmax=2,ilack.max=1,mu0=10)),silent=TRUE)
	if(is.character(mixFit[1])){ mixFit=mixFit_store }else{ mixFit_store=mixFit }
	if(!is.character(mixFit[1]))if(likelihood.poissonRegMixture_twoPolynomial(par)<mixFit$value)mixFit=try(constrOptim.nl(par=par,fn=likelihood.poissonRegMixture_twoPolynomial,gr=grad.poissonRegMixture_twoPolynomial,hin=hin,hin.jac=hin.jac,heq=heq,control.optim=list(maxit=1000,reltol=ml.reltol),control.outer=list(trace=FALSE,mu0=10)),silent=TRUE) 
if(is.character(mixFit[1])){ mixFit=mixFit_store }else{ mixFit_store=mixFit }
	if(!is.character(mixFit[1])){
	if(likelihood.poissonRegMixture_twoPolynomial(par)<mixFit$value){
		mixFit$par=par
		mixFit$value=likelihood.poissonRegMixture_twoPolynomial(par)
			}}
	if(!is.character(mixFit[1]))break
		}	
	if(is.character(mixFit[1])){  
            ui=matrix(0,nrow=5,ncol=degree*2+6);ci=c(0,0,0,0,0)  #If it still fails try regular constrOptim
            ui[1,1:2]=c(-1,-1)  #s1 plus s2 <=0 implies p1+p2<1.
            diag(ui[2:3,degree*2+3:4])=c(-1,-1) #Intercept at 0 constraints - ratio with null should be < 1
            diag(ui[4:5,2*degree+5:6])=c(1,1) #Second intercept should be >0 to stop potential symmetry/identifiability issues
mixFit2=try(constrOptim(theta=par_tweak,f=likelihood.poissonRegMixture_twoPolynomial,grad=grad.poissonRegMixture_twoPolynomial,ui=ui,ci=ci,control=list(maxit=100000000,reltol=ml.reltol)),silent=TRUE) 
if(is.character(mixFit2[1]))mixFit2=try(constrOptim(theta=par_tweak,f=likelihood.poissonRegMixture_twoPolynomial,grad=NULL,ui=ui,ci=ci,control=list(maxit=100000000,reltol=ml.reltol)),silent=TRUE) #Try without gradient function if it fails on previous line
	    if(!is.character(mixFit2[1]))if(likelihood.poissonRegMixture_twoPolynomial(par)>mixFit2$value){
			mixFit=mixFit2
			break
				}}
      par_tweak=par_tweak+rnorm(length(par_tweak),sd=0.01)
      par_tweak[2*degree+3:4][par_tweak[2*degree+3:4]>(-1e-5)]=-1e-05 
      par_tweak[2*degree+5:6][par_tweak[2*degree+5:6]<1e-5]=1e-05 
      if(runif(1)<0.1)par_tweak=c(-0.1,-0.1,rep(0,degree),rep(0,degree),-1e-5,-1e-5,1e-5,1e-5)
      print("Error in likelihood maximisation. Adding random noise to starting parameters and retrying")
    }
    Qjj=-likelihood.poissonRegMixture_twoPolynomial(par)
    x=mixFit$par  #Find core expected likelihood with given the parameters
    s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] 
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    p3=1-p1-p2;p3[p3<0]=0
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
    l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*n
    l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*n
    l3=counts*(1-m1-m2)*log(p3)-p3*null*n
    l1[(l1==(-Inf))|is.na(l1)]=-1e300
    l2[(l2==(-Inf))|is.na(l2)]=-1e300
    l3[(l3==(-Inf))|is.na(l3)]=-1e300
    QjjPlus1=sum(l1+l2+l3)
    m1<-lambdaL*p1/(lambdaU*p2+lambdaL*p1+null*p3)
    m2<-lambdaU*p2/(lambdaU*p2+lambdaL*p1+null*p3)
    m1[(lambdaU*p2+lambdaL*p1+null*p3)==0]=0	#Fix numerical probs with zero denominators 
    m2[(lambdaU*p2+lambdaL*p1+null*p3)==0]=0
    parStore=par
    par=mixFit$par
    if(((QjjPlus1-Qjj)<0)){message="Warning: M step has failed to improve on previous optimum"}else{message=""}
    print(paste0("Core expected likelihood update = ",QjjPlus1,"-",Qjj,"=",QjjPlus1-Qjj," ",message))
    if(((QjjPlus1-Qjj)<EMtol[2])&((QjjPlus1-Qjj)>=0)){breakCounter=breakCounter+1}else{breakCounter=0}
    if(breakCounter==nStop)break
    if(i>=maxEMit[2]){
      sumLambdaU=sum(lambdaU);sumLambdaL=sum(lambdaL)
      message=paste0("Warning: Maximum iterations (",maxEMit[2],") reached")
      outputLog=c(outputLog,message)
      print(message)
      break
    }
    i=i+1
  } #End of repeat loops
  return(list(par=par,sumLambdaU=sum(lambdaU),sumLambdaL=sum(lambdaL)))
}#End of function 'EM'
EM_noNull<-function(degree,counts,lambdaDeriv2,lambdaDeriv2m,ij,nullFactor,midPoints,width,printMessage=NULL){
  #Define functions for ML optimisation of NNPED and NNPE models
  likelihood.poissonRegMixture_twoPolynomial<-function(x){
    s1=x[1];b1=x[2:(degree+1)];b2=x[(degree+2):(2*degree+1)];c1=x[2*degree+2];c2=x[2*degree+3];c1a=x[2*degree+4];c2a=x[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=1-p1
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
    l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*n
    l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*n
    l1[(l1==(-Inf))|is.na(l1)]=-1e300
    l2[(l2==(-Inf))|is.na(l2)]=-1e300
    l=sum(l1+l2)
    return(-l)
  }
  grad.poissonRegMixture_twoPolynomial<-function(x){
    s1=x[1];b1=x[2:(degree+1)];b2=x[(degree+2):(2*degree+1)];c1=x[2*degree+2];c2=x[2*degree+3];c1a=x[2*degree+4];c2a=x[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=1-p1  #0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
    cbindlambdaDeriv2=cbind(1,lambdaDeriv2)
    coeffDerivU=2*cbindlambdaDeriv2*cbindlambdaDeriv2%*%(matrix(rep(c(c1a,b1),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    coeffDerivL=2*cbindlambdaDeriv2*cbindlambdaDeriv2%*%(matrix(rep(c(c2a,b2),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    deriv=c(sum(counts*m1/p1-lambdaL*n)-sum(counts*m2/(1-p1)-lambdaU*n),((counts*m2-p2*lambdaU*n)*lambdaDeriv2[,1])%*%coeffDerivU[,-1],-((counts*m1-p1*lambdaL*n)*lambdaDeriv2[,1])%*%coeffDerivL[,-1],sum(counts*m2-p2*lambdaU*n),sum(counts*m1-p1*lambdaL*n),((counts*m2-p2*lambdaU*n)*lambdaDeriv2[,1])%*%coeffDerivU[,1],-((counts*m1-p1*lambdaL*n)*lambdaDeriv2[,1])%*%coeffDerivL[,1])
    deriv[1]=deriv[1]*0.5*(exp(s1)*(s1<0)+exp(-s1)*(s1>=0)) 
    return(-deriv)
  }
hin<-function(x){
         s1=x[1];b1=x[2:(degree+1)];b2=x[(degree+2):(2*degree+1)];c1=x[2*degree+2];c2=x[2*degree+3];c1a=x[2*degree+4];c2a=x[2*degree+5] 
   	 polyUa=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2
   	 polyLa=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2
	c(-s1,-c1,-c2,c1a,c2a,-polyUa,-polyLa)
	}	
hin.jac<-function(x){
    b1=x[2:(degree+1)];b2=x[(degree+2):(2*degree+1)];c1a=x[2*degree+4];c2a=x[2*degree+5] 
    cbindlambdaDeriv2=cbind(1,lambdaDeriv2)
    coeffDerivU=2*cbindlambdaDeriv2*cbindlambdaDeriv2%*%(matrix(rep(c(c1a,b1),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    coeffDerivL=2*cbindlambdaDeriv2*cbindlambdaDeriv2%*%(matrix(rep(c(c2a,b2),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    jac=matrix(0,nrow=5+length(counts)*2,ncol=degree*2+5) 
    jac[1,1]=-1
    diag(jac[2:3,2*degree+2:3])=-1
    diag(jac[4:5,2*degree+4:5])=1
    upperIndex=5:(4+length(counts));lowerIndex=(length(counts)+5):(2*length(counts)+4)
    jac[upperIndex,1+1:degree]=-lambdaDeriv2[,1]*coeffDerivU[,-1]
    jac[lowerIndex,1+(degree+1):(2*degree)]=lambdaDeriv2[,1]*coeffDerivL[,-1]
    jac[upperIndex,2*degree+2]=-1;jac[lowerIndex,2*degree+3]=-1
    jac[upperIndex,2*degree+4]=-lambdaDeriv2[,1]*coeffDerivU[,1]
    jac[lowerIndex,2*degree+5]=lambdaDeriv2[,1]*coeffDerivL[,1]
    return(jac)
}
heq<-function(x){ #Ensures likelihoods sum to 1
    s1=x[1];b1=x[2:(degree+1)];b2=x[(degree+2):(2*degree+1)];c1=x[2*degree+2];c2=x[2*degree+3];c1a=x[2*degree+4];c2a=x[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=1-p1  #0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
    c(sum(lambdaU),sum(lambdaL))	
	}	
  likelihood.poissonRegMixture_twoPolynomial<-cmpfun(likelihood.poissonRegMixture_twoPolynomial)
  grad.poissonRegMixture_twoPolynomial<-cmpfun(grad.poissonRegMixture_twoPolynomial)
  hin<-cmpfun(hin);hin.jac<-cmpfun(hin.jac)
  n<-sum(counts)
  m1<-rep(0.5,length(counts))   #Lower tail responsibilities 
  m2<-rep(0.5,length(counts))   #Upper tail responsibilities 
  i=1;breakCounter=0;par=c(-0.1,rep(0,degree),rep(0,degree),-1e-5,-1e-5,1e-5,1e-5)
   repeat{
    print(paste0("Two-polynomial EM iteration (",printMessage,") =",i))
    par_tweak=par #Modify par to avoid boundary constraint issues
    par_tweak[2*degree+2:3][par[2*degree+2:3]>(-1e-2)]=-1e-2#Initialise these two parameters comfortably within the boundary region
    par_tweak[2*degree+4:5][par[2*degree+4:5]<(1e-2)]=1e-2#Initialise these two parameters comfortably within the boundary region
    repeat{ 
      #Very occassionally constrOptim will fail at certain initial values of par. This loop catches errors and restarts from a new initialisation point after adding some random noise
mixFit=try(auglag(par=par_tweak,fn=likelihood.poissonRegMixture_twoPolynomial,gr=grad.poissonRegMixture_twoPolynomial,hin=hin,hin.jac=hin.jac,heq=heq,control.optim=list(maxit=100000000,reltol=ml.reltol),control.outer=list(trace=FALSE,itmax=1,ilack.max=1,mu0=10)),silent=TRUE)
      mixFit_store=mixFit
      if(!is.character(mixFit[1])){  #Use constrOptim.nl if auglag gives a worse result than the optimisation point - this happens fairly frequently
if(likelihood.poissonRegMixture_twoPolynomial(par_tweak)<mixFit$value)mixFit=try(constrOptim.nl(par=par_tweak,fn=likelihood.poissonRegMixture_twoPolynomial,gr=grad.poissonRegMixture_twoPolynomial,hin=hin,hin.jac=hin.jac,heq=heq,control.optim=list(maxit=1000,reltol=ml.reltol),control.outer=list(trace=FALSE,itmax=1,ilack.max=1,mu0=10)),silent=TRUE)
        if(is.character(mixFit[1])){ mixFit=mixFit_store }else{ mixFit_store=mixFit }
	if(!is.character(mixFit[1]))if(likelihood.poissonRegMixture_twoPolynomial(par)<mixFit$value)mixFit=try(constrOptim.nl(par=par,fn=likelihood.poissonRegMixture_twoPolynomial,gr=grad.poissonRegMixture_twoPolynomial,hin=hin,hin.jac=hin.jac,heq=heq,control.optim=list(maxit=1000,reltol=ml.reltol),control.outer=list(trace=FALSE,mu0=10)),silent=TRUE) 
       if(is.character(mixFit[1])){ mixFit=mixFit_store }else{ mixFit_store=mixFit }
	if(!is.character(mixFit[1])){
	if(likelihood.poissonRegMixture_twoPolynomial(par)<mixFit$value){
		mixFit$par=par
		mixFit$value=likelihood.poissonRegMixture_twoPolynomial(par)
			}}
	if(!is.character(mixFit[1]))break
		}	
	if(is.character(mixFit[1])){  
            ui=matrix(0,nrow=5,ncol=degree*2+5);ci=c(0,0,0,0,0)  #If it still fails try regular constrOptim
            ui[1,1]=-1
            diag(ui[2:3,degree*2+2:3])=c(-1,-1) #Intercept at 0 constraints - ratio with null should be < 1
            diag(ui[4:5,2*degree+4:5])=c(1,1) #Second intercept should be >0 to stop potential symmetry/identifiability issues
mixFit2=try(constrOptim(theta=par_tweak,f=likelihood.poissonRegMixture_twoPolynomial,grad=grad.poissonRegMixture_twoPolynomial,ui=ui,ci=ci,control=list(maxit=100000000,reltol=ml.reltol)),silent=TRUE) 
if(is.character(mixFit2[1]))mixFit2=try(constrOptim(theta=par_tweak,f=likelihood.poissonRegMixture_twoPolynomial,grad=NULL,ui=ui,ci=ci,control=list(maxit=100000000,reltol=ml.reltol)),silent=TRUE) #Try without gradient function if it fails on previous line
	    if(!is.character(mixFit2[1]))if(likelihood.poissonRegMixture_twoPolynomial(par)>mixFit2$value){
			mixFit=mixFit2
			break
				}}
      par_tweak=par_tweak+rnorm(length(par_tweak),sd=0.01)
      par_tweak[2*degree+2:3][par_tweak[2*degree+2:3]>(-1e-5)]=-1e-05 
      par_tweak[2*degree+4:5][par_tweak[2*degree+4:5]<1e-5]=1e-05 
      if(runif(1)<0.1)par_tweak=c(-0.1,rep(0,degree),rep(0,degree),-1e-5,-1e-5,1e-5,1e-5)
      print("Error in likelihood maximisation. Adding random noise to starting parameters and retrying")
    }
    Qjj=-likelihood.poissonRegMixture_twoPolynomial(par)
    x=mixFit$par  #Find core expected likelihood with given the parameters
    s1=x[1];b1=x[2:(degree+1)];b2=x[(degree+2):(2*degree+1)];c1=x[2*degree+2];c2=x[2*degree+3];c1a=x[2*degree+4];c2a=x[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=1-p1
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-midPoints^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-midPoints^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
    l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*n
    l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*n
    l1[(l1==(-Inf))|is.na(l1)]=-1e300
    l2[(l2==(-Inf))|is.na(l2)]=-1e300
    QjjPlus1=sum(l1+l2)
    m1<-lambdaL*p1/(lambdaU*p2+lambdaL*p1)
    m2<-lambdaU*p2/(lambdaU*p2+lambdaL*p1)
    m1[(lambdaU*p2+lambdaL*p1)==0]=0	#Fix numerical probs with zero denominators 
    m2[(lambdaU*p2+lambdaL*p1)==0]=0
    parStore=par
    par=mixFit$par
    if(((QjjPlus1-Qjj)<0)){message="Warning: M step has failed to improve on previous optimum"}else{message=""}
    print(paste0("Core expected likelihood update = ",QjjPlus1,"-",Qjj,"=",QjjPlus1-Qjj," ",message))
    if(((QjjPlus1-Qjj)<EMtol[2])&((QjjPlus1-Qjj)>=0)){breakCounter=breakCounter+1}else{breakCounter=0}
    if(breakCounter==nStop)break
    if(i>=maxEMit[2]){
      sumLambdaU=sum(lambdaU);sumLambdaL=sum(lambdaL)
      message=paste0("Warning: Maximum iterations (",maxEMit[2],") reached")
      outputLog=c(outputLog,message)
      print(message)
      break
    }
    i=i+1
  } #End of repeat loops
  return(list(par=par,sumLambdaU=sum(lambdaU),sumLambdaL=sum(lambdaL)))
}#End of function 'EM_noNull'
EM<-cmpfun(EM);EM_noNull<-cmpfun(EM_noNull)
EMout=EM_noNull(degree=degree_doublePolynomial,counts=counts,lambdaDeriv2=lambdaDeriv2,lambdaDeriv2m=lambdaDeriv2m,ij=ij,nullFactor=nullFactor,midPoints=midPoints,width=width,printMessage="sign mixture")
par=par_NNPE=EMout$par;sumLambdaU=EMout$sumLambdaU;sumLambdaL=EMout$sumLambdaL
sNeg=par[1]
piNeg=sumLambdaL*((0.5*exp(sNeg)*(sNeg<0)+(1-0.5*exp(-sNeg))*(sNeg>=0)))
piNeg_orig=((0.5*exp(sNeg)*(sNeg<0)+(1-0.5*exp(-sNeg))*(sNeg>=0)))
piPos=1-piNeg
#Check fit
print("Printing graph of fitted sign-probability mixture distribution")
predict.poissonRegMixture_twoPolynomial<-function(x,figStem,n,width,midPoints,scaleMeans,scaleSDs,plotNull=TRUE){
  par=x
  x=seq(min(midPoints),max(midPoints),by=diff(range(midPoints))/10000)
  lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),ncol=degree)
  lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
  lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
  lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
  ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
  nullFactor=dnorm(0,mean=0,sd=sDelta)*width
  null=dnorm(x,mean=0,sd=sDelta)*width #For visualisation purposes this is fine - multiply by width inside the following function
  if(plotNull){
	s1=par[1];s2=par[2];b1=par[3:(degree+2)];b2=par[(degree+3):(2*degree+2)];c1=par[2*degree+3];c2=par[2*degree+4];c1a=par[2*degree+5];c2a=par[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
  	p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
  	p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
  	p3=1-p1-p2;p3[p3<0]=0
		}else{
		   s1=par[1];b1=par[2:(degree+1)];b2=par[(degree+2):(2*degree+1)];c1=par[2*degree+2];c2=par[2*degree+3];c1a=par[2*degree+4];c2a=par[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
    		   p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
   		   p2=1-p1  #0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
		   p3=1 #Only used for storing the null density below
			}
  polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-x^2/2+log(nullFactor)
  polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-x^2/2+log(nullFactor)
  lambdaU=exp(polyU)[,1]
  lambdaL=exp(polyL)[,1]
  lambdaMix=lambdaL*p1+lambdaU*p2+plotNull*null*p3
  fit2=n*cbind(lambdaMix,lambdaL,lambdaU,null,lambdaL*p1,lambdaU*p2,null*p3)
  fit=fit2
  fit[is.na(fit)]=0
  png(paste0(outFileStem,figStem,".png"),width=3000,height=3000,res=600)
  par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
  plot(histogram,main=NULL,ylim=c(0,max(counts))*1.05,ylab="Count")
  points(x,fit[,1],col=2,type="l",lwd=2)
  points(x,fit[,2],col=3,type="l",lwd=2)
  points(x,fit[,3],col=5,type="l",lwd=2)
  if(plotNull)points(x,fit[,4],col=4,type="l",lwd=2)
  dev.off()
  png(paste0(outFileStem,figStem,"_mixingPropScaled.png"),width=3000,height=3000,res=600)
  par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
  plot(histogram,main=NULL,ylim=c(0,max(counts))*1.05,ylab="Count")
  points(x,fit[,1],col=2,type="l",lwd=2)
  points(x,fit[,5],col=3,type="l",lwd=2)
  points(x,fit[,6],col=5,type="l",lwd=2)
  if(plotNull)points(x,fit[,7],col=4,type="l",lwd=2)
  dev.off()
  return(fit2)
}
fit<-predict.poissonRegMixture_twoPolynomial(x=par,figStem="poissonMixtureFit_falseSign",n=sum(histogram$counts),width=width,midPoints=midPoints,scaleMeans=scaleMeans,scaleSDs=scaleSDs,plotNull=FALSE);fit[is.na(fit)]=0
fit_nupe<-fit
sign_probability<-function(x){
  s1=par_NNPE[1];b1=par_NNPE[2:(degree+1)];b2=par_NNPE[(degree+2):(2*degree+1)];c1=par_NNPE[2*degree+2];c2=par_NNPE[2*degree+3];c1a=par_NNPE[2*degree+4];c2a=x[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
  p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
  p2=1-p1
  lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
  lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
  lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
  lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
  ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
  nullFactor=dnorm(0,mean=0,sd=sDelta)
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-x^2/2+log(nullFactor)
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-x^2/2+log(nullFactor)
    lambdaU=exp(polyU)[,1]
    lambdaL=exp(polyL)[,1]
   posProb=lambdaU*p2/(lambdaL*p1+lambdaU*p2) #Ignore the null probability as, when we need to calculate signs, it is conditioned on a SNP being alt
   posProb[is.na(posProb)&(x>0)]=1 #Infinities can occasionally cause problems
   posProb[is.na(posProb)&(x<0)]=0 
  return(cbind(posProb,1-posProb,lambdaU,lambdaL))
}
posProb2=sign_probability(fullData$z)
posProb2[-indexNonOutliers,1][fullData$z[-indexNonOutliers]>0]=1#Replace values that fall outside the model range
posProb2[-indexNonOutliers,2][fullData$z[-indexNonOutliers]>0]=0
posProb2[-indexNonOutliers,1][fullData$z[-indexNonOutliers]<=0]=0
posProb2[-indexNonOutliers,2][fullData$z[-indexNonOutliers]<=0]=1
fullData=data.frame(fullData,positiveEffectProb=posProb2[,1],negativeEffectProb=posProb2[,2]) #Store for later
fullData[-indexNonOutliers,'positiveEffectProb'][fullData$z[-indexNonOutliers]>0]=1#Replace values that fall outside the model range
fullData[-indexNonOutliers,'negativeEffectProb'][fullData$z[-indexNonOutliers]>0]=0
fullData[-indexNonOutliers,'positiveEffectProb'][fullData$z[-indexNonOutliers]<=0]=0
fullData[-indexNonOutliers,'negativeEffectProb'][fullData$z[-indexNonOutliers]<=0]=1
posProb2=posProb2[indexNonOutliers_nonExclude,]
upper_Density<-function(x){
        s1=par_NNPE[1];b1=par_NNPE[2:(degree+1)];b2=par_NNPE[(degree+2):(2*degree+1)];c1=par_NNPE[2*degree+2];c2=par_NNPE[2*degree+3];c1a=par_NNPE[2*degree+4];c2a=par_NNPE[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
        p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
        p2=1-p1
        lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
        lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
        lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
        lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)

ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
        nullFactor=dnorm(0,mean=0,sd=sDelta)*width
        polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-x^2/2+log(nullFactor)
        lambdaU=exp(polyU)[,1]
	lambdaU[is.na(lambdaU)]=0 #Infinities can occasionally cause problems
	return(as.vector(lambdaU))
	}
lower_Density<-function(x){
        s1=par_NNPE[1];b1=par_NNPE[2:(degree+1)];b2=par_NNPE[(degree+2):(2*degree+1)];c1=par_NNPE[2*degree+2];c2=par_NNPE[2*degree+3];c1a=par_NNPE[2*degree+4];c2a=par_NNPE[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
        p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
        p2=1-p1
        lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
        lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
        lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
        lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)

ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
        nullFactor=dnorm(0,mean=0,sd=sDelta)*width
        polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-x^2/2+log(nullFactor)
        lambdaL=exp(polyL)[,1]
	lambdaL[is.na(lambdaL)]=0 #Infinities can occasionally cause problems
	return(as.vector(lambdaL))
	}
distU <-AbscontDistribution(d=upper_Density,low1=min(midPoints),up1=max(midPoints),withStand=FALSE) #signature for a dist with pdf ~ twoGroupDensity
rdistU <- r(distU) # function to create random variates 
distL <-AbscontDistribution(d=lower_Density,low1=min(midPoints),up1=max(midPoints),withStand=FALSE) #signature for a dist with pdf ~ twoGroupDensity
rdistL <- r(distL) # function to create random variates 
#Perform upper NNPED estimation
beta=fullData$beta;se=fullData$se;nSNP=nrow(fullData)
zTweed=rdistU(nSim)
h<-hist(zTweed,breaks=sort(c(negBins,posBins)),plot=FALSE)
if(is.null(degree_all))degree_all=length(h$mids)-1
AIC=rep(Inf,length=degree_all-1)
options(warn=-1)
for(d in 2:degree_all){
	polyReg=try(glm(y~poly(x,d,raw=TRUE),family=poisson(link="log"),data=data.frame(y=h$counts,x=h$mids)),silent=T)
	if(!is.character(polyReg[[1]])){
	  AIC[d-1]=2*(d+1)-2*logLik(polyReg)
	}}
d=(2:degree_all)[order(AIC,decreasing=FALSE)[1]]
polyReg=glm(y~poly(x,d,raw=TRUE),family=poisson(link="log"),data=data.frame(y=h$counts,x=h$mids))
options(warn=0)
coeffs=coef(polyReg)[-1];coeffs[is.na(coeffs)]=0  #Missing coefficients can't be estimated.
derivReg=(((poly(beta/se,degree=d,raw=TRUE)/matrix(rep(beta/se,d),ncol=d,byrow=FALSE))%*%(coeffs*(1:d)))[,1])
derivReg[(sign(derivReg)==sign(beta))&(((beta/se)<sort(h$mids,decreasing=FALSE)[2])|((beta/se)>sort(h$mids,decreasing=TRUE)[2]))]=0
postU=(beta/se+derivReg)*fullData$se #Tweedie estimate of posterior expectation
zNoise=rnorm(length(beta),0,1)
derivRegNoise=(((poly(zNoise,degree=d,raw=TRUE)/matrix(rep(zNoise,d),ncol=d,byrow=FALSE))%*%(coeffs*(1:d)))[,1])
derivRegNoise[(sign(derivRegNoise)==sign(zNoise))&(((zNoise)<sort(h$mids,decreasing=FALSE)[2])|((zNoise)>sort(h$mids,decreasing=TRUE)[2]))]=0
postNoise=(zNoise+derivRegNoise)
postUse=sd(postNoise)
png(paste0(outFileStem,"likelihood_for_Tweedie_upper.png"),width=3000,height=3000,res=600)
plot(h,main=NULL)
x=seq(min(h$mids),max(h$mids),by=diff(range(h$mids))/1000)
options(warn=-1)
points(x=x,y=predict(polyReg,type="response",newdata=data.frame(x=x)),type="l",col=2)
options(warn=0)
dev.off()
zTweed=rdistL(nSim)
h<-hist(zTweed,breaks=sort(c(negBins,posBins)),plot=FALSE)
if(is.null(degree_all))degree_all=length(h$mids)-1
AIC=rep(Inf,length=degree_all-1)
options(warn=-1)
for(d in 2:degree_all){
	polyReg=try(glm(y~poly(x,d,raw=TRUE),family=poisson(link="log"),data=data.frame(y=h$counts,x=h$mids)),silent=T)
	if(!is.character(polyReg[[1]])){
	  AIC[d-1]=2*(d+1)-2*logLik(polyReg)
	}}
d=(2:degree_all)[order(AIC,decreasing=FALSE)[1]]
polyReg=glm(y~poly(x,d,raw=TRUE),family=poisson(link="log"),data=data.frame(y=h$counts,x=h$mids))
options(warn=0)
coeffs=coef(polyReg)[-1];coeffs[is.na(coeffs)]=0  #Missing coefficients can't be estimated.
derivReg=(((poly(beta/se,degree=d,raw=TRUE)/matrix(rep(beta/se,d),ncol=d,byrow=FALSE))%*%(coeffs*(1:d)))[,1])
derivReg[(sign(derivReg)==sign(beta))&(((beta/se)<sort(h$mids,decreasing=FALSE)[2])|((beta/se)>sort(h$mids,decreasing=TRUE)[2]))]=0
postL=(beta/se+derivReg)*fullData$se #Tweedie estimate of posterior expectation
zNoise=rnorm(length(beta),0,1)
derivRegNoise=(((poly(zNoise,degree=d,raw=TRUE)/matrix(rep(zNoise,d),ncol=d,byrow=FALSE))%*%(coeffs*(1:d)))[,1])
derivRegNoise[(sign(derivRegNoise)==sign(zNoise))&(((zNoise)<sort(h$mids,decreasing=FALSE)[2])|((zNoise)>sort(h$mids,decreasing=TRUE)[2]))]=0
postNoise=(zNoise+derivRegNoise)
postLse=sd(postNoise)
png(paste0(outFileStem,"likelihood_for_Tweedie_lower.png"),width=3000,height=3000,res=600)
plot(h,main=NULL)
x=seq(min(h$mids),max(h$mids),by=diff(range(h$mids))/1000)
options(warn=-1)
points(x=x,y=predict(polyReg,type="response",newdata=data.frame(x=x)),type="l",col=2)
options(warn=0)
dev.off()
index1=sample(1:length(beta),size=nSim,replace=T,prob=(1-fullData$fdr)*fullData$positiveEffectProb)
index2=sample(1:length(beta),size=nSim,replace=T,prob=(1-fullData$fdr)*fullData$negativeEffectProb)
se1=se[index1];se2=se[index2]
rand=sample(1:nSim,replace=FALSE,size=nSim)
rand[rand==(1:length(rand))]=rand[rand==(1:length(rand))]+sample(c(-1,1),size=sum(rand==(1:length(rand))),replace=T)
rand[rand>nSim]=nSim;rand[rand<1]=1
deltaMatrix=data.frame(Upper_orig=(beta[index1]-beta[index1][rand])/sqrt(se1^2+se1[rand]^2),Upper_rev=(beta[index1]+beta[index2])/sqrt(se1^2+se2^2))
deltaMatrix=data.frame(deltaMatrix,Lower_orig=(beta[index2]-beta[index2][rand])/sqrt(se2^2+se2[rand]^2),Lower_rev=(beta[index2]+beta[index1])/sqrt(se2^2+se1^2))
betaMatrix=data.frame(Upper_orig=beta[index1],Upper_rev=beta[index1],Lower_orig=beta[index2],Lower_rev=beta[index2])
betaMatrix=data.frame(betaMatrix,Upper_orig_y=beta[index1][rand],Upper_rev_y=-beta[index2],Lower_orig_y=beta[index2][rand],Lower_rev_y=-beta[index1])
seMatrix=data.frame(Upper_orig=se1,Upper_rev=se1,Lower_orig=se2,Lower_rev=se2,Upper_orig_y=se1[rand],Upper_rev_y=se2,Lower_orig_y=se2[rand],Lower_rev=se1)
#Fit the 2 polynomial model, again with poisson regression. First remove any extreme outliers that exist in the deltas - should't be many
#cores<-min(c(4,parallel::detectCores(),threads))
#cluster<-parallel::makeCluster(cores,type="PSOCK")
#doParallel::registerDoParallel(cl=cluster)
#print(paste0("Estimating 4 effect size difference distributions across ",foreach::getDoParWorkers()," threads."))
outputlog=c(outputLog,paste0("Estimating 4 effect size difference distributions across ",foreach::getDoParWorkers()," threads."))
  rand1=sample(1:length(beta),size=nSim,replace=T,prob=(1-fullData$fdr))
  rand2=sample(1:length(beta),size=nSim,replace=T,prob=(1-fullData$fdr))
  delta=(beta[rand1]-beta[rand2])/sqrt(se[rand1]^2+se[rand2]^2)
  #se_r=seMatrix[,r]
  #se_ry=seMatrix[,ncol(deltaMatrix)+r]
  #beta_r=betaMatrix[,r];beta_ry=betaMatrix[,ncol(deltaMatrix)+r]
  if(is.null(diffOutlierSD)){ diffOutlierSD_r=max(c(outlierSD*1.33,quantile(abs(delta),0.9999))) }else{ diffOutlierSD_r=diffOutlierSD }
  deltaMax=diffOutlierSD_r
  bool=((abs(delta)>=deltaMax)>0)
  deltaStore=delta
  delta=delta[!bool]
  deltaRange=range(delta) #Will also use this again later in Bayes theorem computations
  if(sum(bool)>0){message=paste0("Warning: there are ",sum(bool)," extreme outlier deltas (out of total ",nSim,") with absolute values greater than than ",deltaMax," or near-zero null densities. Removing these. Consider lowering the outlierSD parameter if the number of these seems high.")}else{message=NULL}
  outputLog=c(outputLog,message);if(!is.null(message))print(message)
  sDelta=1 #sds of null differences
  by=round(breakSize,digits=2)
  if(sum(delta<0)>0){negBins=seq(-by/2,min(delta)-by,by=-by)}else{negBins=NULL}
  if(sum(delta>0)>0){posBins=seq(by/2,max(delta)+by,by=by)}else{posBins=NULL}
  breakPoints=sort(c(negBins,posBins))
  histogram<-hist(delta,breaks=breakPoints,plot=FALSE)#Call histogram but don't plot
  midPoints<-histogram$mids
  histogram<-hist(delta,breaks=breakPoints,plot=FALSE);histogram_storeDelta=histogram #Call histogram but don't plot
  midPoints<-histogram$mids
  counts<-histogram$counts
  degree<-degree_doublePolynomial
  lambdaDeriv2=matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree)
  scaleMeans=rep(0,ncol(lambdaDeriv2));scaleSDs=apply(lambdaDeriv2,FUN=sd,MAR=2)
  lambdaDeriv2=scale(lambdaDeriv2,center=FALSE,scale=scaleSDs) #Don't center 
  colnames(lambdaDeriv2)=NULL
  lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
  lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
  ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)))
  lambdaDeriv2m=lambdaDeriv2m/ij
  ij=matrix(ij[1,],nrow=degree+1,byrow=FALSE) #For gradient function
  #Use EM algorithm to fit two-polynomial constrained mixture density
  width=midPoints[2]-midPoints[1]
  nullFactor=dnorm(0,mean=0,sd=sDelta)*width
 EMout_NNPED=EM_noNull(degree=degree_doublePolynomial,counts=counts,lambdaDeriv2=lambdaDeriv2,lambdaDeriv2m=lambdaDeriv2m,ij=ij,nullFactor=nullFactor,midPoints=midPoints,width=width,printMessage="Differences")
  par_NNPED=EMout_NNPED$par
  #sumLambdaU=EMout$sumLambdaU;sumLambdaL=EMout$sumLambdaL
  #Check fit
  null=nullFactor #For visualisation purposes this is fine - multiply by width inside the following function
  print("Printing graph of fitted two-polynomial mixture distribution")
  fit<-predict.poissonRegMixture_twoPolynomial(par_NNPED,figStem=paste0("poissonMixtureFit_twoPolynomial_","Difference_from_expectations"),n=nSim,width=width,midPoints=midPoints,scaleMeans=scaleMeans,scaleSDs=scaleSDs,plotNull=FALSE);fit[is.na(fit)]=0
scaleList=append(scaleList,list(NNPED=scaleSDs))
if(!modelOnly){
#For each positive beta, find the probabilities that its deltas are negative, by looking them up on the lower distribution with Bayes' theorem. Then, assuming their effect direction was flipped but the size was the same, look their delta up on the upper distribution.
print(paste0("Sampling ",diffSampleSize," effect size differences for each variant"))
#Return to using raw data for calculating the observed deltas
beta=fullData[,'beta'];se=fullData[,'se'];nSNP=nrow(fullData);fdr=fullData[,'fdr'];posProb=fullData[,'positiveEffectProb'];negProb=fullData[,'negativeEffectProb'];#rm(fullData)
#Run bayes functions above on all SNPs, block by block
threads=min(c(parallel::detectCores(),threads))
print(paste0("Using ", threads," out of a maximum ",parallel::detectCores()," threads"))
outputLog=c(outputLog,paste0("Using ", threads," out of a maximum ",parallel::detectCores()," threads"))
if(is.null(blockSize))blockSize=round(1000000000/(threads*(degree_doublePolynomial+1)^2*(diffSampleSizeHighAcc*highAccThresh+diffSampleSize*(1-highAccThresh))))
if(blockSize>=nSNP){
	blockSize=nSNP
	message="Warning: blockSize is greater or equal to the number of variables. Proceeding using one single block of all variables"
	print(message);outputLog=c(outputLog,message)
		}
snpBlocks=as.list(data.frame(matrix(1:(blockSize*floor(nSNP/blockSize)),ncol=floor(nSNP/blockSize),byrow=FALSE)))
if(max(unlist(snpBlocks))<nSNP){
	remainders=(max(unlist(snpBlocks))+1):nSNP
	snpBlocks=c(snpBlocks,list(remainders))
		}
threadPlural="parallel threads";if(threads==1)threadPlural="thread"
bayesLookup<-function(x,y,direction="positive",beta,betaRand,se,seRand,probRand=rep(1,length(betaRand)),scaleSDs,scaleMeans=rep(0,(length(y)-5)/2),deltaRange,printMessage=NULL,sDelta=1){
   	set.seed(seed+x[1])
	print(paste0("Computing ",direction," posteriors for variables ",x[1],"-",x[length(x)],", ",printMessage))
	degree=(length(y)-5)/2
   	s1=y[1];b1=y[2:(degree+1)];b2=y[(degree+2):(2*degree+1)];c1=y[2*degree+2];c2=y[2*degree+3];c1a=y[2*degree+4];c2a=y[2*degree+5] #c1a and c2a are extra intercept terms corresponding to the constants of integration, allowing the pdf to be lower than null at x=0
    	p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
	p2=1-p1
	#x is the 'block' parameter
	rand2=sample(1:length(betaRand),replace=FALSE,size=diffSampleSize)
	betaRand2=betaRand[rand2];seRand2=seRand[rand2];probRand2=probRand[rand2]
	delta<-findDelta(x=cbind(beta[x],se[x]),betaRand2=betaRand2,seRand2=seRand2,deltaRange=deltaRange)
	lambdaDeriv2=matrix(sapply(X=1:degree,y=delta,FUN=make_poly),ncol=degree)
	lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
	lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
	lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
	ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
	nullFactor=dnorm(0,mean=0,sd=sDelta)
	lambdaU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-delta^2/2+log(nullFactor) #Save memory by calling these lambda rather than 'poly'
	lambdaL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-delta^2/2+log(nullFactor)
	rm(lambdaDeriv2);rm(lambdaDeriv2m)
	null=dnorm(as.vector(delta),mean=0,sd=sDelta)
	lambdaU=exp(lambdaU)[,1]
	lambdaL=exp(lambdaL)[,1]
	lambdaU[lambdaU==Inf]=1;lambdaL[lambdaL==Inf]=1
        lambdaMix=lambdaL*p1+lambdaU*p2
	bayes=(direction=="positive")*lambdaU*p2/lambdaMix + (direction=="negative")*lambdaL*p1/lambdaMix
	bayes[is.na(bayes)&(delta>0)]=1*(direction=="positive")#Just in case there are any variables for which the overall density is zero, giving NAs (shouldn't typically be any)
	bayes[is.na(bayes)&(delta<=0)]=1*(direction=="negative")
	bayesStore=(probRand2%*%matrix(bayes,ncol=length(x),byrow=FALSE))/(diffSampleSize)
	x=x[bayesStore<highAccThresh]
	if((length(x)>0)&(diffSampleSizeHighAcc>0)){
	    rand2=sample(1:length(betaRand),replace=FALSE,size=diffSampleSizeHighAcc)
	    betaRand2=betaRand[rand2];seRand2=seRand[rand2];probRand2=probRand[rand2]
	    delta<-findDelta(x=cbind(beta[x],se[x]),betaRand2=betaRand2,seRand2=seRand2,deltaRange=deltaRange)
	    lambdaDeriv2=matrix(sapply(X=1:degree,y=delta,FUN=make_poly),ncol=degree)
	    lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
	    lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
	    lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
	    nullFactor=dnorm(0,mean=0,sd=sDelta)
	    lambdaU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-delta^2/2+log(nullFactor) #Save memory by calling these lambda rather than 'poly'
	    lambdaL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-delta^2/2+log(nullFactor)
	    rm(lambdaDeriv2);rm(lambdaDeriv2m)
	    null=dnorm(as.vector(delta),mean=0,sd=sDelta)
	    lambdaU=exp(lambdaU)[,1]
	    lambdaL=exp(lambdaL)[,1]
	    lambdaU[lambdaU==Inf]=1;lambdaL[lambdaL==Inf]=1	
            lambdaMix=lambdaL*p1+lambdaU*p2
	    bayes=(direction=="positive")*lambdaU*p2/lambdaMix + (direction=="negative")*lambdaL*p1/lambdaMix
	    bayes[is.na(bayes)&(delta>0)]=1*(direction=="positive")#Just in case there are any variables for which the overall density is zero, giving NAs (shouldn't typically be any)
	    bayes[is.na(bayes)&(delta<=0)]=1*(direction=="negative")
	    bayesStore[bayesStore<highAccThresh]=(probRand2%*%matrix(bayes,ncol=length(x),byrow=FALSE)+bayesStore[bayesStore<highAccThresh]*(diffSampleSize))/(diffSampleSize+diffSampleSizeHighAcc) 
		}
	return(bayesStore)
		}
bayesLookup=cmpfun(bayesLookup)
print(paste0("Computing posteriors for ",nSNP," variable in blocks of ",blockSize,", using ",threads," ",threadPlural))
outputLog=c(outputLog,paste0("Computing posteriors for ",nSNP*diffSampleSize," effect size differences in blocks of ",blockSize,", using ",threads," ",threadPlural))
indexU=sample(1:length(beta),size=nSim,prob=posProb*(1-fdr),replace=T)
indexL=sample(1:length(beta),size=nSim,prob=negProb*(1-fdr),replace=T)
#Randomise order in which these are computed to break up LD structure
rand=sample(1:length(beta),size=length(beta),replace=FALSE)
EP_U=unlist(mclapply(FUN=bayesLookup,X=snpBlocks,y=par_NNPED,direction="negative",beta=beta[rand],betaRand=postU[indexU],se=se[rand],seRand=postUse*se[indexU],deltaRange=deltaRange,scaleSDs=scaleList$NNPED,probRand=rep(1,length(indexU)),printMessage="assuming original effect signs",mc.cores=threads,mc.preschedule=FALSE))
EP_U_rev=unlist(mclapply(FUN=bayesLookup,X=snpBlocks,y=par_NNPED,direction="negative",beta=beta[rand],betaRand=-postL[indexL],se=se[rand],seRand=postLse*se[indexL],deltaRange=deltaRange,scaleSDs=scaleList$NNPED,probRand=rep(1,length(indexL)),printMessage="assuming flipped effect signs",mc.cores=threads,mc.preschedule=FALSE))
EP_L=unlist(mclapply(FUN=bayesLookup,X=snpBlocks,y=par_NNPED,direction="positive",beta=beta[rand],betaRand=postL[indexL],se=se[rand],seRand=postLse*se[indexL],deltaRange=deltaRange,scaleSDs=scaleList$NNPED,probRand=rep(1,length(indexL)),printMessage="assuming original effect signs",mc.cores=threads,mc.preschedule=FALSE))
EP_L_rev=unlist(mclapply(FUN=bayesLookup,X=snpBlocks,y=par_NNPED,direction="positive",beta=beta[rand],betaRand=-postU[indexU],se=se[rand],seRand=postUse*se[indexU],deltaRange=deltaRange,scaleSDs=scaleList$NNPED,probRand=rep(1,length(indexU)),printMessage="assuming flipped effect signs",mc.cores=threads,mc.preschedule=FALSE))
set.seed(seed+2) #Reset seed as we don't know which threaded call (i.e. each new seed) finished last
matchRand=match(1:length(beta),rand)
EP_U=EP_U[matchRand];EP_U_rev=EP_U_rev[matchRand];EP_L=EP_L[matchRand];EP_L_rev=EP_L_rev[matchRand]
EPs=((EP_U*piPos+EP_U_rev*piNeg)*posProb+(EP_L*piNeg+EP_L_rev*piPos)*negProb)
priorityFDRinc=fdr*(1-pNull+0.5*pNull)+EPs*(1-pNull)*(1-fdr)  
rankSamp=priorityFDRinc*length(EPs)+ 1  
priorityFDR=fdr+EPs*(1-fdr)		
if(sum(priorityFDR>1)>0)outputLog=c(outputLog,paste0("Warning: ",sum(priorityFDR>1)," variables with local priorityFDR>1, ranging ",round(range(priorityFDR[priorityFDR>1])[1],4),"-",round(range(priorityFDR[priorityFDR>1])[2],4),". Setting these to 1."))
if(sum(priorityFDRinc>1)>0)outputLog=c(outputLog,paste0("Warning: ",sum(priorityFDRinc>1)," variables with local priorityFDRinc>1, ranging ",round(range(priorityFDRinc[priorityFDRinc>1])[1],4),"-",round(range(priorityFDRinc[priorityFDRinc>1])[2],4),". Setting these to 1."))
if(performOneTail){
	EPs_oneTail=((EP_U*piPos+piNeg)*posProb+(EP_L*piNeg+piPos)*negProb)
	priorityFDRinc_oneTail=EPs_oneTail
	priorityFDRinc_oneTail[priorityFDRinc_oneTail>1]=1
	}else{
		EPs_oneTail=priorityFDRinc_oneTail=NULL
		}
priorityFDR[priorityFDR>1]=1
priorityFDRinc[priorityFDRinc>1]=1
print("Finding tail-area rates") 
findTail<-function(x){
	cuMean=cumsum(sort(x,decreasing=FALSE))/seq_along(x)
	cuMean[rank(x,ties.method="max")]
	}
tailStats=data.frame(tailpriorityFDR=findTail(priorityFDR),tailpriorityFDRinc=findTail(priorityFDRinc),tailFDR=findTail(fdr),tailEffectPriority=findTail(EPs))
outputLog=c(outputLog,"priorityFDR.priorsplitteR analysis complete.")
print("priorityFDR.priorsplitteR analysis complete.")
write.table(outputLog,file=paste0(outFileStem,"log.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
if(performOneTail){
results=data.frame(localpriorityFDR=priorityFDR,localpriorityFDRinc_oneTail=as.numeric(priorityFDRinc_oneTail),localEffectPriority=as.numeric(EPs),tailpriorityFDR=tailpriorityFDRs,tailEffectPriority=tailStats$tailEffectPriority) #Need to incorporate tail-area versions of the one tail version
	}else{
results=data.frame(localpriorityFDR=as.numeric(priorityFDR),localpriorityFDRinc=as.numeric(priorityFDRinc),localFDR=as.numeric(fdr),localEffectPriority=as.numeric(EPs),tailStats,rankings=as.numeric(rankSamp))
	}
return(list(results=results,ANEfit=par_ANE,NNPEDfit=cbind(upper=par_NNPE,upper_reverse=par_NNPE,lower=par_NNPE),NNPEfit=par_NNPE,scaleFactors=scaleList,pi=data.frame(pNull=pNull,pPos=piPos,pNeg=piNeg),outputLog=outputLog))
}else{ #End of 'if(!modelOnly)'
	outputLog=c(outputLog,"priorityFDR.priorsplitteR analysis complete.")
	print("priorityFDR.priorsplitteR analysis complete.")
	write.table(outputLog,file=paste0(outFileStem,"log.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
	return(list(ANEfit=par_ANE,NNPEDfit=cbind(upper=par_NNPE,upper_reverse=par_NNPE,lower=par_NNPE),NNPEfit=par_NNPE,scaleFactors=scaleList,pi=data.frame(pNull=pNull,pPos=piPos,pNeg=piNeg),outputLog=outputLog))
		}
			}#End of priorityFDR function


bayesLookupPair<-function(beta,s,y,deltaRange,scaleMeans,scaleSDs,nSim,sDelta=1){
	library(dplyr)
	degree=(length(y)-6)/2
	s1=y[1];s2=y[2];b1=y[3:(degree+2)];b2=y[(degree+3):(2*degree+2)];c1=y[2*degree+3];c2=y[2*degree+4];c1a=y[2*degree+5];c2a=y[2*degree+6] 
    	p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    	p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
	#x is the 'block' parameter
	delta<-findDelta(x=c(beta[1],s[1]),betaRand2=beta[2],sRand2=s[2],deltaRange=deltaRange)
	lambdaDeriv2=matrix(sapply(X=1:degree,y=delta,FUN=make_poly),ncol=degree)
	lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
	lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=t(cbind(1,lambdaDeriv2)))
	lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
	ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
	nullFactor=dnorm(0,mean=0,sd=sDelta)
	lambdaU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1-delta^2/2+log(nullFactor) #Save memory by calling these lambda rather than 'poly'
	lambdaL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2-delta^2/2+log(nullFactor) 
	rm(lambdaDeriv2);rm(lambdaDeriv2m)
	null=dnorm(as.vector(delta),mean=0,sd=sDelta)
	lambdaU=exp(lambdaU)[,1]
	lambdaL=exp(lambdaL)[,1]
	lambdaU[lambdaU==Inf]=nSim;lambdaL[lambdaL==Inf]=nSim
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	bayes=c(lambdaU*p2/lambdaMix,null*(1-p1-p2)/lambdaMix,lambdaL*p1/lambdaMix) #Null is included in the upper lookup to allow for effect priority=1
	return(bayes)
		}

#' volcanoPlot function by Daniel Crouch
#' Produces volcano plots with points coloured by priorityFDRs and FDRs, using output from priorityFDR.priorsplitteR
#' See Crouch et al. 2022, https://doi.org/10.1101/2021.02.05.429962, for details.
#' @param beta Vector of effect size estimates with normal sampling error. Must be supplied 
#' @param s Vector of standard error estimates, assumed to be close to their true values. Must be supplied
#' @param seed Seed for random number generation. Must be supplied

#' @export

volcanoPlot<-function(beta,se,label,priorityFDR.out,FDRthresh=0.01,priorityFDRthresh=0.01,priorityFDRmeasure="tailpriorityFDR",FDRmeasure="tailFDR"){
png(paste0(label,"_volcanoPlot.png"),width=2500,height=2500,res=350)
par(mar=c(4.1,4.0,0.8,1.1))
P=pchisq((beta/se)^2,df=1,lower.tail=FALSE);P[P==0]=min(c(.Machine$double.xmin,P[P>0]))
plot(beta,-log10(P),pch=19,xlab="beta",ylab="-log10 P-value",col="grey55",ylim=c(0,max(-log10(P))*1.07))
bool=(priorityFDR.out$results[,FDRmeasure]<FDRthresh)
points(beta[bool],-log10(P[bool]),pch=19,xlab="beta",ylab="-log10 P-value",col=4)
bool=(priorityFDR.out$results[,priorityFDRmeasure]<priorityFDRthresh)
points(beta[bool],-log10(P[bool]),pch=19,xlab="beta",ylab="-log10 P-value",col=2)
legend("topleft",col=c(2,4),legend=c(paste0("priorityFDR<",priorityFDRthresh),paste0("FDR<",FDRthresh," and priorityFDR>",priorityFDRthresh)),pch=19,cex=0.8)
dev.off()
	}



