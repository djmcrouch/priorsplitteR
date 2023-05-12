findDelta<-function(x,betaRand2,sRand2,deltaRange){ #Assumes dim(x)>0
	numRows=nrow(x)
	if(is.null(nrow(x)))numRows=1
	if(numRows>1){
		index=rep(1:nrow(x),rep(length(betaRand2),nrow(x)))
		delta=(x[index,1]-rep(betaRand2,nrow(x)))/sqrt(x[index,2]^2+rep(sRand2^2,nrow(x)))
		}else{
			delta=(x[1]-betaRand2)/sqrt(x[2]^2+sRand2^2)
			}
	bool=!between(delta,deltaRange[1],deltaRange[2]) #Squash outliers
	delta[bool]=deltaRange[1*(sign(delta[bool])>0)+1]
	return(delta)
	}

make_poly<-function(x,y){ #Function to make polynomials
	y^x
		}

make_lambdaDeriv2m<-function(x,y){ #y should be a matrix, even if it just has one row
  numCols=ncol(y)
  t(matrix(rep(x,numCols),ncol=numCols,byrow=FALSE)*y)
}

#' BFDR.priorsplitteR() function by Daniel Crouch
#' Estimates BFDR from a set of effect size estimates and standard errors. Effect size estimates are assumed to have normal sampling error.
#' Returns local and tail area BFDRs, BDRs and FDRs. 
#' See Crouch et al. 2022, https://doi.org/10.1101/2021.02.05.429962, for details.
#' @param beta Vector of effect size estimates with normal sampling error. Must be supplied 
#' @param s Vector of standard error estimates, assumed to be close to their true values. Must be supplied
#' @param seed Seed for random number generation. Must be supplied
#' @param excludeFromFit Vector of integers specifying the variables to exclude from distribution modelling
#' @param simulationTruth Vector of true effect sizes if known from a simulation, for plotting
#' @param outFileStem File name stem for output files. Must be supplied
#' @param breaks Number of breaks in histogram for Poisson regression distribution modelling, default 100
#' @param degree_altNull Degree of polynomial model for ANE model, default 6
#' @param degree_doublePolynomial Degree of polynomial-like model for NUPE and NUPE-diff models, default 6
#' @param nSim Number of simulated z-scores for NUPE and NUPE-diff models, default 5 million
#' @param diffSampleSize Number of 'k' variables to sample for estimating pairwise posterior probabilities for SNP k exceeding focal SNP j in effect size, default 1000
#' @param diffSampleSizeHighAcc Number of additional 'k' variables to sample for estimating pairwise posterior probabilities for SNP k exceeding focal SNP j in effect size, when estimated local bfdr_j<highAccThresh, default 9000 giving 9000+1000=10000 total samples
#' @param highAccThresh Threshold for re-sampling variables to improve estimates of bfdr_j based on whether it is low enough to be potentially interesting, default 0.05
#' @param blockSize Number of 'j' variables to vectorise bfdr_j for, default 1000. When more memory is available, blockSize can be increased with potential improvements in speed
#' @param EMtol Vector of length 2 containing tolerances for likelihood convergence in EM algorithm (ANE and NUPE/NUPE-diff respectively), default c(1e-6,1e-6)
#' @param nStop Number of consecutive EM iterations that must satisfy EM likelihood convergence tolerance before convergence is accepted, default 5
#' @param maxEMit Vector of length 2 containing maximum number of EM iterations before convergence is accepted, for ANE and NUPE/NUPE-diff respectively, default c(100000,100000). Estimates obtained when the maximum number of iterations are reached should be treated with caution
#' @param ml.reltol Maximum likelihood tolerance, passed to constrOptim, default 1e-300
#' @param outlierSD Variables with absolute z-scores greater than this threshold are removed from distribution modelling. For producing BFDR estimates, they are squashed to have absolute values no greater than the maximum value used for modelling. Default 15
#' @param diffOutlierSD Pairwise effect estimate differences with absolute z-scores greater than this threshold are removed from distribution modelling. For producing BFDR estimates, they are squashed to have absolute values no greater than the maximum value used for modelling. Default NULL, which sets the threshold to outlierSD*1.5
#' @param truncatedNull Range within which observed estimates are known a priori to exist, a vector of length 2. For example, if it is known that estimates can only be positive, specify c(0,Inf). Default NULL for range truncation.
#' @param nullTol Values with null densities equal to or below this value are omitted from fitting and assumed to have either positive or negative true effects (for positive and negative estimates respectively) with posterior probability 1, and similarly for effect differences. The estimation prodecure fits a model as a ratio against the null distribution, so when null density is recorded as zero no meaningful fit can be produced, and when it is close to zero may be unreliable.
#' @param threads Number of threads for parallelising Bayesian computations for pairwise variable effect size comparisons, default 1. We recommended higher settings when multiple threads are available
#' @param omitDist Distributions are not used for subsequent fitting if their prior probability falls below this value, to save unnecessary computation. Default 1e-5.

#' @export

BFDR.priorsplitteR<-function(beta,s,seed,excludeFromFit=NULL,simulationTruth=NULL,outFileStem,breaks=100,degree_altNull=6,degree_doublePolynomial=6,nSim=5000000,diffSampleSize=1000,diffSampleSizeHighAcc=9000,highAccThresh=0.05,blockSize=500,EMtol=c(1e-6,1e-6),nStop=5,maxEMit=c(100000,100000),ml.reltol=1e-300,outlierSD=30,diffOutlierSD=NULL,truncatedNull=NULL,nullTol=1e-300,threads=1){
library(distr);library(parallel);library(dplyr)
set.seed(seed)
outputLog=paste("outFileStem=",outFileStem,"::breaks=",breaks,"::nSim=",nSim,"::diffSampleSize=",diffSampleSize,"::degree_altNull=",degree_altNull,"::degree_doublePolynomial=",degree_doublePolynomial,"::blockSize=",blockSize,"::EMtol=",EMtol,"::nStop=",nStop,"::ml.reltol=",ml.reltol,"::outlierSD=",outlierSD,"::threads=",threads)
if((nchar(outFileStem)>0)&(substr(outFileStem,nchar(outFileStem),nchar(outFileStem))!="_"))outFileStem=paste0(outFileStem,"_")
z=beta/s ; nSNP=length(beta)
fullData=data.frame(beta=beta,s=s,z=z) #Store full data before removing outliers
#For outliers, can't rely on the polynomials outside the window, so set the outlier z-scores to the outlier threshold (w.r.t null distribution), and recalculate their betas from their s and new z
increment=0
repeat{  #Loop here to keep removing zs until there are no 'midPoints' with 0 null density. A precaution that shouldn't normally be necessary.
	z=fullData[,'z'];beta=fullData[,'beta'];s=fullData[,'s'];nSNP=nrow(fullData)
	indexNonOutliers=(1:nSNP)[(abs(z)<=outlierSD)&(dnorm(z+increment)>0)&(dnorm(z-increment)>0)]
	indexNonOutliers_nonExclude=indexNonOutliers[!indexNonOutliers%in%excludeFromFit]
	z[-indexNonOutliers]=range(z[indexNonOutliers_nonExclude])[1*(sign(z[-indexNonOutliers])>0)+1]
	beta[-indexNonOutliers]=z[-indexNonOutliers]*s[-indexNonOutliers]
	#Remove outliers and excluded SNPs only after storing full data
	fullData=data.frame(fullData,beta_squashed=beta,z_squashed=z) #Store data before removing outliers
	outputLog=c(outputLog,paste0("There are ",nSNP-length(indexNonOutliers)," out of ",nSNP," variables exceeding the outlier Z-score threshold of ",outlierSD," or near-zero null densities. Removing these and proceeding with ",length(indexNonOutliers)," variables."))
	z=z[indexNonOutliers_nonExclude];beta=beta[indexNonOutliers_nonExclude];s=s[indexNonOutliers_nonExclude];nSNP=length(indexNonOutliers_nonExclude);simulationTruth=simulationTruth[indexNonOutliers_nonExclude]
	by=round(diff(range(z))/breaks,digits=2)
	if(sum(z<0)>0){negBins=seq(-by/2,min(z)-by,by=-by)}else{negBins=NULL}
	if(sum(z>0)>0){posBins=seq(by/2,max(z)+by,by=by)}else{posBins=NULL}
	breakPoints=sort(c(negBins,posBins));breakPoints_store=breakPoints
	histogram<-hist(z,breaks=breakPoints,plot=F);histogram_store=histogram #Call histogram but don't plot
	midPoints<-histogram$mids
	if(sum(dnorm(midPoints)<=nullTol)==0)break
	increment=increment+1e-5
		}
counts<-histogram$counts
degree<-degree_altNull
lambdaDeriv1=cbind(1,matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree))
scaleMeans=rep(0,ncol(lambdaDeriv1));scaleSDs=apply(X=lambdaDeriv1,FUN=sd,MARGIN=2)
scaleList=list(ANE=scaleSDs)
lambdaDeriv1[,-1]=scale(lambdaDeriv1[,-1],center=FALSE,scale=scaleSDs[-1])
colnames(lambdaDeriv1)=NULL
width=midPoints[2]-midPoints[1]
null=pnorm(breakPoints)[-1]-pnorm(breakPoints)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
bool=(midPoints>0)
null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],lower.tail=FALSE)-pnorm(breakPoints[-1][bool],lower.tail=FALSE) 
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
m<-rep(0.5,length(counts))   #Distribution membership variables
i=1;breakCounter=0;par=c(0,log(dnorm(0))-1e-5,rep(0,degree))
repeat{
	print(paste0("Alt-Null EM iteration=",i))
	ui=matrix(0,nrow=1,ncol=degree+2);ci=c(-log(dnorm(0))+1e-10) #Intercept Less than log null minus small constant
	ui[1,2]=-1
repeat{ 
      #Very occassionally constrOptim will fail at certain initial values of par. This loop catches errors and restarts from a new initialisation point after adding some random noise
mixFit=try(constrOptim(theta=par,f=likelihood.poissonRegMixture_altNull,grad=grad.poissonRegMixture_altNull,method="BFGS",control=list(maxit=100000000,reltol=ml.reltol),ui=ui,ci=ci),silent=TRUE)
      if(!is.character(mixFit[1]))break
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
#Check fit - add more x 'midpoints' points to give a continuous appearance
x=seq(min(midPoints),max(midPoints),by=diff(range(z))/100000)
lambdaDeriv1=matrix(sapply(X=1:degree,y=x,FUN=make_poly),ncol=degree)
lambdaDeriv1=scale(lambdaDeriv1,center=scaleMeans[-1],scale=scaleSDs[-1])
lambdaDeriv1=cbind(1,lambdaDeriv1)
null=dnorm(x,mean=0,1)*width
if(!is.null(truncatedNull))null=null/diff(pnorm(truncatedNull))
predict.poissonRegMixture_altNull<-function(x){
	s=x[1];b=x[-1]
 	p=0.5*exp(s)*(s<0)+(1-0.5*exp(-s))*(s>=0)
	lambda=exp(lambdaDeriv1%*%b)
	lambdaMix=lambda*(1-p)+null*p
	return(nSNP*cbind(lambdaMix,lambda,null,lambda*(1-p),null*p))
		}
fit<-predict.poissonRegMixture_altNull(par);fit[is.na(fit)]=0
print("Printing graph of fitted alt/null mixture distribution")
#Put error bars on this fit
png(paste0(outFileStem,"poissonMixtureFit_altNull.png"),width=3000,height=3000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
plot(histogram,main=NULL,ylim=c(0,max(fit[,1:3]))*1.05,ylab="Count")
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,2],col=3,type="l",lwd=2)
points(x,fit[,3],col=4,type="l",lwd=2)
dev.off()
png(paste0(outFileStem,"poissonMixtureFit_altNull_mixingProp.png"),width=3000,height=3000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
plot(histogram,main=NULL,ylim=c(0,max(fit[,c(1,4,5)]))*1.05,ylab="Count")
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,4],col=3,type="l",lwd=2)
points(x,fit[,5],col=4,type="l",lwd=2)
dev.off()
#Now simulate from the alternative distribution
par_ANE=par
s=par[1];b=par[-1]
pNull=0.5*exp(s)*(s<0)+(1-0.5*exp(-s))*(s>=0)
altNull_Density<-function(x){
	lambdaDeriv1=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
	lambdaDeriv1=scale(lambdaDeriv1,center=scaleMeans[-1],scale=scaleSDs[-1])
	lambdaDeriv1=cbind(1,matrix(lambdaDeriv1,nrow=nrow(lambdaDeriv1)))
	lambda=exp(lambdaDeriv1%*%b)/width
	lambda[is.na(lambda)]=0 #Infinities can occasionally cause problems
	return(as.vector(lambda))
	}
fullData=data.frame(fullData,fdr=dnorm(fullData[,'z_squashed'])*pNull/(altNull_Density(fullData[,'z_squashed'])*(1-pNull)+dnorm(fullData[,'z_squashed'])*pNull)) #Store for later
fdr=fullData[indexNonOutliers_nonExclude,'fdr']
print(paste0("Simulating ",nSim," variants from estimated distribution function. Printing histogram of simulation"))
#Following few lines from 'https://stackoverflow.com/questions/23570952/simulate-from-an-arbitrary-continuous-probability-distribution'
dist <-AbscontDistribution(d=altNull_Density,low1=min(midPoints),up1=max(midPoints),withStand=TRUE) #signature for a dist with pdf ~ twoGroupDensity
rdist <- r(dist) # function to create random variates 
zSim=rdist(nSim)
png(paste0(outFileStem,"poissonMixtureFit_altNull_altSim.png"),width=3000,height=6000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),mfrow=c(2,1),cex.axis=0.7,cex.lab=0.7)
x<-seq(min(midPoints),max(midPoints),by=0.01)
y<-unlist(sapply(X=x,FUN=altNull_Density))
denom=mean(y)*diff(range(x))
y=y/denom
yNull=dnorm(x,mean=0,sd=1)
hist(zSim,freq=FALSE,breaks=breakPoints,col=0,ylim=c(0,max(c(y,yNull))),main=NULL,xlab="",xlim=range(zSim),ylab="Count")
lines(x,y,type='l',col="green")
lines(x,yNull,type='l',col="blue")
hist(zSim,freq=FALSE,breaks=100,col=0,ylim=c(0,max(c(y,yNull))),main=NULL,xlab="",xlim=range(zSim),ylab="Count")
lines(x,y,type='l',col="green")
lines(x,yNull,type='l',col="blue")
dev.off()
#Simulate differences by drawing nSim further variables.
#Perform NNPE estimation for finding the false sign probabilities. This is the same as the NNPED model above but using SNP effects rather than simulated effect differences
#Define histogram
degree<-degree_doublePolynomial
sDelta=1
histogram<-hist(zSim,breaks=breakPoints,plot=FALSE);histogram_storeDelta=histogram #Call histogram but don't plot
counts<-histogram$counts
lambdaDeriv2=matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree)
scaleMeans=rep(0,ncol(lambdaDeriv2));scaleSDs=apply(lambdaDeriv2,FUN=sd,MAR=2)
scaleList=append(scaleList,list(NUPE=scaleSDs))
lambdaDeriv2=scale(lambdaDeriv2,center=FALSE,scale=scaleSDs) #Don't center first term
colnames(lambdaDeriv2)=NULL
lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)))
lambdaDeriv2m=lambdaDeriv2m/ij
ij=matrix(ij[1,],nrow=degree+1,byrow=FALSE) #For gradient function
width=midPoints[2]-midPoints[1]
nullFactor=dnorm(midPoints,mean=0,sd=sDelta)*width
null=pnorm(breakPoints,sd=sDelta)[-1]-pnorm(breakPoints,sd=sDelta)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
bool=(midPoints>0)
null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],sd=sDelta,lower.tail=FALSE)-pnorm(breakPoints[-1][bool],sd=sDelta,lower.tail=FALSE) 
null=nullFactor
EM<-function(degree,counts,lambdaDeriv2,lambdaDeriv2m,ij,nullFactor,null,width,printMessage=NULL){
  #Define functions for ML optimisation of NNPED and NNPE models
  likelihood.poissonRegMixture_twoPolynomial<-function(x){
    s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    p3=1-p1-p2;p3[p3<0]=0
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
    lambdaU=exp(polyU)[,1]*nullFactor
    lambdaL=exp(polyL)[,1]*nullFactor
    l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*n
    l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*n
    l3=counts*(1-m1-m2)*log(p3)-p3*null*n
    l1[(l1==(-Inf))|is.na(l1)]=-Inf
    l2[(l2==(-Inf))|is.na(l2)]=-Inf
    l3[(l3==(-Inf))|is.na(l3)]=-Inf
    l=sum(l1+l2+l3)
    return(-l)
  }
  grad.poissonRegMixture_twoPolynomial<-function(x){
    s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    p3=1-p1-p2;p3[p3<0]=0
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
    lambdaU=exp(polyU)[,1]*nullFactor
    lambdaL=exp(polyL)[,1]*nullFactor
    coeffDerivU=2*cbind(1,lambdaDeriv2)*cbind(1,lambdaDeriv2)%*%(matrix(rep(c(c1a,b1),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    coeffDerivL=2*cbind(1,lambdaDeriv2)*cbind(1,lambdaDeriv2)%*%(matrix(rep(c(c2a,b2),degree+1),ncol=degree+1,byrow=FALSE)/ij)
    deriv=c(sum(counts*m1/p1-counts*(1-m1-m2)/p3-lambdaL*n+null*n),sum(counts*m2/p2-counts*(1-m1-m2)/p3-lambdaU*n+null*n),((counts*m2-p2*lambdaU*n)*lambdaDeriv2[,1])%*%coeffDerivU[,-1],-((counts*m1-p1*lambdaL*n)*lambdaDeriv2[,1])%*%coeffDerivL[,-1],sum(counts*m2-p2*lambdaU*n),sum(counts*m1-p1*lambdaL*n),((counts*m2-p2*lambdaU*n)*lambdaDeriv2[,1])%*%coeffDerivU[,1],-((counts*m1-p1*lambdaL*n)*lambdaDeriv2[,1])%*%coeffDerivL[,1])
    deriv[1:2]=deriv[1:2]*c(0.5*(exp(s1)*(s1<0)+exp(-s1)*(s1>=0)),0.5*(exp(s2)*(s2<0)+exp(-s2)*(s2>=0))) 
    return(-deriv)
  }
  lNull=log(null) #Helps avoid numerical issues in function below
  n<-sum(counts)
  m1<-rep(0.33,length(counts))   #Lower tail responsibilities 
  m2<-rep(0.33,length(counts))   #Upper tail responsibilities 
  i=1;breakCounter=0;par=c(-0.1,-0.1,rep(0,degree),rep(0,degree),-1e-5,-1e-5,1e-5,1e-5)
  repeat{
    print(paste0("Two-polynomial EM iteration (",printMessage,") =",i))
    ui=matrix(0,nrow=5,ncol=degree*2+6);ci=c(0,0,0,0,0)
    ui[1,1:2]=c(-1,-1)  #s1 plus s2 <=0 implies p1+p2<1.
    diag(ui[2:3,degree*2+3:4])=c(-1,-1) #Intercept at 0 constraints - ratio with null should be < 1
    diag(ui[4:5,2*degree+5:6])=c(1,1) #Second intercept should be >0 to stop potential symmetry/identifiability issues
    par_tweak=par #Modify par to avoid boundary constraint issues
    par_tweak[2*degree+3:4][par[2*degree+3:4]>(-1e-2)]=-1e-2#Initialise these two parameters comfortably within the boundary region
    par_tweak[2*degree+5:6][par[2*degree+5:6]<(1e-2)]=1e-2#Initialise these two parameters comfortably within the boundary region
    repeat{ 
      #Very occassionally constrOptim will fail at certain initial values of par. This loop catches errors and restarts from a new initialisation point after adding some random noise
	mixFit=try(constrOptim(theta=par_tweak,f=likelihood.poissonRegMixture_twoPolynomial,grad=grad.poissonRegMixture_twoPolynomial,method="BFGS",control=list(maxit=100000000,reltol=ml.reltol),ui=ui,ci=ci),silent=TRUE)
      if(!is.character(mixFit[1]))break
      par_tweak=par_tweak+rnorm(length(par_tweak),sd=0.01)
      par_tweak[2*degree+3:4][par_tweak[2*degree+3:4]>(-1e-5)]=-1e-05 
      par_tweak[2*degree+5:6][par_tweak[2*degree+5:6]<1e-5]=1e-05 
      print("Error in likelihood maximisation. Adding random noise to starting parameters and retrying")
    }
    x=par  
    s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] 
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    p3=1-p1-p2;p3[p3<0]=0
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
    lambdaU=exp(polyU)[,1]*nullFactor
    lambdaL=exp(polyL)[,1]*nullFactor
    l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*n
    l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*n
    l3=counts*(1-m1-m2)*log(p3)-p3*null*n
    l1[(l1==(-Inf))|is.na(l1)]=-1e300
    l2[(l2==(-Inf))|is.na(l2)]=-1e300
    l3[(l3==(-Inf))|is.na(l3)]=-1e300
    Qjj=sum(l1+l2+l3)
    x=mixFit$par  #Find core expected likelihood with given the parameters
    s1=x[1];s2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] 
    p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
    p3=1-p1-p2;p3[p3<0]=0
    polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
    polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
    lambdaU=exp(polyU)[,1]*nullFactor
    lambdaL=exp(polyL)[,1]*nullFactor
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
    par=mixFit$par
    sumLambdaU=sum(lambdaU);sumLambdaL=sum(lambdaL)
    if(((QjjPlus1-Qjj)<0)){message="Warning: M step has failed to improve on previous optimum"}else{message=""}
    print(paste0("Core expected likelihood update = ",QjjPlus1,"-",Qjj,"=",QjjPlus1-Qjj," ",message))
    if(((QjjPlus1-Qjj)<EMtol[2])&((QjjPlus1-Qjj)>=0)){breakCounter=breakCounter+1}else{breakCounter=0}
    if(breakCounter==nStop)break
    if(i>=maxEMit[2]){
      message=paste0("Warning: Maximum iterations (",maxEMit[2],") reached")
      outputLog=c(outputLog,message)
      print(message)
      break
    }
    i=i+1
  } #End of repeat loops
  return(list(par=par,sumLambdaU=sum(lambdaU),sumLambdaL=sum(lambdaL)))
}#End of function 'EM'
EMout=EM(degree=degree_doublePolynomial,counts=counts,lambdaDeriv2=lambdaDeriv2,lambdaDeriv2m=lambdaDeriv2m,ij=ij,nullFactor=nullFactor,null=null,width=width,printMessage="sign mixture")
par=par_NNPE=EMout$par;sumLambdaU=EMout$sumLambdaU;sumLambdaL=EMout$sumLambdaL
sPos=par[2];sNeg=par[1]
piPos=0.5*exp(sPos)*(sPos<0)+(1-0.5*exp(-sPos))*(sPos>=0)
piNeg=0.5*exp(sNeg)*(sNeg<0)+(1-0.5*exp(-sNeg))*(sNeg>=0)
#Check fit
print("Printing graph of fitted sign-probability mixture distribution")
predict.poissonRegMixture_twoPolynomial<-function(x,figStem,n,width,midPoints,scaleMeans,scaleSDs){
  par=x
  x=seq(min(midPoints),max(midPoints),by=diff(range(midPoints))/100000)
  lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),ncol=degree)
  lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
  lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
  lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
  ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
  nullFactor=dnorm(x,mean=0,sd=sDelta)*width
  null=nullFactor #For visualisation purposes this is fine - multiply by width inside the following function
  s1=par[1];s2=par[2];b1=par[3:(degree+2)];b2=par[(degree+3):(2*degree+2)];c1=par[2*degree+3];c2=par[2*degree+4];c1a=par[2*degree+5];c2a=par[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
  p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
  p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
  p3=1-p1-p2;p3[p3<0]=0
  polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
  polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
  lambdaU=exp(polyU)[,1]*nullFactor
  lambdaL=exp(polyL)[,1]*nullFactor
  lambdaMix=lambdaL*p1+lambdaU*p2+null*p3
  fit2=n*cbind(lambdaMix,lambdaL,lambdaU,null,lambdaL*p1,lambdaU*p2,null*p3)
  fit=fit2
  fit[is.na(fit)]=0
  png(paste0(outFileStem,figStem,".png"),width=3000,height=3000,res=600)
  par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
  plot(histogram,main=NULL,ylim=c(0,max(fit[,-(5:7)]))*1.05,ylab="Count")
  points(x,fit[,1],col=2,type="l",lwd=2)
  points(x,fit[,2],col=3,type="l",lwd=2)
  points(x,fit[,3],col=5,type="l",lwd=2)
  points(x,fit[,4],col=4,type="l",lwd=2)
  dev.off()
  png(paste0(outFileStem,figStem,"_mixingPropScaled.png"),width=3000,height=3000,res=600)
  par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
  plot(histogram,main=NULL,ylim=c(0,max(fit[,-(2:4)]))*1.05,ylab="Count")
  points(x,fit[,1],col=2,type="l",lwd=2)
  points(x,fit[,5],col=3,type="l",lwd=2)
  points(x,fit[,6],col=5,type="l",lwd=2)
  points(x,fit[,7],col=4,type="l",lwd=2)
  dev.off()
  return(fit2)
}
fit<-predict.poissonRegMixture_twoPolynomial(x=par,figStem="poissonMixtureFit_falseSign",n=nSim,width=width,midPoints=midPoints,scaleMeans=scaleMeans,scaleSDs=scaleSDs);fit[is.na(fit)]=0
fit_nupe<-fit
sign_probability<-function(x){
  s1=par_NNPE[1];s2=par_NNPE[2];b1=par_NNPE[3:(degree+2)];b2=par_NNPE[(degree+3):(2*degree+2)];c1=par_NNPE[2*degree+3];c2=par_NNPE[2*degree+4];c1a=par_NNPE[2*degree+5];c2a=par_NNPE[2*degree+6] 
  p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
  p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
  lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
  lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
  lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
  lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
  ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
  nullFactor=dnorm(x,mean=0,sd=sDelta)
  polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
  polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
  lambdaU=exp(polyU)[,1]*nullFactor
  lambdaL=exp(polyL)[,1]*nullFactor
  posProb=lambdaU*p2/(lambdaL*p1+lambdaU*p2) #Ignore the null probability as, when we need to calculate signs, it is conditioned on a SNP being alt
  posProb[is.na(posProb)]=0 #Infinities can occasionally cause problems
  return(as.vector(posProb))
}
sign_probability2<-function(x){
  s1=par_NNPE[1];s2=par_NNPE[2];b1=par_NNPE[3:(degree+2)];b2=par_NNPE[(degree+3):(2*degree+2)];c1=par_NNPE[2*degree+3];c2=par_NNPE[2*degree+4];c1a=par_NNPE[2*degree+5];c2a=par_NNPE[2*degree+6] 
  p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
  p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
  lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
  lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
  lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
  lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
  ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
  nullFactor=dnorm(x,mean=0,sd=sDelta)
  polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
  polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
  lambdaU=exp(polyU)[,1]*nullFactor
  lambdaL=exp(polyL)[,1]*nullFactor
  posProb=cbind(lambdaU*p2/(lambdaL*p1+lambdaU*p2+nullFactor*(1-p1-p2)),lambdaL*p1/(lambdaL*p1+lambdaU*p2+nullFactor*(1-p1-p2)))
  posProb[is.na(posProb)]=0 #Infinities can occasionally cause problems
  return(posProb)
}
posProb2=sign_probability2(fullData[,'z_squashed'])
fullData=data.frame(fullData,positiveEffectProb=posProb2[,1],negativeEffectProb=posProb2[,2]) #Store for later
posProb2=posProb2[indexNonOutliers_nonExclude,]
upper_Density<-function(x){
  s1=par_NNPE[1];s2=par_NNPE[2];b1=par_NNPE[3:(degree+2)];b2=par_NNPE[(degree+3):(2*degree+2)];c1=par_NNPE[2*degree+3];c2=par_NNPE[2*degree+4];c1a=par_NNPE[2*degree+5];c2a=par_NNPE[2*degree+6] 
  p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
  p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
  lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
  lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
  lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
  lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
  ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
  nullFactor=dnorm(x,mean=0,sd=sDelta)
  polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
  lambdaU=exp(polyU)[,1]*nullFactor
  lambdaU=lambdaU/(sumLambdaU)
  return(as.vector(lambdaU))
}
lower_Density<-function(x){
  s1=par_NNPE[1];s2=par_NNPE[2];b1=par_NNPE[3:(degree+2)];b2=par_NNPE[(degree+3):(2*degree+2)];c1=par_NNPE[2*degree+3];c2=par_NNPE[2*degree+4];c1a=par_NNPE[2*degree+5];c2a=par_NNPE[2*degree+6] 
  p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
  p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
  lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
  lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
  lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
  lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
  ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
  nullFactor=dnorm(x,mean=0,sd=sDelta)
  polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
  lambdaL=exp(polyL)[,1]*nullFactor
  lambdaL=lambdaL/(sumLambdaL)
  return(as.vector(lambdaL))
}
print(paste0("Simulating ",nSim," variants from positive effect SNPs. Printing histogram of simulation"))
#Following few lines from 'https://stackoverflow.com/questions/23570952/simulate-from-an-arbitrary-continuous-probability-distribution'
dist <-AbscontDistribution(d=upper_Density,low1=min(midPoints),up1=max(midPoints),withStand=TRUE) #signature for a dist with pdf ~ twoGroupDensity
rdist <- r(dist) # function to create random variates 
zSimUpper=rdist(nSim)
png(paste0(outFileStem,"poissonMixtureFit_sign_upperSim.png"),width=3000,height=3000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
x<-seq(min(midPoints),max(midPoints),by=0.01)
y<-unlist(sapply(X=x,FUN=upper_Density))
denom=mean(y)*diff(range(x))
y=y/denom
yNull=dnorm(x,mean=0,sd=1)
hist(zSimUpper,freq=FALSE,breaks=breakPoints,col=0,ylim=c(0,max(c(y,yNull))),main=NULL,xlab="",xlim=range(z),ylab="Count")
lines(x,y,type='l',col="green")
lines(x,yNull,type='l',col="blue")
dev.off()
dist <-AbscontDistribution(d=lower_Density,low1=min(midPoints),up1=max(midPoints),withStand=TRUE) #signature for a dist with pdf ~ twoGroupDensity
rdist <- r(dist) # function to create random variates 
zSimLower=rdist(nSim)
png(paste0(outFileStem,"poissonMixtureFit_sign_lowerSim.png"),width=3000,height=3000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
x<-seq(min(midPoints),max(midPoints),by=0.01)
y<-unlist(sapply(X=x,FUN=lower_Density))
denom=mean(y)*diff(range(x))
y=y/denom
yNull=dnorm(x,mean=0,sd=1)
hist(zSimLower,freq=FALSE,breaks=breakPoints_store,col=0,ylim=c(0,max(c(y,yNull))),main=NULL,xlab="",xlim=range(z),ylab="Count")
lines(x,y,type='l',col="green")
lines(x,yNull,type='l',col="blue")
dev.off()
#Perform upper NNPED estimation
index1=sample(rep(1:nSNP,rmultinom(1,size=nSim,prob=(posProb2[,1]/(piPos*nSNP))*(1-fdr)/((1-pNull)*nSNP))[,1]),replace=FALSE,size=nSim) #Sample of observed SNPs as if they were from alt #distribution, as before when choosing alt SNPs for plotting
index2=sample(rep(1:nSNP,rmultinom(1,size=nSim,prob=(1-fdr)/((1-pNull)*nSNP))[,1]),replace=FALSE,size=nSim)
s1=fullData[index1,'s'];s2=fullData[index2,'s']
delta=(zSimUpper*s1-zSim*s2)/sqrt(s1^2+s2^2)
if(!is.null(simulationTruth)){ 
index=sample(rep(1:nSNP,rmultinom(1,size=nSim,prob=(posProb2[,1]/(piPos*nSNP))*(1-fdr)/((1-pNull)*nSNP))[,1]),replace=FALSE,size=nSim) #Sample of observed SNPs as if they were from alt #distribution, as before when choosing alt SNPs for plotting
rand=sample(rep(1:nSNP,rmultinom(1,size=nSim,prob=(1-fdr)/((1-pNull)*nSNP))[,1]),replace=FALSE,size=nSim)
rand=sample(rand,size=length(rand),replace=FALSE) #Shuffle order
deltaTrue=(beta[index]-beta[rand])/sqrt(s[index]^2+s[rand]^2)
trueDiff=(simulationTruth[index]-simulationTruth[rand])/sqrt(s[index]^2+s[rand]^2)
png(paste0(outFileStem,"poissonMixtureFit_altNull_diffs_simVersusEstAlt.png"),width=3000,height=3000,res=600)
par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
deltaHist=hist(c(delta,deltaTrue),breaks=breaks,plot=FALSE)
deltaHist1=hist(delta,breaks=deltaHist$breaks,plot=FALSE);deltaHist2=hist(deltaTrue,breaks=deltaHist$breaks,plot=FALSE) #For count range
hist(delta,breaks=deltaHist$breaks,col=rgb(0,0,1,alpha=0.5),main=NULL,xlab="Effect size differences",xlim=range(c(delta,deltaTrue)),ylim=1.2*c(0,max(c(deltaHist1$counts,deltaHist2$counts))),ylab="Count")
hist(deltaTrue,breaks=deltaHist$breaks,add=T,col=rgb(1,0,0,alpha=0.5))	
legend(x=min(c(delta,deltaTrue)),y=1.2*max(c(deltaHist1$counts,deltaHist2$counts)),col=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5)),legend=c("Drawn directly from estimated alt distribution","Empirically estimated alt distribution based on fdrs"),pch=15)
dev.off()
deltaTrue_altEst=deltaTrue;trueDiff_altEst=trueDiff #Store for later #Store for later
		}
#Fit the 2 polynomial model, again with poisson regression. First remove any extreme outliers that exist in the deltas - should't be many
if(is.null(diffOutlierSD))diffOutlierSD=outlierSD*1.5
deltaMax=diffOutlierSD
bool=(((abs(delta)>=deltaMax)+dnorm(delta)==0)>0)
delta=delta[!bool]
deltaRange_U=range(delta) #Will also use this again later in Bayes theorem computations
if(sum(bool)>0){message=paste0("Warning: there are ",sum(bool)," extreme outlier deltas (out of total ",nSim,") with absolute values greater than than ",deltaMax," or near-zero null densities. Removing these. Consider lowering the outlierSD parameter if the number of these seems high.")}else{message=NULL}
outputLog=c(outputLog,message);if(!is.null(message))print(message)
sDelta=1 #sds of null differences
by=round(diff(range(delta))/breaks,digits=2)
if(sum(delta<0)>0){negBins=seq(-by/2,min(delta)-by,by=-by)}else{negBins=NULL}
if(sum(delta>0)>0){posBins=seq(by/2,max(delta)+by,by=by)}else{posBins=NULL}
breakPoints=sort(c(negBins,posBins))
histogram<-hist(delta,breaks=breakPoints,plot=F)#Call histogram but don't plot
midPoints<-histogram$mids
delta=delta[between(delta,min(breakPoints[-length(breakPoints)][dnorm(midPoints)>nullTol]),max(breakPoints[-1][dnorm(midPoints)>nullTol]))]
breakPoints=sort(unique(c(breakPoints[-1][dnorm(midPoints)>nullTol],breakPoints[-length(breakPoints)][dnorm(midPoints)>nullTol])));breakPoints_storeDelta=breakPoints
histogram<-hist(delta,breaks=breakPoints,plot=F);histogram_storeDelta=histogram #Call histogram but don't plot
midPoints<-histogram$mids
counts<-histogram$counts
degree<-degree_doublePolynomial
lambdaDeriv2=matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree)
scaleMeans=rep(0,ncol(lambdaDeriv2));scaleSDs=apply(lambdaDeriv2,FUN=sd,MAR=2)
scaleList=append(scaleList,list(NUPEd_upper=scaleSDs))
lambdaDeriv2=scale(lambdaDeriv2,center=FALSE,scale=scaleSDs) #Don't center 
colnames(lambdaDeriv2)=NULL
make_lambdaDeriv2m<-function(x,y){
		t(matrix(rep(x,ncol(y)),ncol=ncol(y),byrow=FALSE)*y)
			}
lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)))
lambdaDeriv2m=lambdaDeriv2m/ij
ij=matrix(ij[1,],nrow=degree+1,byrow=FALSE) #For gradient function
#Use EM algorithm to fit two-polynomial constrained mixture density
width=midPoints[2]-midPoints[1]
nullFactor=dnorm(midPoints,mean=0,sd=sDelta)*width
null=pnorm(breakPoints,sd=sDelta)[-1]-pnorm(breakPoints,sd=sDelta)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
bool=(midPoints>0)
null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],sd=sDelta,lower.tail=FALSE)-pnorm(breakPoints[-1][bool],sd=sDelta,lower.tail=FALSE) 
null=nullFactor
EMout=EM(degree=degree_doublePolynomial,counts=counts,lambdaDeriv2=lambdaDeriv2,lambdaDeriv2m=lambdaDeriv2m,ij=ij,nullFactor=nullFactor,null=null,width=width,printMessage="positive effects, observed sign")
par=par_NNPED_upper=EMout$par;sumLambdaU=EMout$sumLambdaU;sumLambdaL=EMout$sumLambdaL
#Check fit
null=nullFactor #For visualisation purposes this is fine - multiply by width inside the following function
print("Printing graph of fitted two-polynomial mixture distribution")
fit<-predict.poissonRegMixture_twoPolynomial(par,figStem="poissonMixtureFit_twoPolynomial_upper",n=nSim,width=width,midPoints=midPoints,scaleMeans=scaleMeans,scaleSDs=scaleSDs);fit[is.na(fit)]=0
x=seq(min(midPoints),max(midPoints),by=diff(range(midPoints))/100000)
if(!is.null(simulationTruth)){ 
  deltaTrue=deltaTrue_altEst;trueDiff=trueDiff_altEst
  rand=1:nSim#rep(1:nSim,rmultinom(1,size=nSim,prob=(1-twoPolynomial_fdr(x=deltaTrue,y=par))/((par[1]+par[2])*nSim))[,1])
  deltaTrue=deltaTrue[rand];trueDiff=trueDiff[rand]
  keep=((deltaTrue<=max(histogram$breaks))&((deltaTrue>=min(histogram$breaks))))  #Avoids graphing problems
  deltaTrue=deltaTrue[keep];trueDiff=trueDiff[keep]
  deltaHist1=hist(delta,breaks=histogram$breaks,plot=FALSE);deltaHist2=hist(deltaTrue[trueDiff<0],breaks=histogram$breaks,plot=FALSE);deltaHist3=hist(deltaTrue[trueDiff>0],breaks=histogram$breaks,plot=FALSE) #For count range
  png(paste0(outFileStem,"poissonMixtureFit_twoPolynomial_mixingPropScaled_withEstAltDiffs.png"),width=3000,height=3000,res=600)
  par(mar=c(3.6,4.1,1.1,1.1),mgp=c(2.5,1,0),cex.axis=0.7,cex.lab=0.7)
  plot(histogram,main=NULL,ylim=1.2*c(0,max(c(deltaHist1$counts,deltaHist2$counts,deltaHist3$counts))),ylab="Count")
  hist(deltaTrue[trueDiff>0],breaks=histogram$breaks,col=rgb(1,0,0,alpha=0.5),add=T)
  hist(deltaTrue[trueDiff<0],breaks=histogram$breaks,col=rgb(0,0,1,alpha=0.5),add=T)
  hist(deltaTrue[trueDiff==0],breaks=histogram$breaks,col=rgb(0,1,0,alpha=0.25),add=T)
  points(x,fit[,1],col=2,type="l",lwd=2)
  points(x,fit[,5],col=3,type="l",lwd=2)
  points(x,fit[,6],col=5,type="l",lwd=2)
  points(x,fit[,7],col=4,type="l",lwd=2)
  legend(x=min(midPoints),y=1.2*max(c(deltaHist1$counts,deltaHist2$counts)),col=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,1,0,alpha=0.25)),legend=c("True effect < 0","True effect > 0","True effect = 0"),pch=15,cex=0.6)
  dev.off()
} #End of 'if greater than omitFit'
#Reverse the sign of the alt distribution and repeat
delta=(zSimUpper*s1+zSim*s2)/sqrt(s1^2+s2^2)
if(is.null(diffOutlierSD))diffOutlierSD=outlierSD*1.5
deltaMax=diffOutlierSD 
bool=(((abs(delta)>=deltaMax)+dnorm(delta)==0)>0)
delta=delta[!bool]
deltaRange_U_rev=range(delta) #Will also use this again later in Bayes theorem computations
if(sum(bool)>0){message=paste0("Warning: there are ",sum(bool)," extreme outlier deltas (out of total ",nSim,") with absolute values greater than than ",deltaMax," or near-zero null densities. Removing these. Consider lowering the outlierSD parameter if the number of these seems high.")}else{message=NULL}
outputLog=c(outputLog,message);if(!is.null(message))print(message)
sDelta=1 #sds of null differences
by=round(diff(range(delta))/breaks,digits=2)
if(sum(delta<0)>0){negBins=seq(-by/2,min(delta)-by,by=-by)}else{negBins=NULL}
if(sum(delta>0)>0){posBins=seq(by/2,max(delta)+by,by=by)}else{posBins=NULL}
breakPoints=sort(c(negBins,posBins))
histogram<-hist(delta,breaks=breakPoints,plot=F)#Call histogram but don't plot
midPoints<-histogram$mids
delta=delta[between(delta,min(breakPoints[-length(breakPoints)][dnorm(midPoints)>nullTol]),max(breakPoints[-1][dnorm(midPoints)>nullTol]))]
breakPoints=sort(unique(c(breakPoints[-1][dnorm(midPoints)>nullTol],breakPoints[-length(breakPoints)][dnorm(midPoints)>nullTol])));breakPoints_storeDelta=breakPoints
histogram<-hist(delta,breaks=breakPoints,plot=F);histogram_storeDelta=histogram #Call histogram but don't plot
midPoints<-histogram$mids
counts<-histogram$counts
degree<-degree_doublePolynomial
lambdaDeriv2=matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree)
scaleMeans=rep(0,ncol(lambdaDeriv2));scaleSDs=apply(lambdaDeriv2,FUN=sd,MAR=2)
scaleList=append(scaleList,list(NUPEd_upper_rev=scaleSDs))
lambdaDeriv2=scale(lambdaDeriv2,center=FALSE,scale=scaleSDs) #Don't center 
colnames(lambdaDeriv2)=NULL
make_lambdaDeriv2m<-function(x,y){
  t(matrix(rep(x,ncol(y)),ncol=ncol(y),byrow=FALSE)*y)
}
lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)))
lambdaDeriv2m=lambdaDeriv2m/ij
ij=matrix(ij[1,],nrow=degree+1,byrow=FALSE) #For gradient function
#Use EM algorithm to fit two-polynomial constrained mixture density
width=midPoints[2]-midPoints[1]
nullFactor=dnorm(midPoints,mean=0,sd=sDelta)*width
null=pnorm(breakPoints,sd=sDelta)[-1]-pnorm(breakPoints,sd=sDelta)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
bool=(midPoints>0)
null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],sd=sDelta,lower.tail=FALSE)-pnorm(breakPoints[-1][bool],sd=sDelta,lower.tail=FALSE) 
null=nullFactor
EMout=EM(degree=degree_doublePolynomial,counts=counts,lambdaDeriv2=lambdaDeriv2,lambdaDeriv2m=lambdaDeriv2m,ij=ij,nullFactor=nullFactor,null=null,width=width,printMessage="positive effects, reversed sign")
par=par_NNPED_upper_rev=EMout$par;sumLambdaU=EMout$sumLambdaU;sumLambdaL=EMout$sumLambdaL
#Check fit
print("Printing graph of fitted two-polynomial mixture distribution")
fit<-predict.poissonRegMixture_twoPolynomial(par,figStem="poissonMixtureFit_twoPolynomial_upper_rev",n=nSim,width=width,midPoints=midPoints,scaleMeans=scaleMeans,scaleSDs=scaleSDs);fit[is.na(fit)]=0
#Perform lower NNPED estimation
index1=sample(rep(1:nSNP,rmultinom(1,size=nSim,prob=(posProb2[,2]/(piNeg*nSNP))*(1-fdr)/((1-pNull)*nSNP))[,1]),replace=FALSE,size=nSim) #Sample of observed SNPs as if they were from alt #distribution, as before when choosing alt SNPs for plotting
s1=fullData[index1,'s']
delta=(zSimLower*s1-zSim*s2)/sqrt(s1^2+s2^2)
#Fit the 2 polynomial model, again with poisson regression. First remove any extreme outliers that exist in the deltas - should't be many
if(is.null(diffOutlierSD))diffOutlierSD=outlierSD*1.5
deltaMax=diffOutlierSD*sd(delta) 
bool=(((abs(delta)>=deltaMax)+dnorm(delta)==0)>0)
delta=delta[!bool]
deltaRange_L=range(delta) #Will also use this again later in Bayes theorem computations
if(sum(bool)>0){message=paste0("Warning: there are ",sum(bool)," extreme outlier deltas (out of total ",nSim,") with absolute values greater than than ",deltaMax," or near-zero null densities. Removing these. Consider lowering the outlierSD parameter if the number of these seems high.")}else{message=NULL}
outputLog=c(outputLog,message);if(!is.null(message))print(message)
sDelta=1 #sds of null differences
by=round(diff(range(delta))/breaks,digits=2)
if(sum(delta<0)>0){negBins=seq(-by/2,min(delta)-by,by=-by)}else{negBins=NULL}
if(sum(delta>0)>0){posBins=seq(by/2,max(delta)+by,by=by)}else{posBins=NULL}
breakPoints=sort(c(negBins,posBins))
histogram<-hist(delta,breaks=breakPoints,plot=F)#Call histogram but don't plot
midPoints<-histogram$mids
delta=delta[between(delta,min(breakPoints[-length(breakPoints)][dnorm(midPoints)>nullTol]),max(breakPoints[-1][dnorm(midPoints)>nullTol]))]
breakPoints=sort(unique(c(breakPoints[-1][dnorm(midPoints)>nullTol],breakPoints[-length(breakPoints)][dnorm(midPoints)>nullTol])));breakPoints_storeDelta=breakPoints
histogram<-hist(delta,breaks=breakPoints,plot=F);histogram_storeDelta=histogram #Call histogram but don't plot
midPoints<-histogram$mids
counts<-histogram$counts
degree<-degree_doublePolynomial
lambdaDeriv2=matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree)
scaleMeans=rep(0,ncol(lambdaDeriv2));scaleSDs=apply(lambdaDeriv2,FUN=sd,MAR=2)
scaleList=append(scaleList,list(NUPEd_lower=scaleSDs))
lambdaDeriv2=scale(lambdaDeriv2,center=FALSE,scale=scaleSDs) #Don't center 
colnames(lambdaDeriv2)=NULL
make_lambdaDeriv2m<-function(x,y){
  t(matrix(rep(x,ncol(y)),ncol=ncol(y),byrow=FALSE)*y)
}
lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)))
lambdaDeriv2m=lambdaDeriv2m/ij
ij=matrix(ij[1,],nrow=degree+1,byrow=FALSE) #For gradient function
#Use EM algorithm to fit two-polynomial constrained mixture density
width=midPoints[2]-midPoints[1]
nullFactor=dnorm(midPoints,mean=0,sd=sDelta)*width
null=pnorm(breakPoints,sd=sDelta)[-1]-pnorm(breakPoints,sd=sDelta)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
bool=(midPoints>0)
null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],sd=sDelta,lower.tail=FALSE)-pnorm(breakPoints[-1][bool],sd=sDelta,lower.tail=FALSE) 
null=nullFactor
EMout=EM(degree=degree_doublePolynomial,counts=counts,lambdaDeriv2=lambdaDeriv2,lambdaDeriv2m=lambdaDeriv2m,ij=ij,nullFactor=nullFactor,null=null,width=width,printMessage="negative effects, observed sign")
par=par_NNPED_lower=EMout$par;sumLambdaU=EMout$sumLambdaU;sumLambdaL=EMout$sumLambdaL
#Check fit
print("Printing graph of fitted two-polynomial mixture distribution")
fit<-predict.poissonRegMixture_twoPolynomial(par,figStem="poissonMixtureFit_twoPolynomial_lower",n=nSim,width=width,midPoints=midPoints,scaleMeans=scaleMeans,scaleSDs=scaleSDs);fit[is.na(fit)]=0
#Reverse sign and repeat
delta=(zSimLower*s1+zSim*s2)/sqrt(s1^2+s2^2)
if(is.null(diffOutlierSD))diffOutlierSD=outlierSD*1.5
deltaMax=diffOutlierSD*sd(delta) 
bool=(((abs(delta)>=deltaMax)+dnorm(delta)==0)>0)
delta=delta[!bool]
deltaRange_L_rev=range(delta) #Will also use this again later in Bayes theorem computations
if(sum(bool)>0){message=paste0("Warning: there are ",sum(bool)," extreme outlier deltas (out of total ",nSim,") with absolute values greater than than ",deltaMax," or near-zero null densities. Removing these. Consider lowering the outlierSD parameter if the number of these seems high.")}else{message=NULL}
outputLog=c(outputLog,message);if(!is.null(message))print(message)
sDelta=1 #sds of null differences
by=round(diff(range(delta))/breaks,digits=2)
if(sum(delta<0)>0){negBins=seq(-by/2,min(delta)-by,by=-by)}else{negBins=NULL}
if(sum(delta>0)>0){posBins=seq(by/2,max(delta)+by,by=by)}else{posBins=NULL}
breakPoints=sort(c(negBins,posBins))
histogram<-hist(delta,breaks=breakPoints,plot=F)#Call histogram but don't plot
midPoints<-histogram$mids
delta=delta[between(delta,min(breakPoints[-length(breakPoints)][dnorm(midPoints)>nullTol]),max(breakPoints[-1][dnorm(midPoints)>nullTol]))]
breakPoints=sort(unique(c(breakPoints[-1][dnorm(midPoints)>nullTol],breakPoints[-length(breakPoints)][dnorm(midPoints)>nullTol])));breakPoints_storeDelta=breakPoints
histogram<-hist(delta,breaks=breakPoints,plot=F);histogram_storeDelta=histogram #Call histogram but don't plot
midPoints<-histogram$mids
counts<-histogram$counts
degree<-degree_doublePolynomial
lambdaDeriv2=matrix(sapply(X=1:degree,y=midPoints,FUN=make_poly),ncol=degree)
scaleMeans=rep(0,ncol(lambdaDeriv2));scaleSDs=apply(lambdaDeriv2,FUN=sd,MAR=2)
scaleList=append(scaleList,list(NUPEd_lower_rev=scaleSDs))
lambdaDeriv2=scale(lambdaDeriv2,center=FALSE,scale=scaleSDs) #Don't center 
colnames(lambdaDeriv2)=NULL
make_lambdaDeriv2m<-function(x,y){
  t(matrix(rep(x,ncol(y)),ncol=ncol(y),byrow=FALSE)*y)
}
lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)))
lambdaDeriv2m=lambdaDeriv2m/ij
ij=matrix(ij[1,],nrow=degree+1,byrow=FALSE) #For gradient function
#Use EM algorithm to fit two-polynomial constrained mixture density
width=midPoints[2]-midPoints[1]
nullFactor=dnorm(midPoints,mean=0,sd=sDelta)*width
null=pnorm(breakPoints,sd=sDelta)[-1]-pnorm(breakPoints,sd=sDelta)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
bool=(midPoints>0)
null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],sd=sDelta,lower.tail=FALSE)-pnorm(breakPoints[-1][bool],sd=sDelta,lower.tail=FALSE) 
null=nullFactor
EMout=EM(degree=degree_doublePolynomial,counts=counts,lambdaDeriv2=lambdaDeriv2,lambdaDeriv2m=lambdaDeriv2m,ij=ij,nullFactor=nullFactor,null=null,width=width,printMessage="negative effects, reversed sign")
par=par_NNPED_lower_rev=EMout$par;sumLambdaU=EMout$sumLambdaU;sumLambdaL=EMout$sumLambdaL
#Check fit
null=nullFactor #For visualisation purposes this is fine - multiply by width inside the following function
print("Printing graph of fitted two-polynomial mixture distribution")
fit<-predict.poissonRegMixture_twoPolynomial(par,figStem="poissonMixtureFit_twoPolynomial_lower_rev",n=nSim,width=width,midPoints=midPoints,scaleMeans=scaleMeans,scaleSDs=scaleSDs);fit[is.na(fit)]=0
twoPolynomial_fdr<-function(x,y){ #Needed for excluding lower power SNPs when plotting simulation truth, below
	p1=y[1];p2=y[2];b1=y[3:(degree+2)];b2=y[(degree+3):(2*degree+2)];c1=y[2*degree+3];c2=y[2*degree+4];c1a=y[2*degree+5];c2a=x[2*degree+6] 
	lambdaDeriv2=matrix(sapply(X=1:degree,y=x,FUN=make_poly),ncol=degree)
	lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
	lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
	lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
	ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
	polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
	polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	null=dnorm(x,mean=0,sd=sDelta);nullFactor=null
	lambdaU=exp(polyU)[,1]*nullFactor
	lambdaL=exp(polyL)[,1]*nullFactor
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	fdr=null*(1-p1-p2)/lambdaMix
	fdr[is.na(lambdaMix)|(lambdaMix==0)]=1 #Infinities can occasionally cause problems
	fdr[lambdaMix==Inf]=0
	return(fdr)
	}
#For each positive beta, find the probabilities that its deltas are negative, by looking them up on the lower distribution with Bayes' theorem. Then, assuming their effect direction was flipped but the size was the same, look their delta up on the upper distribution.
if(diffSampleSize>nSNP){
	diffSampleSize=nSNP
	message="Warning: diffSampleSize is greater than the number of SNPs. Setting to number of SNPs by default"
	print(message)
	outputLog=c(outputLog,message)
	}
print(paste0("Sampling ",diffSampleSize," effect size differences for each variant"))
#Return to using raw data for calculating the observed deltas
beta=fullData[,'beta'];s=fullData[,'s'];nSNP=nrow(fullData);fdr=fullData[,'fdr'];posProb=fullData[,'positiveEffectProb'];negProb=fullData[,'negativeEffectProb'];rm(fullData)
#Run bayes functions above on all SNPs, block by block
if(blockSize>=nSNP){
	blockSize=nSNP
	message="Warning: blockSize is greater or equal to the number of SNPs. Proceeding using one single block of all SNPs"
	print(message);outputLog=c(outputLog,message)
		}
snpBlocks=as.list(data.frame(matrix(1:(blockSize*floor(nSNP/blockSize)),ncol=floor(nSNP/blockSize),byrow=FALSE)))
if(max(unlist(snpBlocks))<nSNP){
	remainders=(max(unlist(snpBlocks))+1):nSNP
	snpBlocks=c(snpBlocks,list(remainders))
		}
threadPlural="parallel threads";if(threads==1)threadPlural="thread"
print(paste0("Computing posteriors for ",nSNP*diffSampleSize," effect size differences in blocks of ",blockSize,", using ",threads," ",threadPlural,". Can take more than about ",ceiling(20/threads)," minutes"))
index=sample(rep((1:nSNP),rmultinom(1,size=nSim,prob=(1-fdr)/((1-pNull)*nSNP))[,1]),replace=FALSE,size=nSim) #Sample of observed SNPs as if they were from alt distribution, as before when choosing alt SNPs for plotting
bayesLookup<-function(x,y,direction="positive",beta,betaRand,s,sRand,probRand,deltaRange,scaleMeans,scaleSDs,printMessage=NULL,diffSampleSize,includeEquality=TRUE){
	set.seed(seed+x[1])
	print(paste0("Computing ",direction," posteriors for SNPs ",x[1],"-",x[length(x)],", ",printMessage))
	degree=(length(y)-6)/2
	s1=y[1];s2=y[2];b1=y[3:(degree+2)];b2=y[(degree+3):(2*degree+2)];c1=y[2*degree+3];c2=y[2*degree+4];c1a=y[2*degree+5];c2a=y[2*degree+6] 
    	p1=0.5*exp(s1)*(s1<0)+(1-0.5*exp(-s1))*(s1>=0)
    	p2=0.5*exp(s2)*(s2<0)+(1-0.5*exp(-s2))*(s2>=0)
	#x is the 'block' parameter
	rand=sample(1:length(betaRand),replace=FALSE,size=diffSampleSize)
	betaRand2=betaRand[rand];sRand2=sRand[rand];probRand2=probRand[rand]
	delta<-findDelta(x=cbind(beta[x],s[x]),betaRand2=betaRand2,sRand2=sRand2,deltaRange=deltaRange)
	lambdaDeriv2=matrix(sapply(X=1:degree,y=delta,FUN=make_poly),ncol=degree)
	lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
	lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
	lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
	ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
	lambdaU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1 #Save memory by calling these lambda rather than 'poly'
	lambdaL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	rm(lambdaDeriv2);rm(lambdaDeriv2m)
	null=dnorm(as.vector(delta),mean=0,sd=sDelta);nullFactor=null
	lambdaU=exp(lambdaU)[,1]*nullFactor
	lambdaL=exp(lambdaL)[,1]*nullFactor
	lambdaU[lambdaU==Inf]=nSim;lambdaL[lambdaL==Inf]=nSim
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	bayes=(direction=="positive")*lambdaU*p2/lambdaMix + (direction=="negative")*lambdaL*p1/lambdaMix + includeEquality*null*(1-p1-p2)/lambdaMix #Null is included in the upper lookup to allow for BDR=1
	bayesStore=probRand2%*%matrix(bayes,ncol=length(x),byrow=FALSE)/diffSampleSize
	x=x[bayesStore<highAccThresh]
	if(length(x)>0){
	rand=sample(1:length(betaRand),replace=FALSE,size=diffSampleSizeHighAcc)
	betaRand2=betaRand[rand];sRand2=sRand[rand];probRand2=probRand[rand]
	delta<-findDelta(x=cbind(beta[x],s[x]),betaRand2=betaRand2,sRand2=sRand2,deltaRange=deltaRange)
	lambdaDeriv2=matrix(sapply(X=1:degree,y=delta,FUN=make_poly),ncol=degree)
	lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
	lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
	lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
	ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
	lambdaU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1 #Save memory by calling these lambda rather than 'poly'
	lambdaL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	rm(lambdaDeriv2);rm(lambdaDeriv2m)
	null=dnorm(as.vector(delta),mean=0,sd=sDelta);nullFactor=null
	lambdaU=exp(lambdaU)[,1]*nullFactor
	lambdaL=exp(lambdaL)[,1]*nullFactor
	lambdaU[lambdaU==Inf]=nSim;lambdaL[lambdaL==Inf]=nSim
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	bayes=(direction=="positive")*lambdaU*p2/lambdaMix + (direction=="negative")*lambdaL*p1/lambdaMix + includeEquality*null*(1-p1-p2)/lambdaMix #Null is included in the upper lookup to allow for BDR=1
	bayesStore[bayesStore<highAccThresh]=(probRand2%*%matrix(bayes,ncol=length(x),byrow=FALSE)+bayesStore[bayesStore<highAccThresh]*diffSampleSize)/(diffSampleSize+diffSampleSizeHighAcc) 
		}
	return(bayesStore)
		}
BDR_U=unlist(mclapply(FUN=bayesLookup,X=snpBlocks,y=par_NNPED_upper,direction="negative",beta=beta,betaRand=beta[index],s=s,sRand=s[index],probRand=rep(1,length(index)),deltaRange=deltaRange_U,scaleMeans=rep(0,degree),scaleSDs=scaleList$NUPEd_upper,includeEquality=TRUE,printMessage="assuming original effect signs",diffSampleSize=diffSampleSize,mc.cores=threads,mc.preschedule=FALSE))
BDR_U_rev=unlist(mclapply(FUN=bayesLookup,X=snpBlocks,y=par_NNPED_upper_rev,direction="negative",beta=beta,betaRand=-beta[index],s=s,sRand=s[index],probRand=rep(1,length(index)),deltaRange=deltaRange_U_rev,scaleMeans=rep(0,degree),scaleSDs=scaleList$NUPEd_upper_rev,includeEquality=FALSE,printMessage="assuming flipped effect signs",diffSampleSize=diffSampleSize,mc.cores=threads,mc.preschedule=FALSE))
BDR_L=unlist(mclapply(FUN=bayesLookup,X=snpBlocks,y=par_NNPED_lower,direction="positive",beta=beta,betaRand=beta[index],s=s,sRand=s[index],probRand=rep(1,length(index)),deltaRange=deltaRange_L,scaleMeans=rep(0,degree),scaleSDs=scaleList$NUPEd_lower,includeEquality=TRUE,printMessage="assuming original effect signs",diffSampleSize=diffSampleSize,mc.cores=threads,mc.preschedule=FALSE))
BDR_L_rev=unlist(mclapply(FUN=bayesLookup,X=snpBlocks,y=par_NNPED_lower_rev,direction="positive",beta=beta,betaRand=-beta[index],s=s,sRand=s[index],probRand=rep(1,length(index)),deltaRange=deltaRange_L_rev,scaleMeans=rep(0,degree),scaleSDs=scaleList$NUPEd_lower_rev,includeEquality=FALSE,printMessage="assuming flipped effect signs",diffSampleSize=diffSampleSize,mc.cores=threads,mc.preschedule=FALSE))
BDRs=(BDR_U+BDR_U_rev)*posProb+(BDR_L+BDR_L_rev)*negProb+(1-posProb-negProb)
BFDR=fdr+BDRs*(1-fdr)
if(sum(BFDR>1)>0)outputLog=c(outputLog,paste0("Warning: ",sum(BFDR>1)," variables with local BFDR>1, ranging ",round(range(BFDR[BFDR>1])[1],4),"-",round(range(BFDR[BFDR>1])[2],4),". Setting these to 1."))
BFDR[BFDR>1]=1
print("Finding tail-area rates")
cuMean=cumsum(sort(BDRs,decreasing=FALSE))/seq_along(BDRs)
tailStats=cuMean[rank(BDRs,ties.method="max")]
cuMean=cumsum(sort(fdr,decreasing=FALSE))/seq_along(fdr)
tailStats=cbind(cuMean[rank(fdr,ties.method="max")],tailStats)
tailBFDRs=tailStats[,1]+tailStats[,2]*(1-tailStats[,1])
print("BFDR.priorsplitteR analysis complete.")
outputLog=c(outputLog,"BFDR.priorsplitteR analysis complete.")
write.table(outputLog,file=paste0(outFileStem,"log.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
return(list(results=data.frame(localBFDR=BFDR,localFDR=fdr,localBDR=BDRs,tailBFDR=tailBFDRs,tailFDR=tailStats[,1],tailBDR=tailStats[,2]),ANEfit=par_ANE,NNPEDfit=cbind(upper=par_NNPED_upper,upper_reverse=par_NNPED_upper_rev,lower=par_NNPED_lower,lower_reverse=par_NNPED_lower_rev),NNPEfit=par_NNPE,scaleFactors=scaleList,outputLog=outputLog))
			}#End of BFDR function


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
	lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
	lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
	ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
	lambdaU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1 #Save memory by calling these lambda rather than 'poly'
	lambdaL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	rm(lambdaDeriv2);rm(lambdaDeriv2m)
	null=dnorm(as.vector(delta),mean=0,sd=sDelta);nullFactor=null
	lambdaU=exp(lambdaU)[,1]*nullFactor
	lambdaL=exp(lambdaL)[,1]*nullFactor
	lambdaU[lambdaU==Inf]=nSim;lambdaL[lambdaL==Inf]=nSim
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	bayes=c(lambdaU*p2/lambdaMix,null*(1-p1-p2)/lambdaMix,lambdaL*p1/lambdaMix) #Null is included in the upper lookup to allow for BDR=1
	return(bayes)
		}

