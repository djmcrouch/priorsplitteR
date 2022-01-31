BFDR.priorsplitteR<-function(beta,s,seed,excludeFromFit=NULL,simulationTruth=NULL,outFileStem,breaks=100,degree_altNull=6,degree_doublePolynomial=6,nSim=5000000,diffSampleSize=1000,diffSampleSizeHighAcc=9000,highAccThresh=0.05,blockSize=10000,EMtol=c(1e-15,1e-15),nStop=5,maxEMit=c(100000,100000),ml.reltol=1e-300,outlierSD=15,diffOutlierSD=NULL,threads=1){
library(distr);library(parallel);library(dplyr)
set.seed(seed)
outputLog=paste("outFileStem=",outFileStem,"::breaks=",breaks,"::nSim=",nSim,"::diffSampleSize=",diffSampleSize,"::degree_altNull=",degree_altNull,"::degree_doublePolynomial=",degree_doublePolynomial,"::blockSize=",blockSize,"::EMtol=",EMtol,"::nStop=",nStop,"::ml.reltol=",ml.reltol,"::outlierSD=",outlierSD,"::threads=",threads)
if((nchar(outFileStem)>0)&(substr(outFileStem,nchar(outFileStem),nchar(outFileStem))!="_"))outFileStem=paste0(outFileStem,"_")
z=beta/s ; nSNP=length(beta)
fullData=data.frame(beta=beta,s=s,z=z) #Store full data before removing outliers
#For outliers, can't rely on the polynomials outside the window, so set the outlier z-scores to the outlier threshold (w.r.t null distribution), and recalculate their betas from their s and new z
indexNonOutliers=(1:nSNP)[abs(z)<=outlierSD]
indexNonOutliers_nonExclude=indexNonOutliers[!indexNonOutliers%in%excludeFromFit]
z[-indexNonOutliers]=range(z[indexNonOutliers_nonExclude])[1*(sign(z[-indexNonOutliers])>0)+1]
beta[-indexNonOutliers]=z[-indexNonOutliers]*s[-indexNonOutliers]
#Remove outliers and excluded SNPs only after storing full data
fullData=data.frame(fullData,beta_squashed=beta,z_squashed=z) #Store data before removing outliers
outputLog=c(outputLog,paste0("There are ",nSNP-length(indexNonOutliers)," out of ",nSNP," variables exceeding the outlier Z-score threshold of ",outlierSD,". Removing these and proceeding with ",length(indexNonOutliers)," variables."))
z=z[indexNonOutliers_nonExclude];beta=beta[indexNonOutliers_nonExclude];s=s[indexNonOutliers_nonExclude];nSNP=length(indexNonOutliers_nonExclude);simulationTruth=simulationTruth[indexNonOutliers_nonExclude]
#Simulate Z-score difference data from the null - we'll use this simulated data for finding the pairwise difference distribution. Use Poisson mixture regression with one distribution fixed as null, to fit the two groups
make_poly<-function(x,y){ #Function to make polynomials
	y^x
		}
by=round(diff(range(z))/breaks,digits=2)
breakPoints=sort(c(seq(-by/2,min(z)-by,by=-by),seq(by/2,max(z)+by,by=by))) 
histogram<-hist(z,breaks=breakPoints,plot=F) #Call histogram but don't plot
midPoints<-histogram$mids
counts<-histogram$counts
degree<-degree_altNull
lambdaDeriv1=cbind(1,sapply(X=1:degree,y=midPoints,FUN=make_poly))
scaleMeans=colMeans(lambdaDeriv1);scaleSDs=apply(lambdaDeriv1,FUN=sd,MAR=2)
lambdaDeriv1[,-1]=scale(lambdaDeriv1[,-1])
colnames(lambdaDeriv1)=NULL
width=midPoints[2]-midPoints[1]
null=pnorm(breakPoints)[-1]-pnorm(breakPoints)[-length(breakPoints)]
#Fix numerical problems by taking the upper rather than lower cdf. Don't want any zeros here for when we take null probs
bool=(midPoints>0)
null[bool]=pnorm(breakPoints[-length(breakPoints)][bool],lower.tail=FALSE)-pnorm(breakPoints[-1][bool],lower.tail=FALSE) 
lNull=log(null) #Helps avoid numerical issues in function below
#Use EM algorithm to fit two group mixture density
likelihood.poissonRegMixture_altNull<-function(x){
	p=x[1];b=x[-1]
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
	p=x[1];b=x[-1]
	lambda=exp(lambdaDeriv1%*%b)[,1]
	deriv=c(sum(lambda*nSNP-counts*(1-m)/(1-p)+counts*m/p-null*nSNP),(counts*(1-m))%*%lambdaDeriv1-(lambda*nSNP*(1-p))%*%lambdaDeriv1)
	return(-deriv)
	}
m<-rep(0.5,length(counts))   #Distribution membership variables
i=1;breakCounter=0;par=c(0.99,lNull[midPoints==0]-1e-5,rep(0,degree))
repeat{
	print(paste0("Alt-Null EM iteration=",i))
	ui=matrix(0,nrow=3,ncol=degree+2);ui[1,1]=1;ui[2,1]=-1;ci=c(0,-1,-lNull[midPoints==0]+1e-10) #Intercept Less than log null minus small constant
	ui[3,2]=-1
mixFit=constrOptim(theta=par,f=likelihood.poissonRegMixture_altNull,grad=grad.poissonRegMixture_altNull,method="BFGS",control=list(maxit=100000000,reltol=ml.reltol),ui=ui,ci=ci)
	x=par
	p=x[1];b=x[-1]
	poly=lambdaDeriv1%*%b
	lambda=exp(poly)[,1]
	l1=counts*(1-m)*(log(1-p)+poly)-(1-p)*lambda*nSNP
	l2=counts*m*(log(p))-null*p*nSNP
	l1[(l1==(-Inf))|is.na(l1)]=-1e300
	l2[(l2==(-Inf))|is.na(l2)]=-1e300
	Qjj=sum(l1+l2) #Don't need to worry about constants in the likelihoods - see tibshirani, hastie and friedman 'elements of statistical learning' EM algorithm section p276/277, equation 8.46 and following short paragraphs of text. The goal is find a value s.t. Q(\theta^{j+1},\theta^j) > Q(\theta^{j},\theta^j), and the denominator is the same in both these terms.
	x=mixFit$par #Find core expected likelihood given new parameters
	p=x[1];b=x[-1]
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
lambdaDeriv1=sapply(X=1:degree,y=x,FUN=make_poly)
lambdaDeriv1=scale(lambdaDeriv1,center=scaleMeans[-1],scale=scaleSDs[-1])
lambdaDeriv1=cbind(1,lambdaDeriv1)
null=dnorm(x,mean=0,1)*width
predict.poissonRegMixture_altNull<-function(x){
	p=x[1];b=x[-1]
	lambda=exp(lambdaDeriv1%*%b)
	lambdaMix=lambda*(1-p)+null*p
	return(nSNP*cbind(lambdaMix,lambda,null,lambda*(1-p),null*p))
		}
fit<-predict.poissonRegMixture_altNull(par);fit[is.na(fit)]=0
print("Printing graph of fitted alt/null mixture distribution")
#Put error bars on this fit
tiff(paste0(outFileStem,"poissonMixtureFit_altNull.tiff"))
plot(histogram,main=NULL,ylim=c(0,max(fit))*1.05)
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,2],col=3,type="l",lwd=2)
points(x,fit[,3],col=4,type="l",lwd=2)
dev.off()
tiff(paste0(outFileStem,"poissonMixtureFit_altNull_mixingProp.tiff"))
plot(histogram,main=NULL,ylim=c(0,max(fit))*1.05)
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,4],col=3,type="l",lwd=2)
points(x,fit[,5],col=4,type="l",lwd=2)
dev.off()
#Now simulate from the alternative distribution
parAltNull=par
pNull=par[1];b=par[-1]
altNull_Density<-function(x){
	lambdaDeriv1=matrix(sapply(X=1:degree,y=x,FUN=make_poly),nrow=length(x),byrow=FALSE)
	lambdaDeriv1=scale(lambdaDeriv1,center=scaleMeans[-1],scale=scaleSDs[-1])
	lambdaDeriv1=cbind(1,matrix(lambdaDeriv1,nrow=nrow(lambdaDeriv1)))
	lambda=exp(lambdaDeriv1%*%b)/width
	lambda[is.na(lambda)]=0 #Infinities can occasionally cause problems
	return(as.vector(lambda))
	}
fullData=data.frame(fullData,fdr=dnorm(fullData[,'z_squashed'])*pNull/(altNull_Density(fullData[,'z_squashed'])*(1-pNull)+dnorm(fullData[,'z_squashed'])*pNull)) #Store for later
print(paste0("Simulating ",nSim," variants from estimated distribution function. Printing histogram of simulation"))
#Following few lines from 'https://stackoverflow.com/questions/23570952/simulate-from-an-arbitrary-continuous-probability-distribution'
dist <-AbscontDistribution(d=altNull_Density,low1=min(midPoints),up1=max(midPoints)) #signature for a dist with pdf ~ twoGroupDensity
rdist <- r(dist) # function to create random variates 
sim=rdist(nSim)
tiff(paste0(outFileStem,"poissonMixtureFit_altNull_altSim.tiff"))
x<-seq(min(midPoints),max(midPoints),by=0.01)
y<-unlist(sapply(X=x,FUN=altNull_Density))
yNull=dnorm(x,mean=0,sd=1)
hist(sim,freq=FALSE,breaks=breakPoints,col=0,ylim=c(0,max(c(y,yNull))),main=NULL,xlab="",xlim=range(z))
lines(x,y,type='l',col="green")
lines(x,yNull,type='l',col="blue")
dev.off()
#Simulate differences by drawing nSim further variables.
fdr=fullData[,'fdr']
index1=rep(1:length(fdr),rmultinom(1,size=nSim,prob=(1-fdr)/((1-pNull)*length(fdr)))[,1]) #Sample of observed SNPs as if they were from alt #distribution, as before when choosing alt SNPs for plotting
index2=rep(1:length(fdr),rmultinom(1,size=nSim,prob=(1-fdr)/((1-pNull)*length(fdr)))[,1]) 
s1=fullData[index1,'s'];s2=fullData[index2,'s']
delta=(sim*s1-rdist(nSim)*s2)/sqrt(s1^2+s2^2)
if(!is.null(simulationTruth)){ 
	#This section estimates the distribution of effect size differences for known alt-SNPs (using truth information e.g. from a simulation). This is used later in the code, and a plot is produced which compares the differences drawn from the estimated distribution of alt-SNP differences, with the equivalent differences estimated from the simulation truth.
	index=sample((1:nSNP)[simulationTruth!=0],size=nSim,replace=TRUE)
	rand=sample(index,size=length(index),replace=TRUE)
	rand[rand==index]=sample(rand[rand!=index],size=sum(rand==index),replace=FALSE) #Replace ones where index and rand are the same, as this introduces a bias towards 'null' differences
	deltaTrue=(beta[index]-beta[rand])/sqrt(s[index]^2+s[rand]^2)
	trueDiff=(simulationTruth[index]-simulationTruth[rand])/sqrt(s[index]^2+s[rand]^2)
	tiff(paste0(outFileStem,"poissonMixtureFit_altNull_diffs_simVersusReal.tiff"))
	deltaHist=hist(c(delta,deltaTrue),breaks=breaks)
	deltaHist1=hist(delta,breaks=deltaHist$breaks);deltaHist2=hist(deltaTrue,breaks=deltaHist$breaks) #For count range
	hist(delta,breaks=deltaHist$breaks,col=rgb(0,0,1,alpha=0.5),main=NULL,xlab="Effect size differences",xlim=range(c(delta,deltaTrue)),ylim=1.2*c(0,max(c(deltaHist1$counts,deltaHist2$counts))),ylab="Counts")
	hist(deltaTrue,breaks=deltaHist$breaks,add=T,col=rgb(1,0,0,alpha=0.5))	
	legend(x=min(c(delta,deltaTrue)),y=1.2*max(c(deltaHist1$counts,deltaHist2$counts)),col=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5)),legend=c("Drawn directly from estimated alt distribution","Simulation truth alt differences"),pch=15)
	dev.off()
	deltaTrue_simTruth=deltaTrue;trueDiff_simTruth=trueDiff #Store for later
#Repeat the above plot inside the condition, but this time using a sample of *presumed* alt SNPs - i.e. sample SNPs according to their power. A measure of power is the expected FDR for each true effect size (see Efron 2007: "Size, power and and false discovery rates"). 
index=rep((1:nSNP),rmultinom(1,size=nSim,prob=(1-fullData[indexNonOutliers,'fdr'])/((1-pNull)*nSNP))[,1])
rand=rep(1:nSNP,rmultinom(1,size=nSim,prob=(1-fullData[indexNonOutliers,'fdr'])/((1-pNull)*nSNP))[,1])
rand=sample(rand,size=length(rand),replace=FALSE) #Shuffle order
rand[rand==index]=sample(rand[rand!=index],size=sum(rand==index),replace=FALSE) #Replace ones where index and rand are the same, as this introduces a bias towards 'null' differences
deltaTrue=(beta[index]-beta[rand])/sqrt(s[index]^2+s[rand]^2)
trueDiff=(simulationTruth[index]-simulationTruth[rand])/sqrt(s[index]^2+s[rand]^2)
tiff(paste0(outFileStem,"poissonMixtureFit_altNull_diffs_simVersusEstAlt.tiff"))
deltaHist=hist(c(delta,deltaTrue),breaks=breaks,plot=FALSE)
deltaHist1=hist(delta,breaks=deltaHist$breaks,plot=FALSE);deltaHist2=hist(deltaTrue,breaks=deltaHist$breaks,plot=FALSE) #For count range
hist(delta,breaks=deltaHist$breaks,col=rgb(0,0,1,alpha=0.5),main=NULL,xlab="Effect size differences",xlim=range(c(delta,deltaTrue)),ylim=1.2*c(0,max(c(deltaHist1$counts,deltaHist2$counts))),ylab="Counts")
hist(deltaTrue,breaks=deltaHist$breaks,add=T,col=rgb(1,0,0,alpha=0.5))	
legend(x=min(c(delta,deltaTrue)),y=1.2*max(c(deltaHist1$counts,deltaHist2$counts)),col=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5)),legend=c("Drawn directly from estimated alt distribution","Empirically estimated alt distribution based on fdrs"),pch=15)
dev.off()
deltaTrue_altEst=deltaTrue;trueDiff_altEst=trueDiff #Store for later #Store for later
		}
#Fit the 2 polynomial model, again with poisson regression. First remove any extreme outliers that exist in the deltas - should't be many
if(is.null(diffOutlierSD))diffOutlierSD=outlierSD*1.5
deltaMax=diffOutlierSD 
bool=(abs(delta)>=deltaMax)
delta=delta[!bool];nSim=length(delta)
deltaRange=range(delta) #Will also use this again later in Bayes theorem computations
if(sum(bool)>0){message=paste0("Warning: there are ",sum(bool)," extreme outlier deltas (out of total ",nSim,") with absolute values greater than than ",diffOutlierSD,". Removing these. Consider lowering the outlierSD parameter if the number of these seems high.")}else{message=NULL}
outputLog=c(outputLog,message);if(!is.null(message))print(message)
sDelta=1 #sds of null differences
by=round(diff(range(delta))/breaks,digits=2)
breakPoints=sort(c(seq(-by/2,min(delta)-by,by=-by),seq(by/2,max(delta)+by,by=by))) 
histogram<-hist(delta,breaks=breakPoints,plot=F) #Call histogram but don't plot
midPoints<-histogram$mids
counts<-histogram$counts
degree<-degree_doublePolynomial
lambdaDeriv2=sapply(X=1:degree,y=midPoints,FUN=make_poly)
scaleMeans=colMeans(lambdaDeriv2);scaleSDs=apply(lambdaDeriv2,FUN=sd,MAR=2)
lambdaDeriv2=scale(lambdaDeriv2) #Don't center first term
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
lNull=log(null) #Helps avoid numerical issues in function below
likelihood.poissonRegMixture_twoPolynomial<-function(x){
	p1=x[1];p2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
	polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
	polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	lambdaU=exp(polyU)[,1]*nullFactor
	lambdaL=exp(polyL)[,1]*nullFactor
	l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*nSim
	l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*nSim
	l3=counts*(1-m1-m2)*(log(1-p1-p2))-(1-p1-p2)*null*nSim
	l1[(l1==(-Inf))|is.na(l1)]=-1e300
	l2[(l2==(-Inf))|is.na(l2)]=-1e300
	l3[(l3==(-Inf))|is.na(l3)]=-1e300
	l=sum(l1+l2+l3)
	return(-l)
	}
grad.poissonRegMixture_twoPolynomial<-function(x){
	p1=x[1];p2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
	polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
	polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	lambdaU=exp(polyU)[,1]*nullFactor
	lambdaL=exp(polyL)[,1]*nullFactor
	coeffDerivU=2*cbind(1,lambdaDeriv2)*cbind(1,lambdaDeriv2)%*%(matrix(rep(c(c1a,b1),degree+1),ncol=degree+1,byrow=FALSE)/ij)
	coeffDerivL=2*cbind(1,lambdaDeriv2)*cbind(1,lambdaDeriv2)%*%(matrix(rep(c(c2a,b2),degree+1),ncol=degree+1,byrow=FALSE)/ij)
	deriv=c(sum(counts*m1/p1-counts*(1-m1-m2)/(1-p1-p2)-lambdaL*nSim+null*nSim),sum(counts*m2/p2-counts*(1-m1-m2)/(1-p1-p2)-lambdaU*nSim+null*nSim),((counts*m2-p2*lambdaU*nSim)*lambdaDeriv2[,1])%*%coeffDerivU[,-1],-((counts*m1-p1*lambdaL*nSim)*lambdaDeriv2[,1])%*%coeffDerivL[,-1],sum(counts*m2-p2*lambdaU*nSim),sum(counts*m1-p1*lambdaL*nSim),((counts*m2-p2*lambdaU*nSim)*lambdaDeriv2[,1])%*%coeffDerivU[,1],-((counts*m1-p1*lambdaL*nSim)*lambdaDeriv2[,1])%*%coeffDerivL[,1])
	return(-deriv)
	}
m1<-rep(0.33,length(counts))   #Lower tail responsibilities 
m2<-rep(0.33,length(counts))   #Upper tail responsibilities 
i=1;breakCounter=0;par=c(0.333,0.333,rep(0,degree),rep(0,degree),-1e-5,-1e-5,1e-5,1e-5)
repeat{
	print(paste0("Two-polynomial EM iteration=",i))
	ui=matrix(0,nrow=9,ncol=degree*2+6);ci=c(0,-1,0,-1,-1,0,0,0,0)
	ui[1:2,1]=c(1,-1) #Mixing proportion constraints - between pMin and 1
	ui[3:4,2]=c(1,-1) #Mixing proportion constraints - between pMax and 0
	ui[5,1:2]=-1 #Sum of mixing props <1 
	diag(ui[6:7,degree*2+3:4])=c(-1,-1) #Intercept at 0 constraints - ratio with null should be < 1
	diag(ui[8:9,2*degree+5:6])=c(1,1) #Second intercept should be >0 to stop potential symmetry/identifiability issues
	par_tweak=par #Modify par to avoid boundary constraint issues
	par_tweak[2*degree+3:4][par[2*degree+3:4]>(-1e-2)]=-1e-2#Initialise these two parameters comfortably within the boundary region
	par_tweak[2*degree+5:6][par[2*degree+5:6]<(1e-2)]=1e-2#Initialise these two parameters comfortably within the boundary region
	if(sum(par[1:2])<1e-2)par_tweak[1:2]=par_tweak[1:2]+1e-2
	if(sum(par[1:2])>0.99)par_tweak[1:2]=par_tweak[1:2]-1e-2
	par_tweak[1:2][par[1:2]<1e-2]=1e-2#Initialise these two parameters comfortably within the boundary region
	par_tweak[1:2][par[1:2]>0.99]=0.99#Initialise these two parameters comfortably within the boundary region
	repeat{ 
		#Very occassionally constrOptim will fail at certain initial values of par. This loop catches errors and restarts from a new initialisation point after adding some random noise
	mixFit=try(constrOptim(theta=par_tweak,f=likelihood.poissonRegMixture_twoPolynomial,grad=grad.poissonRegMixture_twoPolynomial,method="BFGS",control=list(maxit=100000000,reltol=ml.reltol),ui=ui,ci=ci),silent=TRUE)
	if(!is.character(mixFit[1]))break
	par_tweak=par_tweak+rnorm(length(par_tweak),sd=0.01)
	par_tweak[1:2][!between(par_tweak[1:2],0.01,0.99)]=0.33
	if(sum(par_tweak[1:2])>0.99)par_tweak[1:2]=0.33	
	par_tweak[2*degree+3:4][par_tweak[2*degree+3:4]>(-1e-5)]=-1e-05 
	par_tweak[2*degree+5:6][par_tweak[2*degree+5:6]<1e-5]=1e-05 
	print("Error in likelihood maximisation. Adding random noise to starting parameters and retrying")
		}
	x=par
	p1=x[1];p2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] 
	polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
	polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	lambdaU=exp(polyU)[,1]*nullFactor
	lambdaL=exp(polyL)[,1]*nullFactor
	l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*nSim
	l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*nSim
	l3=counts*(1-m1-m2)*(log(1-p1-p2))-(1-p1-p2)*null*nSim
	l1[(l1==(-Inf))|is.na(l1)]=-1e300
	l2[(l2==(-Inf))|is.na(l2)]=-1e300
	l3[(l3==(-Inf))|is.na(l3)]=-1e300
	Qjj=sum(l1+l2+l3)
	x=mixFit$par  #Find core expected likelihood with given the parameters
	p1=x[1];p2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] 
	polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
	polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	lambdaU=exp(polyU)[,1]*nullFactor
	lambdaL=exp(polyL)[,1]*nullFactor
	l1=counts*m1*(log(p1)+polyL)-p1*lambdaL*nSim
	l2=counts*m2*(log(p2)+polyU)-p2*lambdaU*nSim
	l3=counts*(1-m1-m2)*(log(1-p1-p2))-(1-p1-p2)*null*nSim
	l1[(l1==(-Inf))|is.na(l1)]=-1e300
	l2[(l2==(-Inf))|is.na(l2)]=-1e300
	l3[(l3==(-Inf))|is.na(l3)]=-1e300
	QjjPlus1=sum(l1+l2+l3)
	m1<-lambdaL*p1/(lambdaU*p2+lambdaL*p1+null*(1-p1-p2))
	m2<-lambdaU*p2/(lambdaU*p2+lambdaL*p1+null*(1-p1-p2))
	m1[(lambdaU*p2+lambdaL*p1+null*(1-p1-p2))==0]=0	#Fix numerical probs with zero denominators 
	m2[(lambdaU*p2+lambdaL*p1+null*(1-p1-p2))==0]=0
	par=mixFit$par
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
#Check fit
x=seq(min(midPoints),max(midPoints),by=diff(range(midPoints))/100000)
lambdaDeriv2=sapply(X=1:degree,y=x,FUN=make_poly)
lambdaDeriv2=scale(lambdaDeriv2,center=scaleMeans,scale=scaleSDs)
lambdaDeriv2m=apply(X=cbind(1,lambdaDeriv2),FUN=make_lambdaDeriv2m,MARGIN=2,y=cbind(1,lambdaDeriv2))
lambdaDeriv2m=matrix(t(lambdaDeriv2m),ncol=(degree+1)^2,byrow=TRUE)
ij=(1+matrix(rep(rep(0:degree,degree+1)+sort(rep(0:degree,degree+1),decreasing=FALSE),nrow(lambdaDeriv2)),ncol=(degree+1)^2,byrow=TRUE,nrow=nrow(lambdaDeriv2)));lambdaDeriv2m=lambdaDeriv2m/ij
nullFactor=dnorm(x,mean=0,sd=sDelta)*width
null=nullFactor #For visualisation purposes this is fine - multiply by width inside the following function
predict.poissonRegMixture_twoPolynomial<-function(x){
	p1=x[1];p2=x[2];b1=x[3:(degree+2)];b2=x[(degree+3):(2*degree+2)];c1=x[2*degree+3];c2=x[2*degree+4];c1a=x[2*degree+5];c2a=x[2*degree+6] #c is an 'extra' intercept term corresponding to the constant of integration for the upper distribution, to allow it to be lower than null at x=0
	polyU=(lambdaDeriv2m%*%matrix(c(c1a,b1)%*%t(c(c1a,b1)),ncol=1))*lambdaDeriv2[,1]+c1
	polyL=-(lambdaDeriv2m%*%matrix(c(c2a,b2)%*%t(c(c2a,b2)),ncol=1))*lambdaDeriv2[,1]+c2
	lambdaU=exp(polyU)[,1]*nullFactor
	lambdaL=exp(polyL)[,1]*nullFactor
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	return(nSim*cbind(lambdaMix,length(lambdaL)/length(midPoints)*lambdaL/sum(lambdaL),length(lambdaU)/length(midPoints)*lambdaU/sum(lambdaU),null,lambdaL*p1,lambdaU*p2,null*(1-p1-p2)))
		}
print("Printing graph of fitted two-polynomial mixture distribution")
fit<-predict.poissonRegMixture_twoPolynomial(par);fit[is.na(fit)]=0
tiff(paste0(outFileStem,"poissonMixtureFit_twoPolynomial.tiff"))
plot(histogram,main=NULL,ylim=c(0,max(fit[,-(5:7)]))*1.05)
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,2],col=3,type="l",lwd=2)
points(x,fit[,3],col=5,type="l",lwd=2)
points(x,fit[,4],col=4,type="l",lwd=2)
dev.off()
tiff(paste0(outFileStem,"poissonMixtureFit_twoPolynomial_mixingPropScaled.tiff"))
plot(histogram,main=NULL,ylim=c(0,max(fit[,-(2:4)]))*1.05)
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,5],col=3,type="l",lwd=2)
points(x,fit[,6],col=5,type="l",lwd=2)
points(x,fit[,7],col=4,type="l",lwd=2)
dev.off()
twoPolynomial_fdr<-function(x,y){ #Needed for excluding lower power SNPs when plotting simulation truth, below
	p1=y[1];p2=y[2];b1=y[3:(degree+2)];b2=y[(degree+3):(2*degree+2)];c1=y[2*degree+3];c2=y[2*degree+4];c1a=y[2*degree+5];c2a=x[2*degree+6] 
	lambdaDeriv2=sapply(X=1:degree,y=x,FUN=make_poly)
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
if(!is.null(simulationTruth)){ 
deltaTrue=deltaTrue_simTruth;trueDiff=trueDiff_simTruth #Store for later
keep=((deltaTrue<=max(histogram$breaks))&((deltaTrue>=min(histogram$breaks))))  #Avoids graphing problems
deltaTrue=deltaTrue[keep];trueDiff=trueDiff[keep]
deltaHist1=hist(delta,breaks=histogram$breaks,plot=FALSE);deltaHist2=hist(deltaTrue[trueDiff<0],breaks=histogram$breaks,plot=FALSE);deltaHist3=hist(deltaTrue[trueDiff>0],breaks=histogram$breaks,plot=FALSE) #For count range
tiff(paste0(outFileStem,"poissonMixtureFit_twoPolynomial_mixingPropScaled_withSimAltDiffs.tiff"),width=1000,height=1000,res=200)
plot(histogram,main=NULL,ylim=1.2*c(0,max(c(deltaHist1$counts,deltaHist2$counts,deltaHist3$counts))))
hist(deltaTrue[trueDiff>0],breaks=histogram$breaks,col=rgb(1,0,0,alpha=0.5),add=T)
hist(deltaTrue[trueDiff<0],breaks=histogram$breaks,col=rgb(0,0,1,alpha=0.5),add=T)
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,5],col=3,type="l",lwd=2)
points(x,fit[,6],col=5,type="l",lwd=2)
points(x,fit[,7],col=4,type="l",lwd=2)
legend(x=min(midPoints),y=1.2*max(c(deltaHist1$counts,deltaHist2$counts)),col=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5)),legend=c("True effect < 0","True effect > 0"),pch=15,cex=0.6)
dev.off()
deltaTrue=deltaTrue_altEst;trueDiff=trueDiff_altEst
rand=1:nSim#rep(1:nSim,rmultinom(1,size=nSim,prob=(1-twoPolynomial_fdr(x=deltaTrue,y=par))/((par[1]+par[2])*nSim))[,1])
deltaTrue=deltaTrue[rand];trueDiff=trueDiff[rand]
keep=((deltaTrue<=max(histogram$breaks))&((deltaTrue>=min(histogram$breaks))))  #Avoids graphing problems
deltaTrue=deltaTrue[keep];trueDiff=trueDiff[keep]
deltaHist1=hist(delta,breaks=histogram$breaks,plot=FALSE);deltaHist2=hist(deltaTrue[trueDiff<0],breaks=histogram$breaks,plot=FALSE);deltaHist3=hist(deltaTrue[trueDiff>0],breaks=histogram$breaks,plot=FALSE) #For count range
tiff(paste0(outFileStem,"poissonMixtureFit_twoPolynomial_mixingPropScaled_withEstAltDiffs.tiff"),width=1000,height=1000,res=200)
plot(histogram,main=NULL,ylim=1.2*c(0,max(c(deltaHist1$counts,deltaHist2$counts,deltaHist3$counts))))
hist(deltaTrue[trueDiff>0],breaks=histogram$breaks,col=rgb(1,0,0,alpha=0.5),add=T)
hist(deltaTrue[trueDiff<0],breaks=histogram$breaks,col=rgb(0,0,1,alpha=0.5),add=T)
hist(deltaTrue[trueDiff==0],breaks=histogram$breaks,col=rgb(0,1,0,alpha=0.25),add=T)
points(x,fit[,1],col=2,type="l",lwd=2)
points(x,fit[,5],col=3,type="l",lwd=2)
points(x,fit[,6],col=5,type="l",lwd=2)
points(x,fit[,7],col=4,type="l",lwd=2)
legend(x=min(midPoints),y=1.2*max(c(deltaHist1$counts,deltaHist2$counts)),col=c(rgb(0,0,1,alpha=0.5),rgb(1,0,0,alpha=0.5),rgb(0,1,0,alpha=0.25)),legend=c("True effect < 0","True effect > 0","True effect = 0"),pch=15,cex=0.6)
dev.off()
		}
#For each positive beta, find the probabilities that its deltas are negative, by looking them up on the lower distribution with Bayes' theorem. Then, assuming their effect direction was flipped but the size was the same, look their delta up on the upper distribution.
if(diffSampleSize>nSNP){
	diffSampleSize=nSNP
	message="Warning: diffSampleSize is greater than the number of SNPs. Setting to number of SNPs by default"
	print(message)
	outputLog=c(outputLog,message)
	}
print(paste0("Sampling ",diffSampleSize," effect size differences for each variant"))
findDelta<-function(x,betaRand2,sRand2){ #Assumes dim(x)>0
	if(dim(x)[1]>1){
		index=rep(1:nrow(x),rep(length(betaRand2),nrow(x)))
		delta=(x[index,1]-rep(betaRand2,nrow(x)))/sqrt(x[index,2]^2+rep(sRand2^2,nrow(x)))
		}else{
			delta=(x[1]-betaRand2)/sqrt(x[2]^2+sRand2^2)
			}
	bool=!between(delta,deltaRange[1],deltaRange[2]) #Squash outliers
	delta[bool]=deltaRange[1*(sign(delta[bool])>0)+1]
	return(delta)
	}
bayesLookupL<-function(x,y,beta,betaRand,s,sRand,scaleMeans,scaleSDs){
	set.seed(seed+x[1])
	print(paste0("Computing lower posteriors for SNPs ",x[1],"-",x[length(x)]))
	p1=y[1];p2=y[2];b1=y[3:(degree+2)];b2=y[(degree+3):(2*degree+2)];c1=y[2*degree+3];c2=y[2*degree+4];c1a=y[2*degree+5];c2a=y[2*degree+6] 
	#x is the 'block' parameter
	rand=sample(1:length(betaRand),replace=FALSE,size=diffSampleSize)
	betaRand2=betaRand[rand];sRand2=sRand[rand]
	delta<-findDelta(x=cbind(beta[x],s[x]),betaRand2=betaRand2,sRand2=sRand2)
	lambdaDeriv2=sapply(X=1:degree,y=delta,FUN=make_poly)
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
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	bayes=lambdaL*p1/lambdaMix
	bayesStore=colMeans(matrix(bayes,ncol=length(x),byrow=FALSE)) 
	x=x[bayesStore<highAccThresh]
	if(length(x)>0){
	rand=sample(1:length(betaRand),replace=FALSE,size=diffSampleSizeHighAcc)
	betaRand2=betaRand[rand];sRand2=sRand[rand]
	delta<-findDelta(x=cbind(beta[x],s[x]),betaRand2=betaRand2,sRand2=sRand2)
	lambdaDeriv2=sapply(X=1:degree,y=delta,FUN=make_poly)
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
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	bayes=lambdaL*p1/lambdaMix
	bayesStore[bayesStore<highAccThresh]=(colSums(matrix(bayes,ncol=length(x),byrow=FALSE))+bayesStore[bayesStore<highAccThresh]*diffSampleSize)/(diffSampleSize+diffSampleSizeHighAcc) 
		}
	return(bayesStore)
		}
bayesLookupU<-function(x,y,beta,betaRand,s,sRand,scaleMeans,scaleSDs){
	set.seed(seed+x[1])
	print(paste0("Computing upper posteriors for SNPs ",x[1],"-",x[length(x)]))
	p1=y[1];p2=y[2];b1=y[3:(degree+2)];b2=y[(degree+3):(2*degree+2)];c1=y[2*degree+3];c2=y[2*degree+4];c1a=y[2*degree+5];c2a=y[2*degree+6] 
	#x is the 'block' parameter
	rand=sample(1:length(betaRand),replace=FALSE,size=diffSampleSize)
	betaRand2=betaRand[rand];sRand2=sRand[rand]
	delta<-findDelta(x=cbind(beta[x],s[x]),betaRand2=betaRand2,sRand2=sRand2)
	lambdaDeriv2=sapply(X=1:degree,y=delta,FUN=make_poly)
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
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	bayes=(lambdaU*p2+null*(1-p1-p2))/lambdaMix #Null is included in the upper lookup to allow for BDR=1
	bayesStore=colMeans(matrix(bayes,ncol=length(x),byrow=FALSE)) 
	x=x[bayesStore<highAccThresh]
	if(length(x)>0){
	rand=sample(1:length(betaRand),replace=FALSE,size=diffSampleSizeHighAcc)
	betaRand2=betaRand[rand];sRand2=sRand[rand]
	delta<-findDelta(x=cbind(beta[x],s[x]),betaRand2=betaRand2,sRand2=sRand2)
	lambdaDeriv2=sapply(X=1:degree,y=delta,FUN=make_poly)
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
	lambdaMix=lambdaL*p1+lambdaU*p2+null*(1-p1-p2)
	bayes=(lambdaU*p2+null*(1-p1-p2))/lambdaMix #Null is included in the upper lookup to allow for BDR=1
	bayesStore[bayesStore<highAccThresh]=(colSums(matrix(bayes,ncol=length(x),byrow=FALSE))+bayesStore[bayesStore<highAccThresh]*diffSampleSize)/(diffSampleSize+diffSampleSizeHighAcc) 
		}
	return(bayesStore)
		}
#Return to using raw data for calculating the observed deltas
beta=fullData[,'beta'];s=fullData[,'s'];nSNP=nrow(fullData);fdr=fullData[,'fdr'];rm(fullData)
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
print(paste0("Computing posteriors for ",nSNP*diffSampleSize," effect size differences in blocks of ",blockSize,", using ",threads," ",threadPlural,". Can take more than about ",ceiling(10/threads)," minutes"))
index=rep((1:nSNP),rmultinom(1,size=nSim,prob=(1-fdr)/((1-pNull)*nSNP))[,1]) #Sample of observed SNPs as if they were from alt distribution, as before when choosing alt SNPs for plotting
BDRpos=unlist(mclapply(FUN=bayesLookupL,X=snpBlocks,y=par,beta=abs(beta),betaRand=beta[index],s=s,sRand=s[index],scaleMeans=scaleMeans,scaleSDs=scaleSDs,mc.cores=threads,mc.preschedule=FALSE))
BDRneg=unlist(mclapply(FUN=bayesLookupU,X=snpBlocks,y=par,beta=-abs(beta),betaRand=beta[index],s=s,sRand=s[index],scaleMeans=scaleMeans,scaleSDs=scaleSDs,mc.cores=threads,mc.preschedule=FALSE))
BDRs=BDRpos+BDRneg
BFDR=fdr+BDRs*(1-fdr)
print("Finding tail-area rates")
cuMean=cumsum(sort(BDRs,decreasing=FALSE))/seq_along(BDRs)
tailStats=cuMean[rank(BDRs,ties.method="max")]
cuMean=cumsum(sort(fdr,decreasing=FALSE))/seq_along(fdr)
tailStats=cbind(cuMean[rank(fdr,ties.method="max")],tailStats)
tailBFDRs=tailStats[,1]+tailStats[,2]*(1-tailStats[,1])
print("BFDR.priorsplitteR analysis complete.")
outputLog=c(outputLog,"BFDR.priorsplitteR analysis complete.")
write.table(outputLog,file=paste0(outFileStem,"log.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
return(list(results=data.frame(localBFDR=BFDR,localFDR=fdr,localBDR=BDRs,tailBFDR=tailBFDRs,tailFDR=tailStats[,1],tailBDR=tailStats[,2]),altNullFit=parAltNull,polyFit=par,outputLog=outputLog))
			}#End of BFDR function

