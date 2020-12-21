######### Synthetic S4 Chart 
######### Evaluation of the zero-state IC ARL
H<-3 # specify the limit of the CRL chart
M<-(H+1)^2 # TPM dimension
######### Define a function of constants k1 and k2 for the control and the warning limits of the Xbar chart
######### The limit H for the CRL chart equals 3
ARLS4fun<-function(k1,k2){
p0E<-pnorm(k1)-pnorm(k2) # probability for a point in the Xbar to be in (LWL,UCL2]
p0F<-pnorm(k2)-pnorm(-k2) # probability for a point in the Xbar to be in (LWL,UWL]
p0G<-pnorm(-k2)-pnorm(-k1) # probability for a point in the Xbar to be in (LCL2,UWL]
p0D<-1-p0E-p0F-p0G
########### Construction of the Essential TPM Q in the IC state, case H=3 #####
Q0<-matrix(0,nrow=M,ncol=M)
Q0[1,1]<-p0F;Q0[2,3]<-p0F;Q0[3,5]<-p0F;Q0[4,7]<-p0F;
Q0[5,1]<-p0F;Q0[6,15]<-p0F;Q0[7,13]<-p0F;Q0[8,9]<-p0F;
Q0[9,10]<-p0F;Q0[10,1]<-p0F;Q0[11,5]<-p0F;Q0[12,3]<-p0F;
Q0[13,1]<-p0F;Q0[14,11]<-p0F;Q0[15,13]<-p0F;Q0[16,15]<-p0F;
##################################
Q0[1,2]<-p0E;Q0[13,2]<-p0E;Q0[15,12]<-p0E;Q0[16,14]<-p0E;
################################
Q0[1,16]<-p0G;Q0[2,4]<-p0G;Q0[3,6]<-p0G;Q0[5,16]<-p0G
##################################
qvec<-rep(0,M);qvec[8]<-1 # initial probabilities vector
IDm<-diag(M) # identity matrix
l1<-rep(1,M) # vector with 1s
arl<-as.vector(qvec%*%solve(IDm-Q0)%*%l1)
arl # zero-state ARL
}
########### design procedure
k01<-3.5 # specify the value for constant k1
k02<-1.2 # starting value for constant k2
arl0<-370.4 # desired IC ARL value
zsarl<-ARLS4fun(k01,k02)
while(zsarl<=arl0){
k02<-k02+0.0001
zsarl<-ARLS4fun(k01,k02)
}
cat("k1:",k01," k2:",k02," H:",H," zsARL0:",zsarl,"\n")
############ Example
mu0<-60 # IC mean
sigma0<-2 # IC sigma
n<-5 # sample size
LCL2<-mu0-k01*sigma0/sqrt(n)
UCL2<-mu0+k01*sigma0/sqrt(n)
LWL<-mu0-k02*sigma0/sqrt(n)
UWL<-mu0+k02*sigma0/sqrt(n)
cat(" The LCL for the Xbar chart is",LCL2," the LWL is",LWL," the UWL is",UWL," the UCL is",UCL2," and the limit H of the CRL chart is",H,"\n")