######### Synthetic S2 Chart 
######### Evaluation of the zero-state IC ARL,
######### Define a function of constants k1, k2 for the control and warning limits of the Xbar chart and of the limit H for the CRL chart
ARLS2fun<-function(k1,k2,H){
M<-H+1 # TPM dimension
p0E<-pnorm(k1)-pnorm(k2)+pnorm(-k2)-pnorm(-k1) # probability for a point to occur in the (LCL,LWL] or the [UWL,UCL) in Xbar chart
p0F<-pnorm(k2)-pnorm(-k2) # probability for a point to occur in the (LWL,UWL] in Xbar chart
p0D<-1-p0E-p0F
####### Construction of the Essential TPM Q in the IC state ###################
Q0<-matrix(0,nrow=M,ncol=M)
Q0[1,1]<-p0F
Q0[M,1]<-p0F
Q0[1,2]<-p0E
for(j in 2:H){
Q0[j,j+1]<-p0F
}
##################################
qvec<-rep(0,M);qvec[2]<-1 # initial probabilities vector
IDm<-diag(M) # identity matrix
l1<-rep(1,M) # vector with 1s
arl<-as.vector(qvec%*%solve(IDm-Q0)%*%l1)
arl # zero-state ARL
}
########### design procedure
H0<-3 # specify H
k01<-3.5 # specify k1
k02<-1.2 # starting value for constant k2
arl0<-370.4 # desired IC ARL value
zsarl<-ARLS2fun(k01,k02,H0)
while(zsarl<=arl0){
k02<-k02+0.0001
zsarl<-ARLS2fun(k01,k02,H0)
}
cat("k1:",k01," k2:",k02," H:",H0," zsARL0:",zsarl,"\n")
############ Example
mu0<-60 # IC mean
sigma0<-2 # IC sigma
n<-5 # sample size
LWL<-mu0-k02*sigma0/sqrt(n)
UWL<-mu0+k02*sigma0/sqrt(n)
LCL2<-mu0-k01*sigma0/sqrt(n)
UCL2<-mu0+k01*sigma0/sqrt(n)
cat(" The LCL for the Xbar chart is",LCL2," the LWL is",LWL," the UWL is",UWL," the UCL is",UCL2," and the limit H of the CRL chart is",H0,"\n")