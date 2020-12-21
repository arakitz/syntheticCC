######### Synthetic S1 Chart 
######### Evaluation of the zero-state IC ARL,
######### Define a function of constant k for the control limits of the Xbar chart and of the limit H for the CRL chart
ARLS1fun<-function(k2,H){
M<-H+1 # TPM dimension
p0A<-1-pnorm(k2)+pnorm(-k2) # probability for a point to occur beyond the (LCL,UCL) in Xbar chart
p0B<-1-p0A
################ Construction of the Essential TPM Q in the IC state ##########
Q0<-matrix(0,nrow=M,ncol=M)
Q0[1,1]<-p0B
Q0[M,1]<-p0B
Q0[1,2]<-p0A
for(j in 2:H){
Q0[j,j+1]<-p0B
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
k02<-1.2 # starting value for constant k
arl0<-370.4 # desired IC ARL value
zsarl<-ARLS1fun(k02,H0)
while(zsarl<=arl0){
k02<-k02+0.0001
zsarl<-ARLS1fun(k02,H0)
}
cat("k2:",k02," H:",H0," zsARL0:",zsarl,"\n") 
############ Example
mu0<-60 # IC mean
sigma0<-2 # IC sigma
n<-5 # sample size
LCL1<-mu0-k02*sigma0/sqrt(n)
UCL1<-mu0+k02*sigma0/sqrt(n)
cat(" The LCL for the Xbar chart is",LCL1," the UCL is",UCL1," and the limit H of the CRL chart is",H0,"\n")