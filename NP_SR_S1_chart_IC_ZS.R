######### Nonparametric Synthetic S1 Chart with Signed-Rank Chart as sub-chart 
######### Evaluation of the zero-state IC ARL,
######### Define a function of the control limit UCL of the signed-rank chart and of the limit H for the CRL chart
ARLSRS1fun<-function(n,UCL,H){
    M<-H+1 # TPM dimension
    LCL<-(-UCL) # symmetric limits for the two-sided signed-rank sub-chart
    p0A<-1-psignrank(UCL-1,n)+psignrank(LCL,n) # probability for a point in signed-rank chart to be beyond the (LCL,UCL)
    p0B<-1-p0A
###### Construction of the Essential TPM Q in the IC state ###########################
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
########### design
n0<-5 # specify sample size n
UCLmax<-n0*(n0+1)/2 # Maximum possible value for UCL
H0<-3 # specify the control limit for the CRL chart
arl0<-370.4 # desired IC ARL value
for(UCL0 in 0:UCLmax){
  zsarl<-ARLSRS1fun(n0,UCL0,H0)
  if(zsarl>arl0){break}
  else{
    cat("n:",n0," H:",H0," UCL:",UCL0," LCL:",-UCL0," zsARL0:",zsarl,"\n")
  }
}
################
### The results are printed on screen, providing the attained IC zsARL for (LCL,UCL)
### that exceeds the desired IC ARL value, as well as for (LCL+1,UCL-1) which results
### in an IC zsARL value that it is below the desired IC ARL.
cat("n:",n0," H:",H0," UCL:",UCL0," LCL:",-UCL0," zsARL0:",zsarl,"\n")
cat("n:",n0," H:",H0," UCLb:",UCL0-1," LCLb:",-UCL0+1," zsARL0:",ARLSRS1fun(n0,UCL0-1,H0),"\n")
###############







