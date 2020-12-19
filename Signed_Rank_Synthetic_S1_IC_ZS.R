######### NP Synthetic S1 zero-state IC ARL, H=3, Signed-Rank
ARLSRS1fun<-function(n,UCL,H){
  
    M<-H+1
    LCL<-(-UCL)
    p0A<-1-psignrank(UCL-1,n)+psignrank(LCL,n)
    p0B<-1-p0A
    #################################
    Q0<-matrix(0,nrow=M,ncol=M)
    Q0[1,1]<-p0B
    Q0[M,1]<-p0B
    Q0[1,2]<-p0A
    for(j in 2:H){
      Q0[j,j+1]<-p0B
    }
    ##################################
    qvec<-rep(0,M);qvec[2]<-1
    IDm<-diag(M)
    l1<-rep(1,M)
    arl<-as.vector(qvec%*%solve(IDm-Q0)%*%l1)
    arl
}
########### design
n0<-5
UCLmax<-n0*(n0+1)/2
H0<-3
UCL0<-0
arl0<-370.4

for(UCL0 in 0:UCLmax){
  zsarl<-ARLSRS1fun(n0,UCL0,H0)
  if(zsarl>arl0){break}
  else{
    #   print(c(n0,H0,UCL0,-UCL0,zsarl))
    cat("n:",n0," H:",H0," UCL:",UCL0," LCL:",-UCL0," zsARL0:",zsarl,"\n")
  }
}
################
cat("n:",n0," H:",H0," UCL:",UCL0," LCL:",-UCL0," zsARL0:",zsarl,"\n")
cat("n:",n0," H:",H0," UCLb:",UCL0-1," LCLb:",-UCL0+1," zsARL0:",ARLSRS1fun(n0,UCL0-1,H0),"\n")
###############







