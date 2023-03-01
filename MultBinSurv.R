
#####################################################
# R CODE FOR MULTILEVEL JOINT FRAILTY MODEL FOR     # 
# HIERARCHICALLY CLUSTERED SURVIVAL AND BINARY DATA #
#                                                   #
# BY RICHARD TAWIAH                                 #
# DATE: 1ST MARCH 2023                              #
#####################################################

options(max.print=1000000)
setwd("C:/Users/...")
#Read data into R

#Set inital values of variance components, tolerance level and maximum iteration
thetau0=0.5; thetaU0=0.5; rho0=0.5; sigv0=0.5; sigV0=0.5 
eps.reg=1e-2; eps.var=1e-2; itmax=300

jointBinSurv<-function(dat,hospid,patid,thetau0,thetaU0,sigv0,sigV0,rho0,eps.reg,eps.var,itmax){
  
  p<-ncol(dat)-5
  n<-nrow(dat)
  m<-length(unique(hospid))
  
  R1<-matrix(0,ncol=n,nrow=n)
  for(j in 1:n){
    R1[,j]<-ifelse(patid==j,1,0)
  }

  R2<-matrix(0,ncol=m,nrow=n)
  for(i in 1:m){
    R2[,i]<-ifelse(hospid==i,1,0)
  }

  
  XZ<-data.frame(cbind(dat,R1,R2))
  XZ<-as.matrix(XZ)
  
  XZ.r<-XZ[sort.list(XZ[,3]),] #reorder by surv time   
  time<-as.vector(XZ.r[,3])    # reordered surv time
  indi<-as.vector(XZ.r[,4])    #reordered censoring indicator
  y <- as.vector(XZ.r[,5])     #reordered binary response
  
  rand_vec <- as.matrix(XZ.r[,-(1:(5+p))])
  R1 <- as.matrix(rand_vec[,(1:n)])
  R2 <- as.matrix(rand_vec[,-(1:n)])

  Z<-as.matrix(XZ.r[,6:(5+p)])
  W<-as.matrix(cbind(1,Z)) 
  d<-ncol(W)
  
  M1<-diag(0,n)
  for( j in 1:n)
    for(i in 1:j) M1[j,i]<-1
  
  
  #initial values
  alpha0<-as.vector(rep(0,d))
  beta0<-as.vector(rep(0,p))
  u<-as.vector(rep(0,n))
  U<-as.vector(rep(0,n))
  v<-as.vector(rep(0,m))
  V<-as.vector(rep(0,m))
  
  par0<-as.vector(c(alpha0,beta0,u,U,v,V))
  
  flag.reg<-0
  flag.var<-0
  
  K1<-rbind(cbind(diag(n),diag(0,n)),cbind(diag(0,n),diag(0,n)))
  K2<-rbind(cbind(diag(0,n),diag(n)),cbind(diag(n),diag(0,n)))
  K3<-rbind(cbind(diag(0,n),diag(0,n)),cbind(diag(0,n),diag(n)))
  
  for(iter in 1:itmax){ 
    
    #cat("iter=",iter,"\n")
    
    UG<-diag(0,(d+p+n+n+m+m))
    UG[(d+p+1):(d+p+n),(d+p+1):(d+p+n)]<- (1/((thetau0)*(1-rho0^2)))*diag(n)
    UG[(d+p+n+1):(d+p+n+n),(d+p+1):(d+p+n)]<--(rho0*(1/(sqrt(thetau0*thetaU0)*(1-rho0^2))))*diag(n)
    UG[(d+p+1):(d+p+n),(d+p+n+1):(d+p+n+n)]<--(rho0*(1/(sqrt(thetau0*thetaU0)*(1-rho0^2))))*diag(n)
    UG[(d+p+n+1):(d+p+n+n),(d+p+n+1):(d+p+n+n)]<- (1/((thetaU0)*(1-rho0^2)))*diag(n)
    UG[(d+p+n+n+1):(d+p+n+n+m),(d+p+n+n+1):(d+p+n+n+m)]<-(1/sigv0)*diag(m)
    UG[(d+p+n+n+m+1):(d+p+n+n+m+m),(d+p+n+n+m+1):(d+p+n+n+m+m)]<-(1/sigV0)*diag(m)

    for(BLUP in 1:itmax){
      
      xi<-as.vector(W%*%alpha0+R1%*%u+R2%*%v)
      eta<-as.vector(Z%*%beta0+R1%*%U+R2%*%V)
      
      #1st & 2nd order derivative of l w.r.t xi 
      f1.xi<-as.vector(y-exp(xi)/(1+exp(xi)))
      f2.xi<-diag(as.vector((exp(xi)/(1+exp(xi))^2)))  
      
      w<-diag(as.vector(exp(eta)))  
      A<-diag(as.vector(indi/(t(M1)%*%(exp(eta)))))
      B<-diag(as.vector(M1%*%A%*%rep(1,n)))
      
      #1st & 2nd order derivative of l w.r.t eta
      f1.eta<-as.vector(indi-w%*%M1%*%A%*%rep(1,n))
      f2.eta<-w%*%B-w%*%M1%*%A%*%A%*%t(M1)%*%w
      
      #Score functions 
      dl.dalpha<-t(W)%*%f1.xi
      dl.dbeta<-t(Z)%*%f1.eta
      dl.du<-t(R1)%*%f1.xi-(u*thetaU0-U*rho0*sqrt(thetau0*thetaU0))/(thetau0*thetaU0*(1-rho0^(2)))
      dl.dU<-t(R1)%*%f1.eta-(U*thetau0-u*rho0*sqrt(thetau0*thetaU0))/(thetau0*thetaU0*(1-rho0^(2)))
      dl.dv<-t(R2)%*%f1.xi-(1/sigv0)*v
      dl.dV<-t(R2)%*%f1.eta-(1/sigV0)*V
      
      O.np <- matrix(rep(0,n*p), nrow=n,byrow=TRUE)
      O.nm <- matrix(rep(0,n*m), nrow=n,byrow=TRUE)
      O.nd <- matrix(rep(0,n*d), nrow=n,byrow=TRUE)
      O.nn <- matrix(rep(0,n*n), nrow=n,byrow=TRUE)
      
      XX.bin <- cbind(W, O.np, R1, O.nn, R2, O.nm)
      XX.surv<- cbind(O.nd, Z, O.nn, R1, O.nm, R2)
      XX <- rbind(XX.bin,XX.surv)
      
      zeros <- matrix(rep(0,n*n), nrow=n,byrow=TRUE)
      one <- as.matrix(cbind(f2.xi, zeros))
      two <- as.matrix(cbind(zeros, f2.eta))
      
      f2.xi.eta <- as.matrix(rbind(one,two))
      
      G <- t(XX)%*%f2.xi.eta%*%XX+UG
      G_inv <- chol2inv(chol(G)) # Cholesky decomposition 
      
      Score <- as.vector(c(dl.dalpha,
                           dl.dbeta,
                           dl.du,
                           dl.dU,
                           dl.dv,
                           dl.dV))
      
      par <- par0+G_inv%*%Score     
      
      
      if(max(abs(c((par-par0))))<eps.reg){flag.reg<-1;break}
      par0 <- par
      
      alpha0 <- par[1:d]
      beta0 <- par[(d+1):(d+p)]
      u <- par[(d+p+1):(d+p+n)]
      U <- par[(d+p+n+1):(d+p+n+n)]
      q <- par[(d+p+1):(d+p+n+n)]
      v <- par[(d+p+n+n+1):(d+p+n+n+m)]
      V <- par[(d+p+n+n+m+1):(d+p+n+n+m+m)]

      #cat("alpha0=",alpha0, '\n',
      #    "beta0=",beta0,'\n')
      
    } 
    if(flag.reg==0)stop("Convergence not attained")
    flag.reg<-0
    
    # Variance parameters 
    
    tau <- (G_inv[(d+p+1):(d+p+n+n),(d+p+1):(d+p+n+n)] + q%*%t(q))  
    
    BB1<-sum(diag(K1%*%tau))
    BB2<-sum(diag(K2%*%tau))/2
    BB3<-sum(diag(K3%*%tau))
    
    thetau<-BB1/n; 
    thetaU<-BB3/n; 
    rho<-BB2/sqrt(BB1*BB3)

    sigv<-as.vector(t(v)%*%v+sum(diag(G_inv[(d+p+n+n+1):(d+p+n+n+m),(d+p+n+n+1):(d+p+n+n+m)])))/m
    sigV<-as.vector(t(V)%*%V+sum(diag(G_inv[(d+p+n+n+m+1):(d+p+n+n+m+m),(d+p+n+n+m+1):(d+p+n+n+m+m)])))/m

   # cat("thetau=",thetau,"thetaU=",thetaU,"rho=",rho,
   # "sigv=",sigv,"sigV=",sigV,'\n')
    
    if(max(abs(c((thetau-thetau0),(thetaU-thetaU0),(rho0-rho),
    (sigv-sigv0),(sigV-sigV0))))<eps.var){flag.var<-1;break}
    
    thetau0<-thetau
    thetaU0<-thetaU
    sigv0 <- sigv
    sigV0 <- sigV
    rho0<-rho
    
  } 
 if(flag.var==0)stop("Convergence not attained")
 flag.var<-0

  # Std. error for alpha & beta
  alpha<-par[1:d]
  se.alpha<-sqrt(diag(G_inv)[1:d])
  
  beta<-par[(d+1):(d+p)]
  se.beta<-sqrt(diag(G_inv)[(d+1):(d+p)])
  
  #SE for variance components for patient effects
  D1<-thetau0*diag(n); D2<-rho0*sqrt(thetau0*thetaU0)*diag(n);
  D3<-thetaU0*diag(n)
  
  top_row<-cbind(D1,D2)
  bot_row<-cbind(D2,D3)
  Sigma<-as.matrix(rbind(top_row,bot_row))
  
  dif.is.th1<-(1/(2*thetau0^(2)*thetaU0*(1-rho0^2)))*(rbind(cbind(-2*thetaU0*diag(n),
                                                                    rho0*sqrt(thetau0*thetaU0)*diag(n)),
                                                              cbind(rho0*sqrt(thetau0*thetaU0)*diag(n),
                                                                    diag(0,n))))
    
    dif.is.th2<-(1/(2*thetau0*thetaU0^(2)*(1-rho0^2)))*(rbind(cbind(diag(0,n),
                                                                    rho0*sqrt(thetau0*thetaU0)*diag(n)),
                                                              cbind(rho0*sqrt(thetau0*thetaU0)*diag(n),
                                                                    -2*thetau0*diag(n))))
    
    dif.is.rho<-(1/(thetau0*thetaU0*(1-rho0^2)^2))*(rbind(cbind(2*rho0*thetaU0*diag(n),
                                                                -(1+rho0^(2))*sqrt(thetau0*thetaU0)*diag(n)),
                                                          cbind(-(1+rho0^(2))*sqrt(thetau0*thetaU0)*diag(n),
                                                                2*rho0*thetau0*diag(n))))

  
  J1<-G_inv[(d+p+1):(d+p+n+n),(d+p+1):(d+p+n+n)]%*%dif.is.th1
  J2<-Sigma%*%dif.is.th1
  J3<-G_inv[(d+p+1):(d+p+n+n),(d+p+1):(d+p+n+n)]%*%dif.is.th2
  J4<-Sigma%*%dif.is.th2
  J5<-G_inv[(d+p+1):(d+p+n+n),(d+p+1):(d+p+n+n)]%*%dif.is.rho
  J6<-Sigma%*%dif.is.rho
  
  a11<-sum(diag((J1-J2)%*%(J1-J2)))
  a12<-sum(diag(J1%*%J3+J2%*%J4-2*J1%*%J4))
  a13<-sum(diag(J1%*%J5+J2%*%J6-2*J1%*%J6))
  a22<-sum(diag((J3-J4)%*%(J3-J4)))
  a23<-sum(diag(J3%*%J5+J4%*%J6-2*J3%*%J6))
  a33<-sum(diag((J5-J6)%*%(J5-J6)))
  
  varmat_bin<-2*solve(matrix(c(a11,a12,a13,a12,a22,a23,a13,a23,a33),ncol=3))

 stdvar1<-cbind(c(thetau=thetau,thetaU=thetaU,rho=rho),
                sqrt(diag(varmat_bin)),
                2*(1-pnorm(abs(c(thetau=thetau,thetaU=thetaU,rho=rho)/sqrt(diag(varmat_bin))))))
  
  dimnames(stdvar1)<-list(c("thetau","thetaU","rho"),c("Estimate","SE","p-value"))
  stdvar1<-round(stdvar1,3)

  #SE for variance components for hospital effects

  b11<-sum(diag((diag(m)-G_inv[(d+p+n+n+1):(d+p+n+n+m),(d+p+n+n+1):(d+p+n+n+m)]/sigv0)%*%(diag(m)-G_inv[(d+p+n+n+1):(d+p+n+n+m),(d+p+n+n+1):(d+p+n+n+m)]/sigv0)))/sigv0^2
  b12<-sum(diag(G_inv[(d+p+n+n+1):(d+p+n+n+m),(d+p+n+n+m+1):(d+p+n+n+m+m)]%*%G_inv[(d+p+n+n+m+1):(d+p+n+n+m+m),(d+p+n+n+1):(d+p+n+n+m)]))/(sigv0*sigV0)^2
  b22<-sum(diag((diag(m)-G_inv[(d+p+n+n+m+1):(d+p+n+n+m+m),(d+p+n+n+m+1):(d+p+n+n+m+m)]/sigV0)%*%(diag(m)-G_inv[(d+p+n+n+m+1):(d+p+n+n+m+m),(d+p+n+n+m+1):(d+p+n+n+m+m)]/sigV0)))/sigV0^2

  varmat_surv<-2*solve(matrix(c(b11,b12,b12,b22),nrow=2,byrow=TRUE)) 
  se.var.par2<-sqrt(diag(varmat_surv))
  
  stdvar2<-cbind(c(sigv=sigv,sigV=sigV),se.var.par2,2*(1-pnorm(abs(c(sigv=sigv,sigV=sigV)/se.var.par2)))) 
  dimnames(stdvar2)<-list(c("sigmav","sigmaV"),c("Estimate","SE","p-value"))
    
  stdvar2<-round(stdvar2,3)

  #DISPLAY RESULTS
  
  #Binary submodel
  OR <- exp(alpha)
  L.CI_OR <- exp(alpha-1.96*se.alpha)
  U.CI_OR <- exp(alpha+1.96*se.alpha)
  
  ealpha<-cbind(alpha,se.alpha,OR,L.CI_OR,U.CI_OR,
                2*(1-pnorm(abs(alpha/se.alpha))))
  
  names(alpha)<-c("Intercept",colnames(dat[,-(1:5)]))
  dimnames(ealpha)<-list(names(alpha),c("Estimate","SE","OR","L.CI","U.CI","p-value"))
  ealpha<-round(ealpha,3)
  options(digits=3)   
  
  cat('\n',"----------------------------------",'\n')
  
  #Survival submodel
  HR <- exp(beta)
  L.CI_HR <- exp(beta-1.96*se.beta)
  U.CI_HR <- exp(beta+1.96*se.beta)
  
  ebeta<-cbind(beta,se.beta,HR,L.CI_HR,U.CI_HR,
               2*(1-pnorm(abs(beta/se.beta))))
  names(beta)<-c(colnames(dat[,-(1:5)]))
  dimnames(ebeta)<-list(names(beta),c("Estimate","SE","HR","L.CI","U.CI","p-value"))
  ebeta<-round(ebeta,3)
  options(digits=3) 

  cat('\n',"----------------------------------",'\n')
  
  cat('\n',"Multilevel joint frailty model for clustered survival and binary data:\n\n")
  
  cat('\n',"Binary submodel",'\n')
  print(ealpha)

  cat('\n')
  
  cat('\n',"Survival submodel",'\n')
  print(ebeta)
  
  cat('\n')
  
  cat('\n',"Variance parameters (Patient frailty)",'\n')
  print(stdvar1)

  cat('\n')
  
  cat('\n',"Variance parameters (Hospital frailty)",'\n')
  print(stdvar2)
  
  cat('\n',"----------------------------------",'\n')
  return("end of program")
   
}
jointBinSurv(dat,hospid=dat[,1],patid=dat[,2],thetau0,thetaU0,sigv0,sigV0,rho0,eps.reg,eps.var,itmax)




