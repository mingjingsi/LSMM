##### Functions to generate data #####
generate_data <- function(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2){
  # fixed effects
  Z         <- rep(0, M*L)
  indexZ    <- sample(M*L, M*L*Z.perc)
  Z[indexZ] <- 1
  Z         <- matrix(Z, M, L)

  # random effects
  A         <- rep(0, M*K)
  indexA    <- sample(M*K, M*K*A.perc)
  A[indexA] <- 1
  A         <- matrix(A, M, K)

  # eta (latent variable which indicate whether the annotation is relevant to the phenotype)
  eta           <- rep(0, K)
  indexeta      <- sample(K, K*omega)
  eta[indexeta] <- 1

  # beta (coefficients of random effects)
  beta           <- rep(0, K)
  beta[indexeta] <- rnorm(K*omega, 0, sqrt(sigma2))

  # gamma (latent variable which indicate whether the SNP is associated with the phenotype)
  pi1               <- sigma(beta0 + Z %*% b + A %*% beta)
  gamma             <- rep(0, M)
  indexgamma        <- (runif(M) < pi1)
  gamma[indexgamma] <- 1

  # Pvalue (p-values of the phenotype)
  Pvalue             <- runif(M)
  Pvalue[indexgamma] <- rbeta(sum(indexgamma), alpha, 1)

  return( list(Z = Z, A = A, Pvalue = Pvalue, beta = beta, pi1 = pi1, eta = eta,
               gamma = gamma))
}

generate_data_corr <- function(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, q){
  # fixed effects
  Z         <- rep(0, M*L)
  indexZ    <- sample(M*L, M*L*Z.perc)
  Z[indexZ] <- 1
  Z         <- matrix(Z, M, L)

  # random effects
  corr <- matrix(0, K, K)
  corr[1:10, 1:10] <- matrix(q, 10, 10)
  diag(corr) <- 1
  preA <- mvrnorm(M, rep(0, K), corr)
  A <- t(t(preA) < apply(preA, 2, quantile, probs = A.perc))

  # eta (latent variable which indicate whether the annotation is relevant to the phenotype)
  eta           <- rep(0, K)
  indexeta      <- sample(K, K*omega)
  eta[indexeta] <- 1

  # beta (coefficients of random effects)
  beta           <- rep(0, K)
  beta[indexeta] <- rnorm(K*omega, 0, sqrt(sigma2))

  # gamma (latent variable which indicate whether the SNP is associated with the phenotype)
  pi1               <- sigma(beta0 + Z %*% b + A %*% beta)
  gamma             <- rep(0, M)
  indexgamma        <- (runif(M) < pi1)
  gamma[indexgamma] <- 1

  # Pvalue (p-values of the phenotype)
  Pvalue             <- runif(M)
  Pvalue[indexgamma] <- rbeta(sum(indexgamma), alpha, 1)

  return( list(Z = Z, A = A, Pvalue = Pvalue, beta = beta, pi1 = pi1, eta = eta,
               gamma = gamma))
}

generate_data_corr_ZA <- function(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, q){
  # fixed effects & random effects
  corr <- matrix(0, K+L, K+L)
  corr[1:(L+50), 1:(L+50)] <- q
  corr[1:(L+50), 1:(L+50)] <- q
  diag(corr) <- 1
  preZA <- mvrnorm(M, rep(0, K+L), corr)
  ZA <- t(t(preZA) < apply(preZA, 2, quantile, probs = A.perc))
  Z <- ZA[, 1:L]
  A <- ZA[, (L+1):(K+L)]

  # eta (latent variable which indicate whether the annotation is relevant to the phenotype)
  eta           <- rep(0, K)
  indexeta      <- sample(K, K*omega)
  eta[indexeta] <- 1

  # beta (coefficients of random effects)
  beta           <- rep(0, K)
  beta[indexeta] <- rnorm(K*omega, 0, sqrt(sigma2))

  # gamma (latent variable which indicate whether the SNP is associated with the phenotype)
  pi1               <- sigma(beta0 + Z %*% b + A %*% beta)
  gamma             <- rep(0, M)
  indexgamma        <- (runif(M) < pi1)
  gamma[indexgamma] <- 1

  # Pvalue (p-values of the phenotype)
  Pvalue             <- runif(M)
  Pvalue[indexgamma] <- rbeta(sum(indexgamma), alpha, 1)

  return( list(Z = Z, A = A, Pvalue = Pvalue, beta = beta, pi1 = pi1, eta = eta,
               gamma = gamma))
}

generate_data_probit <- function(M, L, K, alpha, Z.perc, A.perc, beta0, b, omega, sigma2, r){
  # fixed effects
  Z         <- rep(0, M*L)
  indexZ    <- sample(M*L, M*L*Z.perc)
  Z[indexZ] <- 1
  Z         <- matrix(Z, M, L)

  # random effects
  A         <- rep(0, M*K)
  indexA    <- sample(M*K, M*K*A.perc)
  A[indexA] <- 1
  A         <- matrix(A, M, K)

  # eta (latent variable which indicate whether the annotation is relevant to the phenotype)
  eta           <- rep(0, K)
  indexeta      <- sample(K, K*omega)
  eta[indexeta] <- 1

  # beta (coefficients of random effects)
  beta           <- rep(0, K)
  beta[indexeta] <- rnorm(K*omega, 0, sqrt(sigma2))

  # gamma (latent variable which indicate whether the SNP is associated with the phenotype)
  sigmae2           <- var(Z %*% b + A %*% beta)/r # r is signal-noise ratio
  y                 <- beta0 + Z %*% b + A %*% beta + sqrt(sigmae2) * rnorm(M)
  gamma             <- rep(0, M)
  indexgamma        <- (y > 0)
  gamma[indexgamma] <- 1

  # Pvalue (p-values of the phenotype)
  Pvalue             <- runif(M)
  Pvalue[indexgamma] <- rbeta(sum(indexgamma), alpha, 1)

  return( list(Z = Z, A = A, Pvalue = Pvalue, beta = beta, eta = eta, gamma = gamma))
}

generate_data_probit_dist <- function(M, L, K, Z.perc, A.perc, beta0, b, omega, sigma2, r, dist){
  # fixed effects
  Z         <- rep(0, M*L)
  indexZ    <- sample(M*L, M*L*Z.perc)
  Z[indexZ] <- 1
  Z         <- matrix(Z, M, L)

  # random effects
  A         <- rep(0, M*K)
  indexA    <- sample(M*K, M*K*A.perc)
  A[indexA] <- 1
  A         <- matrix(A, M, K)

  # eta (latent variable which indicate whether the annotation is relevant to the phenotype)
  eta           <- rep(0, K)
  indexeta      <- sample(K, K*omega)
  eta[indexeta] <- 1

  # beta (coefficients of random effects)
  beta           <- rep(0, K)
  beta[indexeta] <- rnorm(K*omega, 0, sqrt(sigma2))

  # gamma (latent variable which indicate whether the SNP is associated with the phenotype)
  sigmae2           <- var(Z %*% b + A %*% beta)/r # r is signal-noise ratio
  y                 <- beta0 + Z %*% b + A %*% beta + sqrt(sigmae2) * rnorm(M)
  gamma             <- rep(0, M)
  indexgamma        <- (y > 0)
  gamma[indexgamma] <- 1

  # Pvalue (p-values of the phenotype)
  Pvalue <- runif(M)
  fz <- match.fun(dist)
  z <- fz(sum(indexgamma))
  Pvalue[indexgamma] <- pnorm(abs(z), lower.tail = FALSE)*2

  return( list(Z = Z, A = A, Pvalue = Pvalue, beta = beta, eta = eta, gamma = gamma))
}

spiky <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 0.25)
  if(sum(r <= 0.2) != 0)
    z[which(r <= 0.2)] <- rnorm(sum(r <= 0.2), 0, 0.5)
  if(sum(r > 0.2 & r <= 0.4) != 0)
    z[which(r > 0.2 & r <= 0.4)] <- rnorm(sum(r > 0.2 & r <= 0.4), 0, 1)
  if(sum(r > 0.8) != 0)
    z[which(r > 0.8)] <- rnorm(sum(r > 0.8), 0, 2)

  return(z)
}

near_normal <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 1)
  if(sum(r <= 1/3) != 0)
    z[which(r <= 1/3)] <- rnorm(sum(r <= 1/3), 0, 2)

  return(z)
}

skew <- function(N){
  r <- runif(N)
  z <- rnorm(N, -2, 2)
  if(sum(r <= 0.25) != 0)
    z[which(r <= 0.25)] <- rnorm(sum(r <= 0.25), -1, 1.5)
  if(sum(r > 0.25 & r <= (0.25+1/3)) != 0)
    z[which(r > 0.25 & r <= (0.25+1/3))] <- rnorm(sum(r > 0.25 & r <= (0.25+1/3)), 0, 1)
  if(sum(r > 5/6) != 0)
    z[which(r > 5/6)] <- rnorm(sum(r > 5/6), 1, 1)

  return(z)
}

big_normal <- function(N){
  r <- runif(N)
  z <- rnorm(N, 0, 4)

  return(z)
}

# sigmoid function
sigma <- function(x){
  y <- 1/(1+exp(-x))
  return (y)
}



##### Functions in cmfdr #####
run_cmlocfdr=function(Pvalue=1,P,X,bases_X=FALSE,K=2,knots=NULL,
                      nIter=160,thin=1,burnIn=10,SSA=1,SSG=1,MA=3,MG=3,theoNULL=FALSE,mu,inits=NULL){

  library(pgnorm)
  library(locfdr)
  library(splines)

  ## Inputs:
  ## P: N vector of P-values or Z scores based on the input from mainfile.R;
  ## X: N x Q matrix of covariates
  ## nIter: Number of MCMC iterations
  ## thin: thinning rate
  ## burnIn: burn-in number
  ## SSA: increase step-size for alpha draw
  ## SS: increase step-size for gamma draw
  ## MA: num of multiple try for alpha
  ## MG: num of multiple try for gamma
  ## mu: origin of gamma distribution, do not change; not used for generalized normal

  if (Pvalue==1){ # input P value
    N=length(P)
    snpid=1:N
    all.complete=complete.cases(cbind(P,X))
    X=X[all.complete,]
    P=P[all.complete]

    P[P==0]=min(P[P>0]);P[P==1]=max(P[P<1])
    Z=-qnorm(P/2);N=length(Z)
    if(!colnames(X)[1]=="Intcpt"){X=cbind(Intcpt=1,X)}
  }else{
    Z=P; #input Z score;
  }

  X_bases=NULL
  if(!is.null(bases_X)){
    for(p in 1:length(bases_X)){
      X_p=X[,bases_X[p]]
      range=c(min(X_p),max(X_p))
      delta=(range[2]-range[1])/200
      grid=seq(range[1],range[2]+delta,by=delta)
      knots=c(min(grid),quantile(X_p,.5),quantile(X_p,.67),max(grid))
      phi.mat=bs(grid,knots=knots[2:(length(knots)-1)],degree=3,intercept=FALSE,
                 Boundary.knots=c(knots[1],knots[length(knots)]))
      for(k in 1:(K+3)){
        phi.mat[,k]=phi.mat[,k]/(sum(phi.mat[,k])*delta)
      }
      #plot(grid,phi.mat[,1],type="l",ylab="density",xlab="annotation score")
      #title(main="density basis functions")
      #for(k in 2:K){
      #	lines(grid,phi.mat[,k],type="l",col=k+1,lwd=.5)
      #}
      Phi.mat=phi.mat
      for(k in 1:(K+3)){
        for(j in 1:length(grid)){Phi.mat[j,k]=sum(phi.mat[grid<=grid[j],k]*delta)}
      }
      #plot(grid,Phi.mat[,2],type="l")
      #lines(grid,Phi.mat[,3],type="l",col=k)
      X_bases_p=array(NA,dim=c(length(X_p),2))
      for(g in 1:length(grid)){
        #print(c(p,g))
        if(length(X_p[abs(X_p-grid[g])<=delta/2])>0){
          tmp=dim(rbind(X_bases_p[abs(X_p-grid[g])<=delta/2,]))[1]
          X_bases_p[abs(X_p-grid[g])<=delta/2,]=cbind(rep(Phi.mat[g,2],tmp),rep(Phi.mat[g,3],tmp))
        }
      }
      X_bases=cbind(X_bases,X_bases_p)
    }
    name=colnames(X)[bases_X[1]]
    colnames(X_bases)=rep(1:2,length(bases_X))
    colnames(X_bases)[1:2]=c(paste(name,"_1",sep=""),paste(name,"_2",sep=""))
    if(length(bases_X)>1){
      for(p in 2:length(bases_X)){
        name=colnames(X)[bases_X[p]]
        colnames(X_bases)[(2*(p-1)+1):(2*(p-1)+2)]=c(paste(name,"_1",sep=""),paste(name,"_2",sep=""))
      }
    }
    colnames(X_bases)
    X=cbind(Intcpt=1,X_bases,X[,-c(1,bases_X)])
  }

  save(file="data_inputs.R",Z,X,N)

  MCMCfit=cmlFDR_GammaDist(Z,X,nIter=nIter,burnIn=burnIn,thin=thin,SSA=SSA,SSG=SSG,MA=MA,MG=MG,mu=mu,
                           theoNULL=theoNULL,inits=inits)

  ALPHA_array=MCMCfit[[1]]
  BETA_array=MCMCfit[[2]]
  GAMMA_array=MCMCfit[[3]]
  SIGMA_SQ_array=MCMCfit[[4]]
  first=1
  Alpha=ALPHA_array[[first]]
  Beta=BETA_array[[first]]
  Gamma=GAMMA_array[[first]]
  if(theoNULL==FALSE){Sigma_sq=SIGMA_SQ_array[[first]]}
  for(j in (first+1):length(ALPHA_array)){
    Alpha=Alpha+ALPHA_array[[j]]
    Beta=Beta+BETA_array[[j]]
    Gamma=Gamma+GAMMA_array[[j]]
    if(theoNULL==FALSE){Sigma_sq=Sigma_sq+SIGMA_SQ_array[[j]]}
  }
  Alpha=cbind(Alpha/length(first:length(ALPHA_array)))
  Beta=Beta/length(first:length(ALPHA_array))
  Gamma=cbind(Gamma/length(first:length(ALPHA_array)))
  if(theoNULL==FALSE){Sigma_sq=Sigma_sq/length(first:length(ALPHA_array))}
  if(theoNULL==TRUE){Sigma_sq=1}

  f0<-2*dnorm(abs(Z),mean=0,sd=Sigma_sq^.5)
  f1<-f0
  for (i in 1:length(Z)){
    X_i=rbind(X[i,])
    f1[i]<-2*dgamma(abs(Z[i])-.68,shape=exp(X_i%*%Alpha),rate=Beta)
  }
  p0=1;
  p1=exp(X%*%Gamma)
  pi0=p0/(p0+p1)
  pi1=1-pi0
  cmfdr=pi0*f0/(pi0*f0+pi1*f1)

  X_mean=rbind(apply(X,2,mean))
  f1=2*dpgnorm(Z,p=Beta,mean=0,sigma=exp(X_mean%*%Alpha)*sqrt(gamma(3/Beta)/gamma(1/Beta)))
  p1=exp(X_mean%*%Gamma)
  pi0=p0/(p0+p1)
  pi1=1-pi0
  fdr=pi0*f0/(pi0*f0+pi1*f1)

  z.locfdr=locfdr(c(Z,-Z),nulltype=1,plot=0)
  efron_fdr=z.locfdr$fdr[1:length(Z)]

  results=list()
  results[[1]]=MCMCfit[[1]]
  results[[2]]=MCMCfit[[2]]
  results[[3]]=MCMCfit[[3]]
  results[[4]]=MCMCfit[[4]]
  results[[5]]=MCMCfit[[5]]
  results[[6]]=array(NA,dim=c(N,3))
  results[[6]][,1]=efron_fdr
  results[[6]][,2]=fdr
  results[[6]][,3]=cmfdr

  return(results)
}

cmlFDR_GammaDist=function (Z,X,nIter=1100,burnIn=100,thin=5,initNULL=0.95,simulate=FALSE,
                           SSA=1,SSG=1,MA=3,MG=3,mu=0.68,theoNULL=FALSE,inits=NULL)
{
  #SSA:scale the diagnal of the covariance matrix in MH of Alpha draw to increase/decrease the step size
  #SSG:scale the diagnal of the covariance matrix in MH of Gamma draw to increase/decrease the step size


  ###########################
  #Load packages/functions
  ###########################

  library(tmvtnorm)
  library(mnormt)
  library(magic)
  library(locfdr)
  library(arm)
  library(pscl)

  if(colnames(X)[1] != "Intcpt") X=cbind(Intcpt=1,X)
  #print(X[1:10,])
  # hyperparameters

  B01=0.001;B02=0.001; #Prior: P(Beta) ~ Gamma(B01,B02)
  a0=0.001;b0=0.001; #Prior: P(sigma^2) ~ IVG(a0,b0) #b0 is scale;
  Sigma_Gamma=matrix(0,nrow=dim(X)[2],ncol=dim(X)[2])
  diag(Sigma_Gamma)=10000; #Prior: P(Gamma) ~ N(0,Sigma_Gamma)

  Sigma_Alpha=matrix(0,nrow=dim(X)[2],ncol=dim(X)[2])
  diag(Sigma_Alpha)=10000; #Prior: P(Alpha) ~ N(0,Sigma_Alpha)
  df=4 # degrees of freedom for multivariate t proposal

  # data
  N=dim(X)[1]
  M=dim(X)[2]

  # Parameter arrays

  ALPHA_array=list()
  BETA_array=list()
  if(theoNULL==FALSE){
    SIGMA_SQ_array=list()
  }
  GAMMA_array=list()
  Accp_Rate_array=list()  	#Alpha draw accept rate
  Accp_Rate_array_g=list()	#Gamma draw accept rate

  array_ind=0

  # Initialize parameters

  if(is.null(inits)){
    Alpha=cbind(rep(0,M));Alpha_mean=Alpha
    pi0=initNULL #min(.98,locfdr(Z,nulltype=1,plot=0)$fp0[5,3]);
    gamma0=log((1-pi0)/pi0) #intercept for Non-NULL gamma
    Gamma=array(0,dim=c(M,1));Gamma_mean=Gamma; Gamma[1,]=gamma0
    Beta=0.1;Beta_mean=Beta
    Phi=1-as.numeric(abs(Z)< sort(abs(Z))[round(pi0*N)]);Phi_mean=0*Phi
    Sigma_sq=1
  }
  if(!is.null(inits)){
    last=length(inits[[1]])
    Alpha=cbind(inits[[1]][[last]]);Alpha_mean=Alpha
    Beta=inits[[2]][[last]];Beta_mean=Beta
    Gamma=cbind(inits[[3]][[last]]);Gamma_mean=Gamma
    if(theoNULL==FALSE){
      Sigma_sq=inits[[4]][[last]];Sigma_sq_mean=Sigma_sq
    }
    if(theoNULL==TRUE){
      Sigma_sq=1
    }
    Phi=inits[[5]];Phi_mean=0*Phi
  }
  PHI_match_rate=NULL


  for(iter in 1:nIter){

    print(iter)

    Z1=abs(Z[!Phi==0])
    X1=X[!Phi==0,]


    ## Draw ALPHA

    if(det(t(X1)%*%X1) != 0){ #avoid singular
      obj=Draw_Alpha_log_M_mu(Alpha,Z1,X1,Beta,Phi,df,MA,Sigma_Alpha,SSA,mu)
      Alpha=obj$par;
    }
    print("Alpha");
    print(Alpha);


    ## Draw BETA

    Beta=rgamma(1,shape=B01+sum(exp(X1%*%Alpha)),rate=B02+sum(Z1-mu))
    print(paste("Beta",Beta));


    ## Draw GAMMA

    #Gamma draw: Multiple try MH
    objg=Draw_Gamma_log_M(Gamma,Z,X,Phi,df,Sigma_Gamma,SSG,MG)

    Gamma=objg$par;
    print("Gamma");
    print(Gamma);


    if(theoNULL==FALSE){
      ## Draw SIGMA_SQ;
      Z0=abs(Z[Phi==0])
      Sigma_sq=rigamma(1,alpha=a0+length(Z0)/2,beta=b0+(Z0%*%Z0)/2);
      print(paste("Sigma_sq",Sigma_sq));
    }


    ## Draw PHI;

    log_P_phi=cbind(X%*%Gamma,0)
    log_P_phi[abs(Z)>mu,1]=log_P_phi[abs(Z)>mu,1]+((exp(X[abs(Z)>mu,]%*%Alpha)-1)*
                                                     log(abs(Z[abs(Z)>mu])-mu)-Beta*(abs(Z[abs(Z)>mu])-mu))+(exp(X[abs(Z)>mu,]%*%Alpha))*
      log(Beta)-lgamma(exp(X[abs(Z)>mu,]%*%Alpha))
    #log_P_phi[,2]=log_P_phi[,2]+log(2)-0.5*log(2*pi)-0.5*Z^2;

    if(theoNULL==TRUE){
      log_P_phi[,2]=log_P_phi[,2]+log(2)-0.5*log(2*pi)-0.5*Z^2;
    }else{
      log_P_phi[,2]=log_P_phi[,2]+log(2)-0.5*log(2*pi*Sigma_sq)-(1/(2*Sigma_sq))*Z^2;
    }

    P_phi=exp(log_P_phi)/apply(exp(log_P_phi),1,sum) #declare the variable P_phi;
    P_phi[abs(Z)>mu,1]=1/(1+exp(log_P_phi[abs(Z)>mu,2]-log_P_phi[abs(Z)>mu,1]))
    P_phi[,2]=1/(1+exp(log_P_phi[,1]-log_P_phi[,2]))

    P_phi[abs(Z)<=mu,1]=0;P_phi[abs(Z)<=mu,2]=1;

    Phi_new=Phi #declare variable Phi_new, create a vector;
    for(i in 1:N){
      Phi_new[i]=sample(c(1,0),size=1,replace=TRUE,prob=P_phi[i,])
    }
    Phi=Phi_new

    if(simulate == TRUE) print(sum(Phi == Phi_true)/N)


    ## Save results after thin

    if(iter%%thin==0 & iter>=burnIn){
      array_ind=array_ind+1
      ALPHA_array[[array_ind]]=Alpha

      GAMMA_array[[array_ind]]=Gamma

      BETA_array[[array_ind]]=Beta

      if(theoNULL==FALSE){
        SIGMA_SQ_array[[array_ind]]=Sigma_sq
      }

      if(det(t(X1)%*%X1) != 0){
        Accp_Rate_array[[array_ind]]=obj$accp
      }
      else Accp_Rate_array[[array_ind]]=0

      Accp_Rate_array_g[[array_ind]]=objg$accp;

      print("Alpha mean:");
      Alpha_mean=((array_ind-1)*Alpha_mean+ALPHA_array[[array_ind]])/array_ind
      if(simulate==TRUE) print(cbind(Alpha_mean,Alpha_true))
      if(simulate==FALSE) print(Alpha_mean)


      print("Beta mean:");
      if(simulate==TRUE) print(cbind(mean(as.numeric(BETA_array)),BETA_true))
      if(simulate==FALSE) print(mean(as.numeric(BETA_array)))

      if(theoNULL==FALSE){
        print("Sigma_sq mean:");
        if(simulate==TRUE) print(cbind(mean(as.numeric(SIGMA_SQ_array)),SIGMA_SQ_true))
        if(simulate==FALSE) print(mean(as.numeric(SIGMA_SQ_array)))

      }

      print("Gamma mean:");
      Gamma_mean=((array_ind-1)*Gamma_mean+GAMMA_array[[array_ind]])/array_ind
      if(simulate==TRUE) print(cbind(Gamma_mean,Gamma_true))
      if(simulate==FALSE) print(Gamma_mean)

      print(paste("Multiple-try MH Accept Rate for Alpha (mean):",mean(as.numeric(Accp_Rate_array))));
      print(paste("Multiple-try MH Accept Rate for Gamma (mean):",mean(as.numeric(Accp_Rate_array_g))));

      if(simulate==TRUE) {
        PHI_match_rate=rbind(PHI_match_rate,sum(Phi==Phi_true)/N);
        print(paste("Phi matching rate (mean):",mean(PHI_match_rate)))
      }

      #probability of each SNP being Non-NULL, average of Phi over the iterations saved;
      Phi_mean=((array_ind-1)*Phi_mean+Phi)/array_ind;


      results=list()
      results[[1]]=ALPHA_array
      results[[2]]=BETA_array
      results[[3]]=GAMMA_array
      if(theoNULL==TRUE){
        results[[4]]=1
      }
      if(theoNULL==FALSE){
        results[[4]]=SIGMA_SQ_array
      }
      results[[5]]=Phi
      save(X,Z,results,file="mcmc_intermediate_outputs.R")

    }



  }

  #Calculate SD;
  #Alpah
  print("Alpha SD:");
  for(i in 1:dim(X)[2]){
    print(sd(as.numeric(unlist(lapply(ALPHA_array, function(x) x[i])))))
  }

  #beta
  print(paste("Beta SD:",sd(as.numeric(BETA_array))))

  if(theoNULL==FALSE){
    #Sigma_sq
    print(paste("Sigma SD:",sd(as.numeric(SIGMA_SQ_array))))
  }

  #Gamma
  print("Gamma SD:");
  for(i in 1:dim(X)[2]){
    print(sd(as.numeric(unlist(lapply(GAMMA_array, function(x) x[i])))))
  }


  #return results;
  #return results;
  results=list()
  results[[1]]=ALPHA_array
  results[[2]]=BETA_array
  results[[3]]=GAMMA_array
  if(theoNULL==TRUE){
    results[[4]]=1
  }
  if(theoNULL==FALSE){
    results[[4]]=SIGMA_SQ_array
  }
  results[[5]]=Phi


  return(results)
}


log_p_alpha=function(Alpha,Z,W,B,SA,MU){
  gammafcn=lgamma(exp(W%*%Alpha));
  #if (gammafcn == Inf || is.na (gammafcn) ) {print(paste("Gamma function",gammafcn));}
  log_p=sum(exp(W%*%Alpha)*log(Z-MU)-gammafcn)+sum(exp(W%*%Alpha))*log(B)-0.5*t(Alpha)%*%solve(SA)%*%Alpha
  #print(log_p)
  return(log_p)
}

## Draw Alpha,log-scale MH alg:multiple-try;
Draw_Alpha_log_M_mu=function(Alpha,Z,X,Beta,Phi,df,Multiple,Sigma_Alpha,SSA,mu)
{

  log_p_Alpha_star=rep(0,Multiple)
  log_p_Alpha_2star=rep(0,Multiple)
  p=rep(0,Multiple)
  den=0;num=0;

  sigma=solve(t(X)%*%X);
  diag(sigma)=diag(sigma)*SSA;

  #if(det(sigma) == 0) {print(sigma); break}

  logp=log_p_alpha(Alpha,Z,X,B=Beta,SA=Sigma_Alpha,MU=mu);
  if( logp == Inf | logp == -Inf | is.na(logp)) {alpha=Alpha}
  else {alpha <-optim(Alpha,Z=Z,W=X,B=Beta,SA=Sigma_Alpha,MU=mu,log_p_alpha,method="Nelder-Mead",
                      hessian=FALSE,control=list(maxit=10,fnscale=-1))$par}

  Alpha_star=t(rmvt(n=Multiple,alpha,sigma=sigma,df=df))

  for (i in 1:Multiple){
    log_p_Alpha_star[i]=log_p_alpha(Alpha_star[,i],Z,X,Beta,Sigma_Alpha,mu)
  }

  #control overfloat, -max(log_p_Alpha_star);
  p=exp(log_p_Alpha_star-max(log_p_Alpha_star))/sum(exp(log_p_Alpha_star - max(log_p_Alpha_star)))

  #in case there is still overfloat;
  p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
  j=sample(c(1:Multiple),1,prob=p);

  Alpha_2star=t(rmvt(n=Multiple-1,Alpha_star[,j],sigma=sigma,df=df))
  Alpha_2star <-cbind(Alpha_2star,Alpha)

  for (i in 1:Multiple){
    log_p_Alpha_2star[i]=log_p_alpha(Alpha_2star[,i],Z,X,Beta,Sigma_Alpha,mu)
  }

  #control overfloat
  num=sum(exp(log_p_Alpha_star -max(log_p_Alpha_star)))
  den=sum(exp(log_p_Alpha_2star -max(log_p_Alpha_star)))

  rho=min(1,num/den)

  #in case overfloat again
  if(is.na(rho)) {rho=0.5};

  accp=0;
  u=runif(1)
  if(u<rho){
    Alpha=Alpha_star[,j]
    accp=1;
  }


  return(list(par=Alpha, accp=accp))
}

log_p_gamma=function(Gamma,Z,X,SG,phi){
  log_p=sum(phi*(X%*%Gamma)-log(1+exp(X%*%Gamma)))-0.5*t(Gamma)%*%solve(SG)%*%Gamma
  return(log_p)
}

##Draw Gamma, log scale, multiple-try MH;
Draw_Gamma_log_M=function(Gamma,Z,X,Phi,df,Sigma_Gamma,SSG,Multiple)

{
  log_p_Gamma_star=rep(0,Multiple)
  log_p_Gamma_2star=rep(0,Multiple)
  p=rep(0,Multiple)

  #if(det(sigma) == 0) {print(sigma); break}

  gamma_opt <-optim(Gamma,Z=Z,X=X,SG=Sigma_Gamma,phi=Phi,log_p_gamma,method="Nelder-Mead",
                    hessian=TRUE,control=list(maxit=10,fnscale=-1))

  gamma=gamma_opt$par
  sigma=solve(-gamma_opt$hessian)
  diag(sigma)=diag(sigma)*SSG

  Gamma_star=t(rmvt(n=Multiple,gamma,sigma=sigma,df=df))

  for (i in 1:Multiple){
    log_p_Gamma_star[i]=log_p_gamma(Gamma_star[,i],Z,X,Sigma_Gamma,Phi)

  }


  #control overfloat, -max(log_p_Gamma_star);
  p=exp(log_p_Gamma_star-max(log_p_Gamma_star))/sum(exp(log_p_Gamma_star - max(log_p_Gamma_star)))

  #in case there is still overfloat;
  p[is.na(p)]=(1-sum(p[!is.na(p)]))/sum(is.na(p))
  j=sample(c(1:Multiple),1,prob=p);


  Gamma_2star=t(rmvt(n=Multiple-1,Gamma_star[,j],sigma=sigma,df=df))
  Gamma_2star <-cbind(Gamma_2star,Gamma)

  for (i in 1:Multiple){
    log_p_Gamma_2star[i]=log_p_gamma(Gamma_2star[,i],Z,X,Sigma_Gamma,Phi)
  }


  #control overfloat
  num=sum(exp(log_p_Gamma_star -max(log_p_Gamma_star)))
  den=sum(exp(log_p_Gamma_2star -max(log_p_Gamma_star)))

  rho=min(1,num/den)

  #in case overfloat again
  if(is.na(rho)) {rho=0.5};


  accp=0;
  u=runif(1)
  if(u<rho){
    Gamma=Gamma_star[,j]
    accp=1;

  }


  return(list(par=Gamma, accp=accp))
}


##### Function to measure performance #####
library(pROC)

performance <- function(true, est, fdr){

  fdr <- as.numeric(fdr)
  t <- table(true, est)
  if (sum(est)==0){
    FDR.fit <- 0
    power   <- 0
  }
  else if (sum(est)==length(est)){
    FDR.fit <- t[1]/(t[1]+t[2])
    power   <- 1
  }
  else{
    FDR.fit <- t[1,2]/(t[1,2]+t[2,2])
    power   <- t[2,2]/(t[2,1]+t[2,2])
  }
  AUC <- as.numeric(roc(true, fdr)$auc)
  pAUC <- as.numeric(roc(true, fdr, partial.auc = c(1, 0.8), parial.auc.correct = TRUE)$auc)

  return( list( FDR.fit = FDR.fit, power = power, AUC = AUC, pAUC = pAUC))
}
