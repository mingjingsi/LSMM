##### detect associated SNPs #####
assoc.SNP <- function(fit, FDRset = 0.1, fdrControl){

  M <- length(fit$pi1)
  
  if(length(fit) == 6){
    est              <- NULL
    est$gamma.stage1 <- rep(0, M)
    est$gamma        <- rep(0, M)
    
    if (fdrControl == "global"){
      FDR.stage1 <- post2FDR(fit$pi1.stage1)
      FDR        <- post2FDR(fit$pi1)
      
      est$gamma.stage1[which(FDR.stage1 <= FDRset)] <- 1
      est$gamma[which(FDR <= FDRset)]               <- 1
    }
    if (fdrControl == "local"){
      est$gamma.stage1[which((1-fit$pi1.stage1) <= FDRset)] <- 1
      est$gamma[which((1-fit$pi1) <= FDR.set)]              <- 1
    }
  }
  else {
    est              <- NULL
    est$gamma.stage1 <- rep(0, M)
    est$gamma.stage2 <- rep(0, M)
    est$gamma        <- rep(0, M)
    
    if (fdrControl == "global"){
      FDR.stage1 <- post2FDR(fit$pi1.stage1)
      FDR.stage2 <- post2FDR(fit$pi1.stage2)
      FDR        <- post2FDR(fit$pi1)
      
      est$gamma.stage1[which(FDR.stage1 <= FDRset)] <- 1
      est$gamma.stage2[which(FDR.stage2 <= FDRset)] <- 1
      est$gamma[which(FDR <= FDRset)]               <- 1
    }
    if (fdrControl == "local"){
      est$gamma.stage1[which((1-fit$pi1.stage1) <= FDRset)] <- 1
      est$gamma.stage2[which((1-fit$pi1.stage2) <= FDRset)] <- 1
      est$gamma[which((1-fit$pi1) <= FDR.set)]              <- 1
    }
  }

  return(est)
}

##### detect relevant annotations #####
relev.Anno <- function(fit, FDRset = 0.1, fdrControl){

  K   <- length(fit$omegak)
  est <- rep(0, K)

  if (fdrControl == "global"){
    FDR <- post2FDR(fit$omegak)
    est[which(FDR <= FDRset)] <- 1
  }
  if (fdrControl == "local"){
    est[which((1-fit$omegak) <= FDRset)] <- 1
  }

  return(est)
}

##### transform posterior to FDR #####
post2FDR <- function(posterior){

  M          <- length(posterior)
  fdr        <- 1 - posterior
  rank.fdr   <- rank(fdr)
  sort.fdr   <- sort(fdr)
  cumsum.fdr <- cumsum(sort.fdr)
  sort.FDR   <- cumsum.fdr/seq(1, M, 1)
  FDR        <- sort.FDR[rank.fdr]

  return(FDR)
}
