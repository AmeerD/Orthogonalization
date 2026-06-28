dzip1 <- function(x, pi, lambda, sig){
  urange <- -sig:sig
  map(urange, ~dzip(x-.x, lambda, pi)) %>% 
    do.call(cbind, .) %>% 
    rowMeans()
}

zipll <- function(pars, dat) {
  if (pars[1] < 0 | pars[1] > 1 | pars[2] < 0) {
    print(pars)
    return(-100000)
  }
  
  sum(dzip(dat, pars[2], pars[1], log=T))
}

truncnorm3 <- function(pars, Asamp, Anaive, sig, nMC) {
  pmat <- matrix(pars, ncol=2, byrow=T)
  
  pA <- map(1:p, ~log(dzip1(Asamp[,.x], pmat[.x,1], pmat[.x,2], sig))) %>% 
    do.call(cbind, .) %>% 
    apply(1, sum) 
  
  return(log(sum(exp(pA - Anaive))/nMC))
}

## pars must be organized as (pi1, lambda1, pi2, lambda2, ..., pip, lambdap)
zipcll3 <- function(pars, dat, Asamp, Anaive, sig, nMC) {
  pmat <- matrix(pars, ncol=2, byrow=T)
  
  ## First compute P(X1 in clusters) for normalising
  pcl <- truncnorm3(pars, Asamp, Anaive, sig, nMC)
  
  ## Compute log-density for the observed X values
  pX <- map(1:nrow(pmat), function(x){
    dzip(dat[,x], pmat[x,2], pmat[x,1], log=T)
  }) %>% 
    unlist %>% 
    sum
  
  return(pX-nrow(dat)*pcl)
}

## pars must be organized as (pi1, lambda1, pi2, lambda2, ..., pip, lambdap)
zipcll3TE <- function(pars, dat, base, tnbase, tngrad) {
  pmat <- matrix(pars, ncol=2, byrow=T)
  
  ## First compute P(X1 in clusters) for normalising
  diff <- pars - base
  pcl <- tnbase + sum(tngrad*diff)
  pcl <- pmin(pmax(pcl, log(0.01)), log(1))
  
  ## Compute log-density for the observed X values
  pX <- map(1:nrow(pmat), function(x){
    dzip(dat[,x], pmat[x,2], pmat[x,1], log=T)
  }) %>% 
    unlist %>% 
    sum
  
  if (is.infinite(pX-nrow(dat)*pcl)) {
    print(pmat)
    print(pcl)
  }
  
  return(pX-nrow(dat)*pcl)
}

zipscore <- function(x, pi, lambda) {
  z <- as.numeric(x == 0) 
  
  pzero <- pi + (1-pi)*exp(-lambda)
  
  piscore <- ifelse(x < 0, 0, z*(1-exp(-lambda))/pzero - (1-z)/(1-pi))
  lscore <- ifelse(x < 0, 0, -z*(1-pi)*exp(-lambda)/pzero + (1-z)*(x/lambda - 1))
  
  return(cbind(piscore, lscore))
}

zipscore1 <- function(x, pi, lambda, sig){
  urange <- -sig:sig
  uweight <- map(urange, ~dzip(x-.x, lambda, pi)) %>% 
    do.call(cbind, .) %>% 
    sweep(1, rowSums(.), FUN="/")
  uscores <- map(urange, ~zipscore(x-.x, pi, lambda))
  map(1:length(urange), ~sweep(uscores[[.x]], 1, uweight[,.x], "*")) %>% 
    abind::abind(along=3) %>% 
    apply(1:2, sum)
}

ziphess <- function(x, pi, lambda) {
  z <- as.numeric(x == 0) 
  
  pzero <- pi + (1-pi)*exp(-lambda)
  
  d2pi <- ifelse(x < 0, 0, -z*((1-exp(-lambda))/pzero)^2 - (1-z)/((1-pi)^2))
  d2l <- ifelse(x < 0, 0, z*pi*(1-pi)*exp(-lambda)/(pzero^2) - (1-z)*x/(lambda^2))
  dpidl <- ifelse(x < 0, 0, z*exp(-lambda)/(pzero^2))
  
  res <- array(NA, dim=c(length(x), 2, 2))
  res[,1,1] <- d2pi
  res[,1,2] <- res[,2,1] <- dpidl
  res[,2,2] <- d2l
  
  return(res)
}