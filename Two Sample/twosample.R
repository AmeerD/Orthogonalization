library(dplyr)
library(ggplot2)
library(mvtnorm)
library(tidyr)
library(purrr)
library(Matrix)
library(argparse)

options(dplyr.summarise.inform = FALSE)

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--nreps", type = "double", default = 50,
                    help = "number of replicates for each set of params")
parser$add_argument("--dim", type="double", default=5,
                    help = "dimension")

args <- parser$parse_args()
jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(jobid)

nreps <- args$nreps
p <- args$dim

set.seed(jobid)

filename <- paste("res/twosamp_", p, "_", jobid, ".txt", sep="")

## -----------------------------------------
## Set data parameters
## -----------------------------------------

# Set data parameters
mu1 <- rep(0, p)
Sig1 <- diag(p)

## -----------------------------------------
## Simulation
## -----------------------------------------

wvar <- c(2, 5, 10) 
nps <- c(250, 500, 1000, 2500) 

gfun <- c("identity", "l2", "first")
targets <- c("mean", "variance", "covariance")

for (target in targets) {
  print(paste0("Starting ", target))
  if (target == "mean") {
    dlist <- seq(0, 4, length.out=5)
  } else if (target == "variance") {
    dlist <- seq(2, 5, length.out=4)
  } else if (target == "covariance") {
    dlist <- seq(0.2, 0.8, length.out=4)
  }
  
  for (i in 1:nreps) {
    for (ns in nps) {
      for (d in dlist) {
        ## Set alternative parameters
        if (target == "mean") {
          mu2 <- rep(d, p)
          Sig2 <- Sig1
        } else if (target == "variance") {
          mu2 <- mu1
          Sig2 <- d*Sig1
        } else if (target == "covariance") {
          mu2 <- mu1
          Sig2 <- matrix(d, nrow=p, ncol=p) + (1-d)*diag(p)
        }
        
        ## Generate data
        X <- data.frame(rmvnorm(ns, mu1, Sig1))
        Y <- data.frame(rmvnorm(ns, mu2, Sig2))
        
        ntest <- nrow(X) + nrow(Y)
        
        for (wv in wvar) {
          # Introduce noise
          WX <- data.frame(rmvnorm(ns, sigma=wv*diag(p)))
          WY <- data.frame(rmvnorm(ns, sigma=wv*diag(p)))
          
          X1 <- X + WX
          X2 <- X - WX
          Y1 <- Y + WY
          Y2 <- Y - WY
          
          rm(WX, WY)
          gc()
          
          # Estimate the conditional means
          ## Start with E[X2|X1] computed using Y observations
          ### Rows correspond to Y and columns correspond to X1
          fR.X <- apply(X1, 1, function(x) {dmvnorm(Y, x, wv*diag(p))})
          N.X <- t(apply(fR.X, 2, function(x) {colMeans(Y*matrix(x, nrow=ns, ncol=p, byrow=F))}))
          D.X <- colMeans(fR.X)
          cmean.X <- 2*N.X/matrix(D.X, nrow=ns, ncol=p, byrow=F) - X1
          ## Next do E[Y2|Y1] computed using X observations
          ### Rows correspond to X and columns correspond to Y1
          fR.Y <- apply(Y1, 1, function(x) {dmvnorm(X, x, wv*diag(p))})
          N.Y <- t(apply(fR.Y, 2, function(x) {colMeans(X*matrix(x, nrow=ns, ncol=p, byrow=F))}))
          D.Y <- colMeans(fR.Y)
          cmean.Y <- 2*N.Y/matrix(D.Y, nrow=ns, ncol=p, byrow=F) - Y1
          
          # Next compute the influence function corrections
          ## Represents errors X makes towards estimating Cn
          IFND.X <- array(NA, dim=c(ns, ns, p)) 
          for (idx in 1:ns) {
            ## For each x, compute the components of the influence function correction that don't involve g
            temp <- (matrix(unlist(X[idx,]), nrow=ns, ncol=p, byrow=T) - (N.Y/matrix(D.Y, nrow=ns, ncol=p, byrow=F)))
            IFND.X[idx,,] <- temp * matrix(fR.Y[idx,]/D.Y, nrow=ns, ncol=p, byrow=F)
          }
          ## Represents errors Y makes towards estimating Cn
          IFND.Y <- array(NA, dim=c(ns, ns, p)) 
          for (idx in 1:ns) {
            ## For each y, compute the components of the influence function correction that don't involve g
            temp <- (matrix(unlist(Y[idx,]), nrow=ns, ncol=p, byrow=T) - (N.X/matrix(D.X, nrow=ns, ncol=p, byrow=F)))
            IFND.Y[idx,,] <- temp * matrix(fR.X[idx,]/D.X, nrow=ns, ncol=p, byrow=F)
          }
          
          Orth <- bind_rows(
            X2 - cmean.X,
            Y2 - cmean.Y
          )
          
          for (gf in gfun) {
            if (gf == "identity") {
              g <- bind_rows(X1, Y1)
            } else if (gf == "l2") {
              g <- data.frame(c(sqrt(rowSums(X1^2)), sqrt(rowSums(Y1^2))))
            } else if (gf == "first") {
              g <- data.frame(bind_rows(X1, Y1)[,1])
            } else {
              g <- data.frame(rep(0, nrow(Orth)))
            }
            
            tdim <- ncol(Orth)*ncol(g)
            
            ## Conduct the chi-squared test with the cross-fitted estimates
            tautemp <- as.matrix(t(KhatriRao(t(Orth), t(g))))
            Cn <- sqrt(ntest) * colMeans(tautemp)
            
            ## Add in g to the influence function correction
            IFcorX <- matrix(NA, nrow=ns, ncol=tdim)
            for (idx in 1:ns) {
              IFcorX[idx,] <- colMeans(t(KhatriRao(t(IFND.X[idx,,]), t(g[(ns+1):(2*ns),]), sparseY=F)))
            }
            
            IFcorY <- matrix(NA, nrow=ns, ncol=tdim)
            for (idx in 1:ns) {
              IFcorY[idx,] <- colMeans(t(KhatriRao(t(IFND.Y[idx,,]), t(g[1:ns,]), sparseY=F)))
            }
            
            IFcor <- rbind(IFcorX, IFcorY)
            
            tauhat <- var(tautemp - 2*IFcor)
            Torth <- Cn %*% solve(tauhat) %*% t(t(Cn))
            porth <- pchisq(Torth, tdim, lower.tail=F)
            
            write(c(jobid, i, p, target, d, ns, wv, gf, "Orth", Torth, porth), file = filename, append=TRUE, ncolumns=11)
          }
          
          rm(X1, X2, Orth, cmean)
          gc()
        }
        
        rm(X)
        gc()
      }
    }
  }
}













