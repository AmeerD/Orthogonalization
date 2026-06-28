library(dplyr)
library(ggplot2)
library(extraDistr)
library(tidyr)
library(argparse)
library(changepoint)
library(DescTools)
library(purrr)

options(dplyr.summarise.inform = FALSE)

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--nreps", type = "double", default = 50,
                    help = "number of replicates for each set of params")
parser$add_argument("--nU", type = "double", default = 5000,
                    help = "number of U samples to generate")
args <- parser$parse_args()
jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(jobid)

nreps <- args$nreps
nU0 <- args$nU
nfold <- 10

nu <- 5

set.seed(jobid)

filename <- paste("res/CP_", jobid, ".txt", sep="")

## -----------------------------------------
## Simulation
## -----------------------------------------

wvar <- c(0.5, 1, 2) 
nps <- c(100, 250, 500) 
dlist <- seq(0, 2.5, by=0.25)
gfun <- c("segment", "identity")

for (i in 1:nreps) {
  for (ns in nps) {
    for (d in dlist) {
      mu <- c(-d, d) + 2
      
      ## Generate null data
      X <- c(mu[1] + rt(ns, df=nu), mu[2] + rt(ns, df=nu)) 
      museq <- c(rep(mu[1], ns), rep(mu[2], ns))
      
      for (wv in wvar) {
        ## Introduce noise and generate U replicates
        W <- rnorm(length(X), sd=sqrt(wv))
        u0 <- data.frame(u1 = rnorm(nU0, sd=sqrt(wv)))
        
        X1 <- X + W
        X2 <- X - W
        
        cphat <- cpts(cpt.mean(X1, method="BinSeg", penalty="None", Q=4, minseglen=30))
        cphat <- cphat[which(cphat > 0)]
        
        seghat <- c(cphat, 2*ns) - c(0, cphat)
        nullseg <- seghat[1:(length(seghat)-1)] + seghat[2:length(seghat)]
        ntest <- max(nullseg)
        nullidx <- which(nullseg == ntest)[1]
        
        segidx <- rep(1:length(seghat), seghat)
        Delta <- sqrt((mean(museq[segidx == nullidx]) - mean(museq[segidx == (nullidx + 1)]))^2)
        dprop <- sum(segidx == nullidx)/sum(segidx %in% c(nullidx, nullidx + 1))
        
        ## Subset data for testing
        Xsub <- X[segidx %in% c(nullidx, nullidx + 1)]
        X1sub <- X1[segidx %in% c(nullidx, nullidx + 1)]
        X2sub <- X2[segidx %in% c(nullidx, nullidx + 1)]
        
        cfolds <- rep(1:nfold, length.out=ntest)
        mucf <- map(1:nfold, function(x){mean(Xsub[which(cfolds != x)])}) %>% unlist

        ## Compute conditional mean and as much of the IF as possible without g and IFlambda
        cmean <- merge(data.frame(x1=X1sub, fold=cfolds) %>% mutate(idx = row_number()), u0, by=NULL) %>%
          mutate(muhat = mucf[fold]) %>%
          group_by(idx, x1) %>%
          mutate(xdens = dt(u1+x1-muhat, df=nu),
                 xgrad = xdens*(nu+1)*(u1+x1-muhat)/(nu+(u1+x1-muhat)^2)) %>%
          summarise(N = mean((u1+x1)*xdens), D = mean(xdens), 
                    EUnabla = mean((u1+x1)*xgrad), Enabla = mean(xgrad)) %>%
          ungroup %>%
          mutate(cm = 2*N/D - x1, 
                 IFpartial = (EUnabla - (N/D)*Enabla)/D) %>%
          select(-idx)

        Orth <- data.frame(orth=X2sub - cmean$cm)
        
        ## Compute quantities for conditonal correction
        IFtheta <- Xsub - mucf[cfolds]
        EIFtheta <- cmean$N/cmean$D - mucf[cfolds]
        
        ## Compute B using a fresh sample as if mucf is the true value
        Btemp <- map(1:nfold, function(x){
          merge(data.frame(x1=mucf[x] + rt(ntest, df=nu)) %>% mutate(idx = row_number()), u0, by=NULL) %>%
            mutate(muhat = mucf[x]) %>%
            group_by(idx, x1) %>%
            mutate(xdens = dt(u1+x1-muhat, df=nu),
                   xgrad = xdens*(nu+1)*(u1+x1-muhat)/(nu+(u1+x1-muhat)^2)) %>%
            summarise(N = mean((u1+x1)*xdens), D = mean(xdens),
                      EUnabla = mean((u1+x1)*xgrad), Enabla = mean(xgrad)) %>%
            ungroup %>%
            mutate(IFpartial = (EUnabla - (N/D)*Enabla)/D)
        })
        Bmc <- map(Btemp, function(x){
          x %>%
            summarise(B = mean(IFpartial) - 1) %>% 
            pull(B)
        }) %>% unlist
        
        ## Compute conditional mean with true parameters as a sanity check
        cmean.known <- merge(data.frame(x1=X1sub, muhat=museq[segidx %in% c(nullidx, nullidx + 1)]) %>% mutate(idx = row_number()), u0, by=NULL) %>%
          mutate(muhat = muhat) %>%
          group_by(idx, x1) %>%
          mutate(xdens = dt(u1+x1-muhat, df=nu),
                 xgrad = xdens*(nu+1)*(u1+x1-muhat)/(nu+(u1+x1-muhat)^2)) %>%
          summarise(N = mean((u1+x1)*xdens), D = mean(xdens),
                    EUnabla = mean((u1+x1)*xgrad), Enabla = mean(xgrad)) %>%
          ungroup %>%
          mutate(cm = 2*N/D - x1,
                 IFpartial = (EUnabla - (N/D)*Enabla)/D) %>%
          select(-idx)

        Orth.known <- data.frame(orth=X2sub - cmean.known$cm)
        
        for (gf in gfun) {
          if (gf == "segment") {
            g <- matrix(as.numeric(segidx[segidx %in% c(nullidx, nullidx + 1)] > nullidx, ncol=1))
            Amc <- map(Btemp, function(x){
              bind_cols(x, g=g[,1]) %>%
                summarise(A = mean(IFpartial*g)) %>% 
                pull(A)
            }) %>% unlist
          } else if (gf == "identity") {
            g <- matrix(X1sub, ncol=1)
            Amc <- map(Btemp, function(x){
              x %>%
                summarise(A = mean(IFpartial*x1)) %>% 
                pull(A)
            }) %>% unlist
          } else {
            g <- rep(0, length(X1sub))
            Amc <- rep(0, nfold)
          }
          
          ntest <- nrow(Orth)
          tdim <- ncol(g)
          
          ## Conduct orthogonalization test with the true conditional mean
          tautemp.known <- matrix(Orth.known$orth, ncol=1, byrow=F) * g
          Cn <- colMeans(tautemp.known)
          tauhat <- var(tautemp.known)
          Torth <- ntest*Cn^2/tauhat
          porth <- pchisq(Torth, df=1, lower.tail=F)
          write(c(jobid, i, d, 2*ns, nU0, "known", Delta, wv, gf, dprop, ntest, Torth, porth, NA), file = filename, append=TRUE, ncolumns=14)
          
          ## Conduct conditional orthogonalization test
          tautemp <- matrix(Orth$orth, ncol=1, byrow=F) * g
          
          Dn <- colMeans(tautemp + 2*(Amc[cfolds]/(1+Bmc[cfolds]))*EIFtheta)
          IFDn <- tautemp - 2*(Amc[cfolds]/(1+Bmc[cfolds]))*(IFtheta - EIFtheta)
          tauhat1 <- var(IFDn[which(segidx[segidx %in% c(nullidx, nullidx+1)] == nullidx),])
          tauhat2 <- var(IFDn[which(segidx[segidx %in% c(nullidx, nullidx+1)] != nullidx),])
          tauhat <- dprop*tauhat1 + (1-dprop)*tauhat2
          Torth <- ntest*Dn^2/tauhat
          porth <- pchisq(Torth, df=1, lower.tail=F)
          write(c(jobid, i, d, 2*ns, nU0, "conditional", Delta, wv, gf, dprop, ntest, Torth, porth, mean(Bmc)), file = filename, append=TRUE, ncolumns=14)
        } 
      }
    }
  }
}











