library(dplyr)
library(ggplot2)
library(extraDistr)
library(tidyr)
library(argparse)
library(changepoint)
library(DescTools)
library(purrr)
library(numDeriv)

options(dplyr.summarise.inform = FALSE)

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--nreps", type = "double", default = 100,
                    help = "number of replicates for each set of params")
parser$add_argument("--p", type = "double", default = 2,
                    help = "dimension of X")
parser$add_argument("--nU", type = "double", default = 5000,
                    help = "number of U samples to generate")
args <- parser$parse_args()
jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(jobid)

nreps <- args$nreps
p <- args$p
nfold <- 5
nclust <- 3
nmc <- args$nU*p

source("helpers.R")

set.seed(jobid)

filename <- paste("res/zipclust_", p, "_", jobid, ".txt", sep="")

## -----------------------------------------
## Simulation
## -----------------------------------------

wvar <- 3
nps <- p*c(50, 125, 250)
pi <- 0.2+0.05*(1:p) %% 3 #rep_len(c(0.2, 0.3, 0.25), p)
lambda <- 3 + (1:p) %% 3 #rep_len(c(4, 3, 5), p)
dlist <- seq(0, 4, by=1)
gfun <- c("cluster", "l2", "sup")

for (wv in wvar) {
  for (d in dlist) {
    ## Set number of replicates to one for non-null simulations
    nreps2 <- ifelse(d > 0, 1, nreps)
    
    condtab1 <- map(1:p, function(x){
      pi0 <- pi[x]
      l0 <- lambda[x]
      
      expand_grid(x1=(0-wv):(qzip(0.999999, l0, pi0)+2*wv), x2=(0-wv):(qzip(0.999999, l0, pi0)+2*wv)) %>% 
        mutate(cprob = dzip((x1+x2)/2, l0, pi0)*ddunif((x1-x2)/2, -wv, wv)) %>% 
        filter(cprob > 0) %>% 
        group_by(x1) %>% 
        mutate(cprob = cprob/sum(cprob),
               ccdf = cumsum(cprob)) %>% 
        ungroup
    })
    condtab2 <- map(1:p, function(x){
      pi0 <- pi[x]
      l0 <- lambda[x]+d
      
      expand_grid(x1=(0-wv):(qzip(0.999999, l0, pi0)+2*wv), x2=(0-wv):(qzip(0.999999, l0, pi0)+2*wv)) %>% 
        mutate(cprob = dzip((x1+x2)/2, l0, pi0)*ddunif((x1-x2)/2, -wv, wv)) %>% 
        filter(cprob > 0) %>% 
        group_by(x1) %>% 
        mutate(cprob = cprob/sum(cprob),
               ccdf = cumsum(cprob)) %>% 
        ungroup
    })
    
    for (ns in nps) {
      print(paste0("Distance: ", d, "; Sample size: ", 2*ns))
      ## Generate X1 data and identify clusters
      X1 <- map(1:p, function(x){
        temp <- c(rzip(ns, lambda[x], pi[x]), rzip(ns, lambda[x]+d, pi[x]))
        temp + rdunif(2*ns, -wv, wv)
      }) %>% do.call(cbind, .) 
      museq <- map(1:p, function(x){
        c(rep((1-pi[x])*lambda[x], ns), rep((1-pi[x])*(lambda[x]+d), ns))
      }) %>% do.call(cbind, .)
      
      clust <- kmeans(X1, centers=nclust, nstart=100)
      Z <- clust$cluster
      ctrs <- clust$centers
      
      top2 <- sort(table(Z), decreasing=T)[1:2]
      dprop <- top2[1]/sum(top2)
      ntest <- sum(top2)
      
      mu1 <- colMeans(museq[which(Z == names(top2)[1]),])
      mu2 <- colMeans(museq[which(Z == names(top2)[2]),])
      Delta <- sqrt(sum((mu1-mu2)^2))
      
      ## Subset X1 based on clustering
      X1sub <- X1[which(Z %in% names(top2)),]
      idxsub <- c(rep(0, ns), rep(1, ns))[which(Z %in% names(top2))]
      
      ## Compute E[X2|X1] for this draw of X1
      ctemp.known <- map(1:ntest, function(x) {
        ## Implement each dimension separately to save on computation - follows from independence of dimensions
        marg <- map(1:p, function(y){
          Utemp <- X1sub[x,y] + (-wv:wv)
          fX <- dzip(Utemp, lambda[y]+d*idxsub[x], pi[y])
          N <- mean(Utemp*fX)
          D <- mean(fX)
          cmean <- 2*N/D - X1sub[x,y]
          
          sc <- zipscore(Utemp, pi[y], lambda[y])
          EUnabla <- apply(sc, 2, function(x){mean(Utemp*x*fX)})
          Enabla <- apply(sc, 2, function(x){mean(x*fX)})
          
          IFpartial <- (EUnabla - t(t(N/D)) %*% Enabla)/D
          list(cmean, IFpartial)
        })
        
        list(map(marg, ~ .x[[1]]) %>% unlist,
             map(marg, ~ .x[[2]]) %>% Matrix::bdiag())
      })
      
      cmean.known <- do.call(rbind, map(ctemp.known, ~.x[[1]]))
      
      for (i in 1:nreps2) {
        ## Generate X2 using X2|X1
        X2 <- map(1:p, function(x){
          c(data.frame(x1=X1[1:ns, x], q=runif(ns)) %>% 
              mutate(idx = row_number()) %>% 
              left_join(condtab1[[x]], by=c("x1"), relationship="many-to-many") %>% 
              filter(ccdf > q) %>% 
              group_by(idx) %>% 
              filter(ccdf == min(ccdf)) %>% 
              pull(x2),
            data.frame(x1=X1[(ns+1):(2*ns), x], q=runif(ns)) %>% 
              mutate(idx = row_number()) %>% 
              left_join(condtab2[[x]], by=c("x1"), relationship="many-to-many") %>% 
              filter(ccdf > q) %>% 
              group_by(idx) %>% 
              filter(ccdf == min(ccdf)) %>% 
              pull(x2))
        }) %>% do.call(cbind, .)
        
        ## Subset X2 for estimation and compute X
        X2sub <- X2[which(Z %in% names(top2)),]
        Xsub <- (X1sub + X2sub)/2
        
        cfolds <- rep(1:nfold, length.out=ntest)
        ## Compute an initial estimate ignoring truncation for seeding the MLE
        ## and drawing samples for importance sampling
        naive <- map(1:nfold, function(x){
          starting <- map(1:p, function(y){
            optim(c(0.5, 3), zipll, dat=Xsub[which(cfolds != x),y],
                  lower=rep(1e-3,2), upper=c(1-1e-3,Inf), method="L-BFGS-B",
                  control=list(fnscale = -1), hessian=F)$par
          }) %>% do.call(rbind, .)
        })

        ## Secondary stage in which we account for the truncation
        parscf <- map(1:nfold, function(x){
          ## Start with an initial truncated estimated with a Taylor approximation to the truncation term
          naivepars <- naive[[x]]
          ## Draw a sample from the distribution implied by the naive estimate
          X1naive <- map(1:p, function(y){
            rzip(nmc, naivepars[y,2], naivepars[y,1]) + extraDistr::rdunif(nmc, -wv, wv)
          }) %>% do.call(cbind, .)
          Asamp <- X1naive[which(max.col(-as.matrix(pdist::pdist(X1naive, ctrs))) %in% names(top2)),]
          Alogdens <- map(1:p, ~log(dzip1(Asamp[,.x], naivepars[.x,1], naivepars[.x,2], wv))) %>%
            do.call(cbind, .) %>%
            apply(1, sum)
          ## Compute components of the linearization step
          Anorm <- nrow(Asamp)/nmc
          Ascore <- map(1:p, ~zipscore1(Asamp[,.x], naivepars[.x,1], naivepars[.x,2], wv)) %>%
            do.call(cbind, .)
          Agrad <- (colSums(Ascore)/nmc)/Anorm
          ## Compute the approximate MLE
          temp <- optim(as.vector(t(naivepars)), zipcll3TE, dat=Xsub[which(cfolds != x),],
                        base=as.vector(t(naivepars)), tnbase=log(Anorm), tngrad=Agrad,
                        lower=rep(1e-3,2*p), upper=rep(c(1-1e-3,Inf), p), method="L-BFGS-B",
                        control=list(fnscale = -1, trace=1, REPORT=1), hessian=F)
          taylor <- matrix(temp$par, ncol=2, byrow=T)

          ## Redo with the approximate MLE with a better starting point
          ## Draw a sample from the distribution implied by the Taylor approximation estimate
          X1taylor <- map(1:p, function(y){
            rzip(nmc, taylor[y,2], taylor[y,1]) + extraDistr::rdunif(nmc, -wv, wv)
          }) %>% do.call(cbind, .)
          Asamp <- X1taylor[which(max.col(-as.matrix(pdist::pdist(X1taylor, ctrs))) %in% names(top2)),]
          Alogdens <- map(1:p, ~log(dzip1(Asamp[,.x], taylor[.x,1], taylor[.x,2], wv))) %>%
            do.call(cbind, .) %>%
            apply(1, sum)
          Anorm <- nrow(Asamp)/nmc
          Ascore <- map(1:p, ~zipscore1(Asamp[,.x], taylor[.x,1], taylor[.x,2], wv)) %>%
            do.call(cbind, .)
          Agrad <- (colSums(Ascore)/nmc)/Anorm
          ## Compute the approximate MLE
          temp <- optim(as.vector(t(taylor)), zipcll3TE, dat=Xsub[which(cfolds != x),],
                        base=as.vector(t(taylor)), tnbase=log(Anorm), tngrad=Agrad,
                        lower=rep(1e-3,2*p), upper=rep(c(1-1e-3,Inf), p), method="L-BFGS-B",
                        control=list(fnscale = -1, trace=1, REPORT=1), hessian=F)

          # ## Compute the MLE - DIRECT OPTIMIZATION NOT FEASIBLE IN HIGH DIMENSIONS
          # temp <- optim(as.vector(t(taylor)), zipcll3, dat=Xsub[which(cfolds != x),],
          #               Asamp=Asamp, Anaive=Alogdens, sig=wv, nMC=nmc,
          #               lower=rep(1e-4,2*p), upper=rep(c(1-1e-4,Inf), p), method="L-BFGS-B",
          #               control=list(fnscale = -1), hessian=F)
          
          temp$par <- matrix(temp$par, ncol=2, byrow=T)

          temp$sX <- map(1:p, function(y){
            zipscore(Xsub[which(cfolds != x),y], temp$par[y,1], temp$par[y,2])
          }) %>% do.call(cbind, .)
          temp$score.recentre <- colMeans(temp$sX)
          # temp$score.recentre <- grad(truncnorm3, c(t(temp$par)), Asamp=Asamp, Anaive=Alogdens, sig=wv, nMC=nmc)

          obsinf <- map(1:p, function(y){
            apply(ziphess(Xsub[which(cfolds != x),y], temp$par[y,1], temp$par[y,2]), 2:3, mean)
          }) %>% Matrix::bdiag() %>% as.matrix
          # temp$hess.shift <- hessian(truncnorm3, c(t(temp$par)), Asamp=Asamp, Anaive=Alogdens, sig=wv, nMC=nmc)
          # J <- obsinf - temp$hess.shift
          J <- -crossprod(sweep(temp$sX, 2, temp$score.recentre, FUN="-"))/nrow(Xsub[which(cfolds != x),])
          temp$hess.shift <- obsinf - J
          temp$Jinv <- solve(-J)

          temp
        })

        ## Estimate conditional mean and various intermediate quantities for each unit
        ctemp <- map(1:ntest, function(x) {
          phat <- parscf[[cfolds[x]]]$par

          ## Implement each dimension separately to save on computation - follows from independence of dimensions
          marg <- map(1:p, function(y){
            Utemp <- X1sub[x,y] + (-wv:wv)
            fX <- dzip(Utemp, phat[y,2], phat[y,1])
            N <- mean(Utemp*fX)
            D <- mean(fX)
            cmean <- 2*N/D - X1sub[x,y]

            sc <- zipscore(Utemp, phat[y,1], phat[y,2])
            EUnabla <- apply(sc, 2, function(x){mean(Utemp*x*fX)})
            Enabla <- apply(sc, 2, function(x){mean(x*fX)})

            IFpartial <- (EUnabla - t(t(N/D)) %*% Enabla)/D

            ## Compute an estimate of the influence function for each unit WITHOUT RECENTRING
            stilde <- zipscore(Xsub[x,y], phat[y,1], phat[y,2])

            ## Compute an estimate of the expected value of the influence function for each unit WITHOUT RECENTRING
            sUtilde <- t(Enabla/D)

            ## Compute the gradient of E[IFtheta|X1] for computing B
            cIFgrad1a <- apply(ziphess(Utemp, phat[y,1], phat[y,2]), 2:3, function(x){mean(x*fX)})
            cIFgrad1b <- (t(sc) %*% diag(fX) %*% sc)/length(Utemp)
            cIFgrad2 <- Enabla %*% t(Enabla)
            cIFgrad <- (cIFgrad1a+cIFgrad1b)/D - cIFgrad2/D^2


            list(cmean, IFpartial, stilde, sUtilde, cIFgrad)
          })

          cmean <- map(marg, ~ .x[[1]]) %>% unlist
          IFpartial <- map(marg, ~ .x[[2]]) %>% Matrix::bdiag() %>% as.matrix()
          IFtheta <- t(parscf[[cfolds[x]]]$Jinv %*% t(do.call(cbind, map(marg, ~.x[[3]])) - parscf[[cfolds[x]]]$score.recentre))

          cIFtheta <- t(parscf[[cfolds[x]]]$Jinv %*% t(do.call(cbind, map(marg, ~.x[[4]])) - parscf[[cfolds[x]]]$score.recentre))
          cIFgrad <- parscf[[cfolds[x]]]$Jinv %*% (as.matrix(Matrix::bdiag(map(marg, ~.x[[5]]))) - parscf[[cfolds[x]]]$hess.shift)


          list(cmean, IFpartial, IFtheta, cIFtheta, cIFgrad)
        })

        cmean <- do.call(rbind, map(ctemp, ~.x[[1]]))
        IFtheta <- do.call(rbind, map(ctemp, ~.x[[3]]))
        cIFtheta <- do.call(rbind, map(ctemp, ~.x[[4]]))
        cIFgrad <- abind::abind(map(ctemp, ~.x[[5]]), along=3)
        B <- map(1:nfold, ~apply(cIFgrad[,,which(cfolds != .x)], 1:2, mean))
        Binv <- map(B, ~ solve(diag(2*p) + .x))

        Orth <- X2sub - cmean
        Orth.known <- X2sub - cmean.known
        
        for (gf in gfun) {
          if (gf == "cluster") {
            g <- as.numeric(Z[which(Z %in% names(top2))] == names(top2)[1])
          } else if (gf == "l2") {
            g <- apply(X1sub, 1, function(x){sqrt(sum(x^2))})
          } else if (gf == "sup") {
            g <- apply(X1sub, 1, function(x){max(abs(x))})
          }
          
          ## Conduct orthogonalization test with the true conditional mean
          tautemp.known <- Orth.known * matrix(g, nrow=ntest, ncol=p, byrow=F)
          Cn <- colMeans(tautemp.known)
          tauhat1 <- var(tautemp.known[which(Z[which(Z %in% names(top2))] == names(top2)[1]),])
          tauhat2 <- var(tautemp.known[which(Z[which(Z %in% names(top2))] == names(top2)[2]),])
          tauhat <- dprop*tauhat1 + (1-dprop)*tauhat2
          Torth <- ntest*(Cn %*% solve(tauhat) %*% Cn)
          porth <- pchisq(Torth, df=p, lower.tail=F)
          write(c(jobid, i, 2*ns, p, wv, d, Delta, gf, ntest, "known", Torth, porth), file = filename, append=TRUE, ncolumns=12)
          
          ## Conduct conditional orthogonalization test with estimated conditional IF variance correction and debiasing
          IFpartial2 <- map(1:ntest, ~ctemp[[.x]][[2]]*g[.x]) %>% abind::abind(along=3)
          A <- map(1:nfold, function(x){apply(IFpartial2[,,which(cfolds != x)], 1:2, mean)})

          tautemp <- Orth * matrix(g, nrow=ntest, ncol=p, byrow=F)
          debias <- map(1:ntest, function(x){2*t(A[[cfolds[x]]] %*% Binv[[cfolds[x]]] %*% cIFtheta[x,])}) %>% do.call(rbind, .)
          Dn <- colMeans(tautemp + debias)

          IFdebias <- map(1:ntest, function(x){2*t(A[[cfolds[x]]] %*% Binv[[cfolds[x]]] %*% (IFtheta[x,] - cIFtheta[x,]))}) %>% do.call(rbind, .)
          IFDn <- tautemp - IFdebias

          tauhat1 <- var(IFDn[which(Z[which(Z %in% names(top2))] == names(top2)[1]),])
          tauhat2 <- var(IFDn[which(Z[which(Z %in% names(top2))] == names(top2)[2]),])
          tauhat <- dprop*tauhat1 + (1-dprop)*tauhat2
          Torth <- ntest*(Dn %*% solve(tauhat) %*% Dn)
          porth <- pchisq(Torth, df=p, lower.tail=F)
          write(c(jobid, i, 2*ns, p, wv, d, Delta, gf, ntest, "conditional", Torth, porth), file = filename, append=TRUE, ncolumns=12)
        }
      }
    }
  }
}
