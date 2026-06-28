library(dplyr)
library(tidyr)
library(ggplot2)
library(extraDistr)
library(scry)
library(DescTools)
library(datathin)
library(purrr)
library(latex2exp)
library(GGally)
library(matrixStats)
library(geometry)
library(kableExtra)
library(DuoClustering2018)
library(patchwork)

set.seed(2026)
fpath <- "Data Analysis/"
source("Clustering/helpers.R")

nclust <- 4
nkeep <- 50
nfold <- 5

nmc <- nkeep*1000

# Load data
sce <- sce_full_Zhengmix4eq()
## Extract useful items
umi <- counts(sce) # UMI count matrix
cm <- as.data.frame(colData(sce)) ## cell data
nUMI <- cm$total_counts #Total number of UMIs per cell
ncells <- ncol(umi)
rm(sce)
gc()

# Orthogonalization
## Baseline pre-processing - subset to the nkeep genes with the highest variance
umi.var <- apply(umi, 1, var)
varbound <- sort(umi.var, decreasing=T)[nkeep]
umi.orth <- umi[which(umi.var >= varbound),] %>% t

## Add and subtract discrete uniform noise to construct train/test data
wv <- ceiling((sqrt(12*varbound+1)-1)/2/2)
W <- matrix(rdunif(ncells*nkeep, -wv, wv), nrow=ncells)

X1 <- umi.orth + W 
X2 <- umi.orth - W

## Selection
X1ra <- sweep(X1, 1, rowSums(X1), "/")
clust <- kmeans(X1ra, centers=nclust, nstart=100)
ctrs <- clust$centers
ct.true <- cm$phenoid
ct.est <- paste0("Cluster ", clust$cluster)
ctf.true <- factor(ct.true, levels=c(unique(ct.true), unique(ct.est)))
ctf.est <- factor(ct.est, levels=c(unique(ct.true), unique(ct.est)))
ctab <- caret::confusionMatrix(ctf.est, ctf.true)$table
ctab2 <- ctab[rowSums(ctab) > 0, colSums(ctab) > 0]
ctab2

## Testing
pmat <- matrix(NA, nrow=nclust, ncol=nclust)
parhat <- array(NA, dim=c(nclust, nclust, nfold, nkeep, 2))
orthlist <- list()
for (i in 1:nclust) {
  for (j in 1:nclust) {
    if (i > j) {
      ## Filter X to data in the grid
      Xsub <- umi.orth[which(clust$cluster %in% c(i,j)),]
      ## Filter data to the top two clusters
      X1sub <- X1[which(clust$cluster %in% c(i,j)),]
      X2sub <- X2[which(clust$cluster %in% c(i,j)),]
      
      ntest <- nrow(Xsub)
      cfolds <- rep(1:nfold, length.out=ntest)
      ## Compute an initial estimate ignoring truncation for seeding the MLE 
      ## and drawing samples for importance sampling
      naive <- map(1:nfold, function(x){
        starting <- map(1:nkeep, function(y){
          Xcf <- Xsub[which(cfolds != x),y]
          init <- c(mean(Xcf == 0), mean(Xcf)/(1-mean(Xcf == 0)))
          optim(init, zipll, dat=Xcf, 
                lower=rep(1e-4,2), upper=c(1-1e-4,Inf), method="L-BFGS-B",
                control=list(fnscale = -1), hessian=F)$par
        }) %>% do.call(rbind, .)
      })
      
      ## Secondary stage in which we account for the truncation
      parscf <- map(1:nfold, function(x){
        print(paste0(i, "-", j, "; fold ", x))
        naivepars <- naive[[x]]
        ## Draw a sample from the distribution implied by the naive estimate
        X1naive <- map(1:nkeep, function(y){
          rzip(nmc, naivepars[y,2], naivepars[y,1]) + extraDistr::rdunif(nmc, -wv, wv)
        }) %>% do.call(cbind, .)
        Asamp <- X1naive[which(max.col(-as.matrix(pdist::pdist(X1naive, ctrs))) %in% c(i,j)),]
        Alogdens <- map(1:nkeep, ~log(dzip1(Asamp[,.x], naivepars[.x,1], naivepars[.x,2], wv))) %>%
          do.call(cbind, .) %>%
          apply(1, sum)
        ## Compute components of the linearization step
        Anorm <- nrow(Asamp)/nmc
        Ascore <- map(1:nkeep, ~zipscore1(Asamp[,.x], naivepars[.x,1], naivepars[.x,2], wv)) %>%
          do.call(cbind, .)
        Agrad <- (colSums(Ascore)/nmc)/Anorm
        ## Compute the approximate MLE
        temp <- optim(as.vector(t(naivepars)), zipcll3TE, dat=Xsub[which(cfolds != x),],
                      base=as.vector(t(naivepars)), tnbase=log(Anorm), tngrad=Agrad,
                      lower=rep(1e-3,2*nkeep), upper=rep(c(1-1e-3,Inf), nkeep), method="L-BFGS-B",
                      control=list(fnscale = -1, trace=1, REPORT=10, maxit=200), hessian=F)
        taylor <- matrix(temp$par, ncol=2, byrow=T)
        
        ## Redo with the approximate MLE with a better starting point
        ## Draw a sample from the distribution implied by the Taylor approximation estimate
        X1taylor <- map(1:nkeep, function(y){
          rzip(nmc, taylor[y,2], taylor[y,1]) + extraDistr::rdunif(nmc, -wv, wv)
        }) %>% do.call(cbind, .)
        Asamp <- X1taylor[which(max.col(-as.matrix(pdist::pdist(X1taylor, ctrs))) %in% c(i,j)),]
        Alogdens <- map(1:nkeep, ~log(dzip1(Asamp[,.x], taylor[.x,1], taylor[.x,2], wv))) %>%
          do.call(cbind, .) %>%
          apply(1, sum)
        Anorm <- nrow(Asamp)/nmc
        Ascore <- map(1:nkeep, ~zipscore1(Asamp[,.x], taylor[.x,1], taylor[.x,2], wv)) %>%
          do.call(cbind, .)
        Agrad <- (colSums(Ascore)/nmc)/Anorm
        ## Compute the approximate MLE
        temp <- optim(as.vector(t(taylor)), zipcll3TE, dat=Xsub[which(cfolds != x),],
                      base=as.vector(t(taylor)), tnbase=log(Anorm), tngrad=Agrad,
                      lower=rep(1e-3,2*nkeep), upper=rep(c(1-1e-3,Inf), nkeep), method="L-BFGS-B",
                      control=list(fnscale = -1, trace=1, REPORT=10, maxit=200), hessian=F)
        temp$par <- matrix(temp$par, ncol=2, byrow=T)
        
        
        temp$sX <- map(1:nkeep, function(y){
          zipscore(Xsub[which(cfolds != x),y], temp$par[y,1], temp$par[y,2])
        }) %>% do.call(cbind, .)
        temp$score.recentre <- colMeans(temp$sX)
        
        obsinf <- map(1:nkeep, function(y){
          apply(ziphess(Xsub[which(cfolds != x),y], temp$par[y,1], temp$par[y,2]), 2:3, mean)
        }) %>% Matrix::bdiag() %>% as.matrix
        J <- -crossprod(sweep(temp$sX, 2, temp$score.recentre, FUN="-"))/nrow(Xsub[which(cfolds != x),])
        temp$hess.shift <- obsinf - J
        temp$Jinv <- pracma::pinv(-J)
        
        temp
      })
      for (x in 1:nfold) {
        parhat[i,j,x,,] <- parscf[[x]]$par
      }
      
      
      ## Estimate conditional mean and various intermediate quantities for each unit
      ctemp <- map(1:ntest, function(x) {
        phat <- parscf[[cfolds[x]]]$par
        
        ## Implement each dimension separately to save on computation - follows from independence of dimensions
        marg <- map(1:nkeep, function(y){
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
      Binv <- map(B, ~ solve(diag(2*nkeep) + .x))
      
      Orth <- X2sub - cmean
      orthlist[[paste0(i,j)]] <- Orth
      
      g <- as.numeric(clust$cluster[which(clust$cluster %in% c(i,j))] == i)
      dprop <- 1 - mean(g)
      
      IFpartial2 <- map(1:ntest, ~ctemp[[.x]][[2]]*g[.x]) %>% abind::abind(along=3)
      A <- map(1:nfold, function(x){apply(IFpartial2[,,which(cfolds != x)], 1:2, mean)})
      
      tautemp <- Orth * matrix(g, nrow=ntest, ncol=nkeep, byrow=F)
      debias <- map(1:ntest, function(x){2*t(A[[cfolds[x]]] %*% Binv[[cfolds[x]]] %*% cIFtheta[x,])}) %>% do.call(rbind, .)
      Dn <- colMeans(tautemp + debias)
      
      IFdebias <- map(1:ntest, function(x){2*t(A[[cfolds[x]]] %*% Binv[[cfolds[x]]] %*% (IFtheta[x,] - cIFtheta[x,]))}) %>% do.call(rbind, .)
      IFDn <- tautemp - IFdebias
      
      tauhat1 <- var(IFDn[which(g == 0),])
      tauhat2 <- var(IFDn[which(g == 1),])
      tauhat <- dprop*tauhat1 + (1-dprop)*tauhat2
      Torth <- ntest*(Dn %*% solve(tauhat) %*% Dn)
      porth <- pchisq(Torth, df=nkeep, lower.tail=F)
      
      if (is.na(porth)) {
        print("NA P-VALUE!!!")
        return()
      }
      
      pmat[i,j] <- porth
      print(paste0("Done ", i, "-", j))
    }
  }
}

pmat[lower.tri(pmat)] <- p.adjust(pmat[lower.tri(pmat)])
save(X1, X2, clust, orthlist, ctab2, pmat, file=paste0(fpath, "Plots/umi_orth.rda"))

if (sum(is.na(pmat[lower.tri(pmat)])) == 0) {
  ## Table 1
  relab <- c(3,4,1,2)
  ctab3 <- ctab2[paste0("Cluster ", relab),]
  colnames(ctab3) <- c("B-cells", "Naive Cytotoxic T-cells", "CD14 Monocytes", "Regulatory T-cells")
  rownames(ctab3) <- paste0("Cluster ", 1:4)
  ktab <- kable(ctab3, booktabs=T) %>% kable_styling(full_width=F)
  save_kable(ktab, paste0(fpath, "Plots/ctab.png"), zoom=2)
  writeLines(kable(ctab3, booktabs=T, format="latex"), paste0(fpath, "Plots/ctab.tex"))
  
  target <- c(rownames(ctab2)[which(ctab2[,"cd14.monocytes"] == max(ctab2[,"cd14.monocytes"]))],
              rownames(ctab2)[which(ctab2[,"b.cells"] == max(ctab2[,"b.cells"]))])
  target <- sort(as.numeric(gsub("Cluster ", "", target)))
  tarlab <- sort(relab[target])
  clsub <- clust$cluster[which(clust$cluster %in% target)]
  
  X2ra <- sweep(X2, 1, rowSums(X2), "/")
  
  Ysub <- orthlist[[paste(sort(target, decreasing=T), collapse="")]]
  Ysubra <- sweep(Ysub, 1, rowSums(X2)[which(clust$cluster %in% target)], "/")
  
  ## 2PC Plots for paper
  X1pca2D <- prcomp(X1ra, rank.=2)
  X2proj2D <- predict(X1pca2D, X2ra)
  Yproj2D <- predict(X1pca2D, Ysubra)
  ctype <- cm %>% 
    mutate(phenoid = case_when(
      phenoid == "b.cells" ~ "B-cell",
      phenoid == "naive.cytotoxic" ~ "Naive cytotoxic T-cell",
      phenoid == "cd14.monocytes" ~ "CD14 monocyte",
      phenoid == "regulatory.t" ~ "Regulatory T-cell",
      TRUE ~ ""
    )) %>% pull(phenoid)
  
  ## Figure 6
  p1 <- data.frame(X1pca2D$x) %>%
    ggplot(aes(x=PC1, y=PC2, colour=as.factor(relab[clust$cluster]), shape=ctype)) + 
    geom_point(alpha=0.5) +
    scale_colour_manual(values=scales::viridis_pal()(4)[relab]) +
    scale_fill_manual(values=scales::viridis_pal()(4)[relab]) +
    theme_minimal() +
    theme(legend.position="none") +
    labs(colour="Cluster")
  
  p2 <- data.frame(X2proj2D) %>%
    ggplot(aes(x=PC1, y=PC2, colour=as.factor(relab[clust$cluster]), shape=ctype)) + 
    geom_point(alpha=0.5) +
    scale_colour_manual(values=scales::viridis_pal()(4)[relab]) +
    scale_fill_manual(values=scales::viridis_pal()(4)[relab]) +
    theme_minimal() +
    theme(legend.position="bottom") +
    labs(colour="Cluster", shape="Cell subtype")
  
  p3 <- data.frame(Yproj2D) %>%
    ggplot(aes(x=PC1, y=PC2, colour=as.factor(relab[clsub]), shape=ctype[which(clust$cluster %in% target)])) + 
    geom_point(alpha=0.5) +
    scale_colour_manual(values=scales::viridis_pal()(4)[relab[tarlab]]) +
    scale_fill_manual(values=scales::viridis_pal()(4)[relab[tarlab]]) +
    theme_minimal() +
    theme(legend.position="none") +
    labs(colour="Cluster")
  pscrna <- p1 + p2 + p3
  ggsave(paste0(fpath, "Plots/scrna.png"), pscrna, width=9, height=4)
}
  


