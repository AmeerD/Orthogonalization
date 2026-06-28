library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(patchwork)

## Clustering simulation plots

## Replace with path to results file
load("Paper/Orthresults.rda")

alph <- 0.05
nps <- sort(unique(clustdf$n/clustdf$p))

# Overall T1E (Figure 4a)
clustt1e <- clustdf %>% 
  filter(Delta == 0) %>% 
  filter(i == 1) %>% 
  filter(method == "conditional") %>% 
  mutate(n = factor(paste0("n=",n/p,"p"), levels=paste0("n=", nps, "p")),
         p = factor(paste0("p=",p), levels=paste0("p=", sort(unique(clustdf$p))))) %>% 
  ggplot(aes(sample=pval, colour=gfun)) +
  geom_abline(intercept=0, slope=1) +
  stat_qq(distribution=qunif, geom="line") +
  theme(legend.position="bottom") +
  labs(colour=TeX(r"($g(X^{(1)}_i)$)")) +
  xlab("Significance cutoff") + ylab("P-value quantiles") +
  facet_grid(n~p) +
  scale_colour_discrete(labels=c(TeX(r"($I(i\in\hat{C}_1)$)"),
                                 TeX(r"($||X^{(1)}_i||_2 I(i\in\hat{C}_1\cup\hat{C}_2)$)"),
                                 TeX(r"($||X^{(1)}_i||_\infty I(i\in\hat{C}_1\cup\hat{C}_2)$)"))) +
  scale_x_continuous(breaks=c(0,0.5,1), limits=c(0,1)) +
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1)) 
clustt1e
ggsave("Paper/Clustering/Plots/clust_t1e.png", clustt1e, width=6, height=7)

# Conditional T1E (Figure 5)
ct1e <- clustdf %>%
  filter(Delta == 0, p == 2) %>% 
  filter(method == "conditional") %>% 
  group_by(n, p, wvar, gfun, job) %>%
  summarise(t1e = mean(pval <= alph, na.rm=T)) %>%
  mutate(n = factor(paste0("n=",n), levels=paste0("n=", sort(unique(clustdf$n)))),
         p = factor(paste0("p=",p), levels=paste0("p=", sort(unique(clustdf$p))))) %>% 
  ggplot(aes(x=job, y=t1e, colour=gfun)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~n) +
  geom_hline(yintercept=alph) +
  labs(colour=TeX(r"($g(X^{(1)}_i)$)")) +
  scale_colour_discrete(labels=c(TeX(r"($I(i\in\hat{C}_1)$)"),
                                 TeX(r"($||X^{(1)}_i||_2I(i\in\hat{C}_1\cup\hat{C}_2)$)"),
                                 TeX(r"($||X^{(1)}_i||_\infty I(i\in\hat{C}_1\cup\hat{C}_2)$)"))) +
  xlab("Sample ") + ylab("Conditional type I error rate") +
  theme(legend.position="bottom")
ct1e
ggsave("Paper/Clustering/Plots/clust_ct1e.png", ct1e, width=8, height=3)

# Power (Figure 4b)
cpower <- clustdf %>%
  filter(method == "conditional") %>% 
  mutate(rej = pval <= 0.05,
         n = factor(paste0("n=",n/p,"p"), levels=paste0("n=", nps, "p")),
         p = factor(paste0("p=",p), levels=paste0("p=", sort(unique(clustdf$p))))) %>% 
  ggplot() +
  geom_smooth(aes(x=Delta, y=as.numeric(rej), colour=gfun),
              method="gam", method.args = list(family="binomial"), se=F) +
  labs(colour=TeX(r"($g(X^{(1)}_i)$)")) +
  scale_colour_discrete(labels=c(TeX(r"($I(i\in\hat{C}_1)$)"),
                                 TeX(r"($||X^{(1)}_i||_2 I(i\in\hat{C}_1\cup\hat{C}_2)$)"),
                                 TeX(r"($||X^{(1)}_i||_\infty I(i\in\hat{C}_1\cup\hat{C}_2)$)"))) +
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) +
  ylab(expression("Power at"~alpha~"= 0.05")) +
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1)) +
  facet_grid(n~p, scales="free_x") +
  geom_hline(yintercept=c(0,alph,1)) +
  theme(legend.position="bottom")
cpower
ggsave("Paper/Clustering/Plots/clust_pow.png", cpower, width=6, height=7)

