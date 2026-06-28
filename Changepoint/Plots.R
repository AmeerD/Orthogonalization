library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(patchwork)

## Replace with path to results data frame
load("Paper/Orthresults.rda")

# Changepoint detection simulation plots

## Over Type I error (Figure S1)
cpt1e <- cpdf %>% 
  filter(Delta == 0) %>% 
  filter(setting == "conditional") %>% 
  mutate(n = factor(paste0("n=",n), levels=paste0("n=", sort(unique(cpdf$n)))),
         wvar = factor(wvar, levels=sort(unique(cpdf$wvar)))) %>% 
  ggplot(aes(sample=pval, colour=gfun, linetype=wvar)) +
  geom_abline(intercept=0, slope=1) +
  stat_qq(distribution=qunif, geom="line") +
  theme(legend.position="bottom") +
  labs(linetype=TeX(r"($Var(W)$)"), colour=TeX(r"($g(X^{(1)}_i)$)")) +
  xlab("Significance cutoff") + ylab("P-value quantiles") +
  facet_wrap(~n) +
  scale_colour_discrete(labels=c(TeX(r"($x^{(1)}_iI(i\in\hat{S}_{pre}\cup\hat{S}_{post})$)"), 
                                 TeX(r"($I(i\in\hat{S}_{post})$)"))) +
  scale_x_continuous(breaks=c(0,0.5,1), limits=c(0,1)) +
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1))
cpt1e
ggsave("Paper/Changepoint/Plots/cpt1e.png", cpt1e, width=8, height=3)

## Power plot (Figure S2)
cppower <- cpdf %>%
  filter(setting == "conditional") %>% 
  mutate(rej = pval <= 0.05) %>% 
  mutate(n = factor(paste0("n=",n), levels=paste0("n=", sort(unique(cpdf$n)))),
         wvar = factor(wvar, levels=sort(unique(cpdf$wvar)))) %>% 
  ggplot() +
  geom_smooth(aes(x=Delta, y=as.numeric(rej), colour=gfun, linetype=wvar),
              method="gam", method.args = list(family="binomial"), se=F) +
  xlab(expression(paste("Effect size (", Delta, ")", sep=""))) +
  ylab(expression("Power at"~alpha~"= 0.05")) +
  labs(linetype=TeX(r"($Var(W)$)"), colour=TeX(r"($g(X^{(1)}_i)$)")) +
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1)) +
  scale_colour_discrete(labels=c(TeX(r"($x^{(1)}_iI(i\in\hat{S}_{pre}\cup\hat{S}_{post})$)"), 
                                 TeX(r"($I(i\in\hat{S}_{post})$)"))) +
  facet_wrap(~n) +
  geom_hline(yintercept=c(0,0.05,1)) +
  theme(legend.position="bottom")
cppower
ggsave("Paper/Changepoint/Plots/cppower.png", cppower, width=8, height=3)
