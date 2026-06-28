library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(patchwork)

## Two sample testing simulation plots

## Replace with path to results file
load("Paper/Orthresults.rda")

alph <- 0.05

# Overall T1E (Figure 2)
tst1e <- tsdf %>% 
  filter(d == 0) %>% 
  mutate(ns = factor(paste0("n=", ns), levels=paste0("n=", sort(unique(tsdf$ns)))),
         dim = factor(paste0("p=", dim), levels=paste0("p=", sort(unique(tsdf$dim)))),
         wvar = factor(wvar, levels=sort(unique(tsdf$wvar)))
  ) %>% 
  ggplot(aes(sample=pval, colour=gfun, linetype=wvar)) +
  stat_qq(distribution=qunif, geom="line") +
  geom_abline(intercept=0, slope=1) +
  theme(legend.position="bottom") + 
  labs(linetype=TeX(r"($Var(W)$)"), colour=TeX(r"($g(X^{(1)}_i)$)")) +
  facet_wrap(~ns, nrow=1) +
  xlab("Significance cutoff") + ylab("P-value quantiles") +
  scale_x_continuous(breaks=c(0,0.5,1), limits=c(0,1)) +
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1)) +
  scale_colour_discrete(labels=c(TeX(r"($X^{(1)}_{i1}$)"),
                                 TeX(r"($X^{(1)}_i$)"),
                                 TeX(r"($||X^{(1)}_i||_2$)")))
tst1e
ggsave("Paper/Two Sample/Plots/ts_t1e.png", tst1e, width=8, height=3)

# Power (Figure 3)
tspowdf <- tsdf %>%
  mutate(rej = pval <= 0.05,
         ns = factor(paste0("n=", ns), levels=paste0("n=", sort(unique(tsdf$ns)))),
         dim = factor(paste0("p=", dim), levels=paste0("p=", sort(unique(tsdf$dim)))),
         wvar = factor(wvar, levels=sort(unique(tsdf$wvar)))
  ) %>% 
  group_by(dim, target, d, ns, wvar, gfun) %>% 
  summarise(power=mean(rej))

tspower1 <- bind_rows(
  tspowdf,
  tspowdf %>% filter(d == 0) %>% mutate(target = "variance", d = 1),
  tspowdf %>% filter(d == 0) %>% mutate(target = "covariance")
) %>% 
  filter(ns %in% c("n=250", "n=500")) %>% 
  mutate(target = factor(target, levels=c("mean", "variance", "covariance")),
         d = ifelse(target == "variance", d-1, d)) %>% 
  ggplot(aes(x=d, y=power, colour=gfun, shape=wvar, linetype=wvar)) +
  geom_point() +
  geom_line() +
  labs(linetype=TeX(r"($Var(W)$)"), shape=TeX(r"($Var(W)$)"), colour=TeX(r"($g(X^{(1)}_i)$)")) +
  scale_colour_discrete(labels=c(TeX(r"($X^{(1)}_{i1}$)"),
                                 TeX(r"($X^{(1)}_i$)"),
                                 TeX(r"($||X^{(1)}_i||_2$)"))) +
  xlab(expression(paste("Effect size (", delta, ")", sep=""))) +
  ylab(expression("Power at"~alpha~"= 0.05")) +
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1)) +
  facet_grid(ns~target, scales="free_x") +
  geom_hline(yintercept=c(0,alph,1)) +
  theme(legend.position="bottom")
tspower2 <- bind_rows(
  tspowdf,
  tspowdf %>% filter(d == 0) %>% mutate(target = "variance", d = 1),
  tspowdf %>% filter(d == 0) %>% mutate(target = "covariance")
) %>% 
  filter(ns %in% c("n=1000", "n=2500")) %>% 
  mutate(target = factor(target, levels=c("mean", "variance", "covariance")),
         d = ifelse(target == "variance", d-1, d)) %>% 
  ggplot(aes(x=d, y=power, colour=gfun, shape=wvar, linetype=wvar)) +
  geom_point() +
  geom_line() +
  labs(linetype=TeX(r"($Var(W)$)"), shape=TeX(r"($Var(W)$)"), colour=TeX(r"($g(X^{(1)}_i)$)")) +
  scale_colour_discrete(labels=c(TeX(r"($X^{(1)}_{i1}$)"),
                                 TeX(r"($X^{(1)}_i$)"),
                                 TeX(r"($||X^{(1)}_i||_2$)"))) +
  xlab(expression(paste("Effect size (", delta, ")", sep=""))) +
  ylab(expression("Power at"~alpha~"= 0.05")) +
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(0,1)) +
  facet_grid(ns~target, scales="free_x") +
  geom_hline(yintercept=c(0,alph,1)) +
  guides(colour="none", shape="none", linetype="none")
tspower <- (tspower1 + tspower2) + plot_layout(axis_titles = "collect", guides="collect") &
  theme(legend.position="bottom")
tspower
ggsave("Paper/Two Sample/Plots/ts_pow.png", tspower, width=8, height=4)

