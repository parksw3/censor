library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("param.R")
load("seir_sim.rda")

r <- 1/2 * (-(sigma+gamma)+sqrt((sigma-gamma)^2 + 4 * beta * sigma))

incdata <- data_frame(
  tmeasure=seir_sim$t_infectious[order(seir_sim$t_infectious)],
  tdiff=(seir_sim$t_infectious-seir_sim$t_infected)[order(seir_sim$t_infectious)],
  keep=!is.na(tdiff)
) %>%
  filter(keep) %>%
  mutate(
    csum=cumsum(tdiff),
    csum2=cumsum(tdiff^2),
    cmean=csum/1:n(),
    csd=sqrt((csum2-csum^2/1:n())/(1:n()-1)),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  )

infdata <- data_frame(
  tmeasure=seir_sim$t_recovered[order(seir_sim$t_recovered)],
  tdiff=(seir_sim$t_recovered-seir_sim$t_infectious)[order(seir_sim$t_recovered)],
  keep=!is.na(tdiff)
) %>%
  filter(keep) %>%
  mutate(
    csum=cumsum(tdiff),
    csum2=cumsum(tdiff^2),
    cmean=csum/1:n(),
    csd=sqrt((csum2-csum^2/1:n())/(1:n()-1)),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  )

gendata <- data_frame(
  tmeasure=seir_sim$t_infected[order(seir_sim$t_infected)],
  tdiff=(seir_sim$t_infected-seir_sim$t_infected[seir_sim$infected_by])[order(seir_sim$t_infected)],
  keep=!is.na(tdiff)
) %>%
  filter(keep) %>%
  mutate(
    csum=cumsum(tdiff),
    csum2=cumsum(tdiff^2),
    cmean=csum/1:n(),
    csd=sqrt((csum2-csum^2/1:n())/(1:n()-1)),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  )

serdata <- data_frame(
  tmeasure=seir_sim$t_infectious[order(seir_sim$t_infectious)],
  tdiff=(seir_sim$t_infectious-seir_sim$t_infectious[seir_sim$infected_by])[order(seir_sim$t_infectious)],
  keep=!is.na(tdiff)
) %>%
  filter(keep) %>%
  mutate(
    csum=cumsum(tdiff),
    csum2=cumsum(tdiff^2),
    cmean=csum/1:n(),
    csd=sqrt((csum2-csum^2/1:n())/(1:n()-1)),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  )

g1 <- ggplot(incdata) +
  geom_hline(yintercept=1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(sigma+r), col='gray', lty=2, lwd=2) +
  geom_line(aes(tmeasure, cmean)) +
  geom_line(aes(tmeasure, lwr)) +
  geom_line(aes(tmeasure, upr)) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_y_continuous("Mean latent period (days)", expand=c(0, 0), limits=c(0, 1/sigma+0.5)) +
  ggtitle("A") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g2 <-ggplot(infdata) +
  geom_hline(yintercept=1/gamma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r), col='gray', lwd=2, lty=2) +
  geom_line(aes(tmeasure, cmean)) +
  geom_line(aes(tmeasure, lwr)) +
  geom_line(aes(tmeasure, upr)) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_y_continuous("Mean infectious period (days)", expand=c(0, 0), limits=c(0, 1/gamma+0.5)) +
  ggtitle("B") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g3 <-ggplot(gendata) +
  geom_hline(yintercept=1/gamma+1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r) + 1/(sigma+r),, col='gray', lwd=2, lty=2) +
  geom_line(aes(tmeasure, cmean)) +
  geom_line(aes(tmeasure, lwr)) +
  geom_line(aes(tmeasure, upr)) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_y_continuous("Mean generation interval (days)", expand=c(0, 0), limits=c(0, 1/sigma + 1/gamma+1)) +
  ggtitle("C") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g4 <-ggplot(serdata) +
  geom_hline(yintercept=1/gamma+1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r) + 1/(sigma+r), col='gray', lwd=2, lty=2) +
  geom_line(aes(tmeasure, cmean)) +
  geom_line(aes(tmeasure, lwr)) +
  geom_line(aes(tmeasure, upr)) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_y_continuous("Mean serial interval (days)", expand=c(0, 0), limits=c(0, 1/sigma + 1/gamma+1)) +
  ggtitle("D") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

gtot <- arrangeGrob(g1, g2, g3, g4, nrow=2)

ggsave("figure_seir.pdf", gtot, width=10, height=6)
