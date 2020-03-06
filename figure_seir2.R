library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("param.R")
load("seir_sim.rda")

r <- 1/2 * (-(sigma+gamma)+sqrt((sigma-gamma)^2 + 4 * beta * sigma))

stepfun <- function(x) {
  nn <- 1
  ww <- 1
  mm <- x[1]
  
  while (nn <= length(x)) {
    if (mm >= x[nn]) {
      x[nn] <- mm
    } else {
      ww <- nn
      mm <- x[ww]
    }
    
    nn <- nn + 1
  }
    
  x
}

incdata <- data_frame(
  tmeasure=seir_sim$t_infectious[order(seir_sim$t_infectious)],
  tdiff=(seir_sim$t_infectious-seir_sim$t_infected)[order(seir_sim$t_infectious)],
  keep=!is.na(tdiff),
  weight=exp(r*tdiff)
) %>%
  filter(keep) %>%
  mutate(
    sw=1:n(),
    wsum=cumsum(weight),
    csum=cumsum(tdiff*weight),
    csum2=cumsum((tdiff)^2*weight),
    cmean=csum/wsum,
    cvar=(csum2-csum^2/wsum)/(wsum-1),
    csd=sqrt(cvar),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  )

incdata2 <- data_frame(
  tmeasure=seir_sim$t_infectious[order(seir_sim$t_infected)],
  tdiff=(seir_sim$t_infectious-seir_sim$t_infected)[order(seir_sim$t_infected)],
  keep=!is.na(tdiff)
) %>%
  filter(keep) %>%
  mutate(
    tmeasure=stepfun(tmeasure),
    csum=cumsum(tdiff),
    csum2=cumsum(tdiff^2),
    cmean=csum/1:n(),
    csd=sqrt((csum2-csum^2/1:n())/(1:n()-1)),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  ) %>%
  group_by(tmeasure) %>%
  summarize(
    cmean=tail(cmean, 1),
    lwr=tail(lwr, 1),
    upr=tail(upr, 1)
  )

infdata <- data_frame(
  tmeasure=seir_sim$t_recovered[order(seir_sim$t_recovered)],
  tdiff=(seir_sim$t_recovered-seir_sim$t_infectious)[order(seir_sim$t_recovered)],
  keep=!is.na(tdiff),
  weight=exp(r*tdiff)
) %>%
  filter(keep) %>%
  mutate(
    sw=1:n(),
    wsum=cumsum(weight),
    csum=cumsum(tdiff*weight),
    csum2=cumsum((tdiff)^2*weight),
    cmean=csum/wsum,
    cvar=(csum2-csum^2/wsum)/(wsum-1),
    csd=sqrt(cvar),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  )

infdata2 <- data_frame(
  tmeasure=seir_sim$t_recovered[order(seir_sim$t_infectious)],
  tdiff=(seir_sim$t_recovered-seir_sim$t_infectious)[order(seir_sim$t_infectious)],
  keep=!is.na(tdiff)
) %>%
  filter(keep) %>%
  mutate(
    tmeasure=stepfun(tmeasure),
    csum=cumsum(tdiff),
    csum2=cumsum(tdiff^2),
    cmean=csum/1:n(),
    csd=sqrt((csum2-csum^2/1:n())/(1:n()-1)),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  ) %>%
  group_by(tmeasure) %>%
  summarize(
    cmean=tail(cmean, 1),
    lwr=tail(lwr, 1),
    upr=tail(upr, 1)
  )

gendata <- data_frame(
  tmeasure=seir_sim$t_infected[order(seir_sim$t_infected)],
  tdiff=(seir_sim$t_infected-seir_sim$t_infected[seir_sim$infected_by])[order(seir_sim$t_infected)],
  keep=!is.na(tdiff),
  weight=exp(r*tdiff)
) %>%
  filter(keep) %>%
  mutate(
    sw=1:n(),
    wsum=cumsum(weight),
    csum=cumsum(tdiff*weight),
    csum2=cumsum((tdiff)^2*weight),
    cmean=csum/wsum,
    cvar=(csum2-csum^2/wsum)/(wsum-1),
    csd=sqrt(cvar),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  )

serdata <- data_frame(
  tmeasure=seir_sim$t_infectious[order(seir_sim$t_infectious)],
  tdiff=(seir_sim$t_infectious-seir_sim$t_infectious[seir_sim$infected_by])[order(seir_sim$t_infectious)],
  keep=!is.na(tdiff),
  weight=exp(r*tdiff)
) %>%
  filter(keep) %>%
  mutate(
    sw=1:n(),
    wsum=cumsum(weight),
    csum=cumsum(tdiff*weight),
    csum2=cumsum((tdiff)^2*weight),
    cmean=csum/wsum,
    cvar=(csum2-csum^2/wsum)/(wsum-1),
    csd=sqrt(cvar),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  )

g1 <- ggplot(incdata2) +
  geom_hline(yintercept=1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(sigma+r), col='gray', lty=2, lwd=2) +
  geom_line(data=incdata, aes(tmeasure, cmean, col="exponential adjustment")) +
  geom_line(data=incdata, aes(tmeasure, lwr, col="exponential adjustment")) +
  geom_line(data=incdata, aes(tmeasure, upr, col="exponential adjustment")) +
  geom_line(aes(tmeasure, cmean, col="cohort-based adjustment")) +
  geom_line(aes(tmeasure, lwr, col="cohort-based adjustment")) +
  geom_line(aes(tmeasure, upr, col="cohort-based adjustment")) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_y_continuous("Mean latent period (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_color_manual(values=c(2, 1)) +
  ggtitle("A") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.title = element_blank(),
    legend.position = c(0.25,0.9)
  )

g2 <-ggplot(infdata2) +
  geom_hline(yintercept=1/gamma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r), col='gray', lty=2, lwd=2) +
  geom_line(data=infdata, aes(tmeasure, cmean)) +
  geom_line(data=infdata, aes(tmeasure, lwr)) +
  geom_line(data=infdata, aes(tmeasure, upr)) +
  geom_line(aes(tmeasure, cmean), col="red") +
  geom_line(aes(tmeasure, lwr), col="red") +
  geom_line(aes(tmeasure, upr), col="red") +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_y_continuous("Mean infectious period (days)", expand=c(0, 0), limits=c(0, NA)) +
  ggtitle("B") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g3 <-ggplot(gendata) +
  geom_hline(yintercept=1/gamma+1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r) + 1/(sigma+r), col='gray', lty=2, lwd=2) +
  geom_line(aes(tmeasure, cmean)) +
  geom_line(aes(tmeasure, lwr)) +
  geom_line(aes(tmeasure, upr)) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_y_continuous("Mean generation interval (days)", expand=c(0, 0), limits=c(0, NA)) +
  ggtitle("C") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

g4 <-ggplot(serdata) +
  geom_hline(yintercept=1/gamma+1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r) + 1/(sigma+r), col='gray', lty=2, lwd=2) +
  geom_line(aes(tmeasure, cmean)) +
  geom_line(aes(tmeasure, lwr)) +
  geom_line(aes(tmeasure, upr)) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_y_continuous("Mean serial interval (days)", expand=c(0, 0), limits=c(0, NA)) +
  ggtitle("D") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

gtot <- arrangeGrob(g1, g2, g3, g4, nrow=2)

ggsave("figure_seir2.pdf", gtot, width=10, height=6)
