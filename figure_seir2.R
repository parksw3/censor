library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(zoo)
source("param.R")
load("seir_sim.rda")
load("seir_gamma.rda")

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
  filter(keep, tmeasure < 120) %>%
  mutate(
    sw=1:n(),
    wsum=cumsum(weight),
    csum=cumsum(tdiff*weight),
    csum2=cumsum((tdiff)^2*weight),
    cmean=csum/wsum,
    cvar=(csum2-csum^2/wsum)/(wsum),
    csd=sqrt(cvar),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse,
    group=120,
    group=ifelse(tmeasure<110, 110, group),
    group=ifelse(tmeasure<100, 100, group),
    group=ifelse(tmeasure<90, 90, group),
    group=ifelse(tmeasure<80, 80, group),
    group=ifelse(tmeasure<70, 70, group),
    group=ifelse(tmeasure<60, 60, group),
    group=ifelse(tmeasure<50, 50, group),
    group=ifelse(tmeasure<40, 40, group)
  ) %>%
  group_by(group) %>%
  filter(tmeasure==max(tmeasure)) %>%
  mutate(
    tmeasure=group,
    type="growth-based mean"
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
  ) %>%
  mutate(
    group=120,
    group=ifelse(tmeasure <110, 110, group),
    group=ifelse(tmeasure <100, 100, group),
    group=ifelse(tmeasure <90, 90, group),
    group=ifelse(tmeasure <80, 80, group),
    group=ifelse(tmeasure <70, 70, group),
    group=ifelse(tmeasure <60, 60, group),
    group=ifelse(tmeasure <50, 50, group),
    group=ifelse(tmeasure <40, 40, group)
  ) %>%
  group_by(group) %>%
  filter(tmeasure==max(tmeasure)) %>%
  mutate(
    tmeasure=group,
    type="cohort-based mean"
  )

incdata2 <- data.frame(
  tmeasure=4:12*10
) %>%
  merge(incdata2, by="tmeasure", all=TRUE) %>%
  na.locf

incdata3 <- incdata3 %>%
  mutate(
    type="likelihood-based mean"
  )

infdata <- data_frame(
  tmeasure=seir_sim$t_recovered[order(seir_sim$t_recovered)],
  tdiff=(seir_sim$t_recovered-seir_sim$t_infectious)[order(seir_sim$t_recovered)],
  keep=!is.na(tdiff),
  weight=exp(r*tdiff)
) %>%
  filter(keep, tmeasure < 120) %>%
  mutate(
    sw=1:n(),
    wsum=cumsum(weight),
    csum=cumsum(tdiff*weight),
    csum2=cumsum((tdiff)^2*weight),
    cmean=csum/wsum,
    cvar=(csum2-csum^2/wsum)/(wsum),
    csd=sqrt(cvar),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse,
    group=120,
    group=ifelse(tmeasure <110, 110, group),
    group=ifelse(tmeasure <100, 100, group),
    group=ifelse(tmeasure <90, 90, group),
    group=ifelse(tmeasure <80, 80, group),
    group=ifelse(tmeasure <70, 70, group),
    group=ifelse(tmeasure <60, 60, group),
    group=ifelse(tmeasure <50, 50, group),
    group=ifelse(tmeasure <40, 40, group)
  ) %>%
  group_by(group) %>%
  filter(tmeasure==max(tmeasure)) %>%
  mutate(
    tmeasure=group,
    type="growth-based mean"
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
  ) %>%
  mutate(
    group=120,
    group=ifelse(tmeasure <110, 110, group),
    group=ifelse(tmeasure <100, 100, group),
    group=ifelse(tmeasure <90, 90, group),
    group=ifelse(tmeasure <80, 80, group),
    group=ifelse(tmeasure <70, 70, group),
    group=ifelse(tmeasure <60, 60, group),
    group=ifelse(tmeasure <50, 50, group),
    group=ifelse(tmeasure <40, 40, group)
  ) %>%
  group_by(group) %>%
  filter(tmeasure==max(tmeasure)) %>%
  mutate(
    tmeasure=group,
    type="cohort-based mean"
  )

infdata2 <- data.frame(
  tmeasure=4:12*10
) %>%
  merge(infdata2, by="tmeasure", all=TRUE) %>%
  na.locf

infdata3 <- infdata3 %>%
  mutate(
    type="likelihood-based mean"
  )

gendata <- data_frame(
  tmeasure=seir_sim$t_infected[order(seir_sim$t_infected)],
  tdiff=(seir_sim$t_infected-seir_sim$t_infected[seir_sim$infected_by])[order(seir_sim$t_infected)],
  keep=!is.na(tdiff),
  weight=exp(r*tdiff)
) %>%
  filter(keep, tmeasure < 120) %>%
  mutate(
    sw=1:n(),
    wsum=cumsum(weight),
    csum=cumsum(tdiff*weight),
    csum2=cumsum((tdiff)^2*weight),
    cmean=csum/wsum,
    cvar=(csum2-csum^2/wsum)/(wsum),
    csd=sqrt(cvar),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse,
    group=120,
    group=ifelse(tmeasure <110, 110, group),
    group=ifelse(tmeasure <100, 100, group),
    group=ifelse(tmeasure <90, 90, group),
    group=ifelse(tmeasure <80, 80, group),
    group=ifelse(tmeasure <70, 70, group),
    group=ifelse(tmeasure <60, 60, group),
    group=ifelse(tmeasure <50, 50, group),
    group=ifelse(tmeasure <40, 40, group)
  ) %>%
  group_by(group) %>%
  filter(tmeasure==max(tmeasure)) %>%
  mutate(
    tmeasure=group,
    type="growth-based mean"
  )

gendata3 <- gendata3 %>%
  mutate(
    type="likelihood-based mean"
  )

serdata <- data_frame(
  tmeasure=seir_sim$t_infectious[order(seir_sim$t_infectious)],
  tdiff=(seir_sim$t_infectious-seir_sim$t_infectious[seir_sim$infected_by])[order(seir_sim$t_infectious)],
  keep=!is.na(tdiff),
  weight=exp(r*tdiff)
) %>%
  filter(keep, tmeasure < 120) %>%
  mutate(
    sw=1:n(),
    wsum=cumsum(weight),
    csum=cumsum(tdiff*weight),
    csum2=cumsum((tdiff)^2*weight),
    cmean=csum/wsum,
    cvar=(csum2-csum^2/wsum)/(wsum),
    csd=sqrt(cvar),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse,
    group=120,
    group=ifelse(tmeasure <110, 110, group),
    group=ifelse(tmeasure <100, 100, group),
    group=ifelse(tmeasure <90, 90, group),
    group=ifelse(tmeasure <80, 80, group),
    group=ifelse(tmeasure <70, 70, group),
    group=ifelse(tmeasure <60, 60, group),
    group=ifelse(tmeasure <50, 50, group),
    group=ifelse(tmeasure <40, 40, group)
  ) %>%
  group_by(group) %>%
  filter(tmeasure==max(tmeasure)) %>%
  mutate(
    tmeasure=group,
    type="growth-based mean"
  )
  
serdata3 <- serdata3 %>%
  mutate(
    type="likelihood-based mean"
  )

incdataall <- bind_rows(incdata, incdata2, incdata3) %>%
  mutate(
    type=factor(type, levels=c("cohort-based mean", "growth-based mean", "likelihood-based mean"))
  )

g1 <- ggplot(incdataall) +
  geom_hline(yintercept=1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(sigma+r), col='gray', lty=2, lwd=2) +
  geom_errorbar(aes(tmeasure, ymin=lwr, ymax=upr, col=type), width=0, position=position_dodge(width=4)) +
  geom_point(aes(tmeasure, cmean, col=type, shape=type), position=position_dodge(width=4)) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(36, 124)) +
  scale_y_continuous("Mean latent period (days)", expand=c(0, 0), limits=c(0, NA)) +
  scale_color_manual(values=c("#D55E00", "#0072B2", "#009E73")) +
  ggtitle("A") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.title = element_blank(),
    legend.position = c(0.75,0.25)
  )

infdataall <- bind_rows(infdata, infdata2, infdata3)

g2 <- ggplot(infdataall) +
  geom_hline(yintercept=1/gamma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r), col='gray', lty=2, lwd=2) +
  geom_errorbar(aes(tmeasure, ymin=lwr, ymax=upr, col=type), width=0, position=position_dodge(width=4)) +
  geom_point(aes(tmeasure, cmean, col=type, shape=type), position=position_dodge(width=4)) +
  scale_color_manual(values=c("#D55E00", "#0072B2", "#009E73")) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(36, 124)) +
  scale_y_continuous("Mean infectious period (days)", expand=c(0, 0), limits=c(0, NA)) +
  ggtitle("B") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

gendataall <- bind_rows(gendata, gendata3)

g3 <- ggplot(gendataall) +
  geom_hline(yintercept=1/gamma+1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r) + 1/(sigma+r), col='gray', lty=2, lwd=2) +
  geom_errorbar(aes(tmeasure, ymin=lwr, ymax=upr, col=type), width=0, position=position_dodge(width=4)) +
  geom_point(aes(tmeasure, cmean, col=type, shape=type), position=position_dodge(width=4)) +
  scale_color_manual(values=c("#0072B2", "#009E73")) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(36, 124)) +
  scale_y_continuous("Mean generation interval (days)", expand=c(0, 0), limits=c(0, NA)) +
  ggtitle("C") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

serdataall <- bind_rows(serdata, serdata3)

g4 <- ggplot(serdataall) +
  geom_hline(yintercept=1/gamma+1/sigma, col='gray', lwd=2) +
  geom_hline(yintercept=1/(gamma+r) + 1/(sigma+r), col='gray', lty=2, lwd=2) +
  geom_errorbar(aes(tmeasure, ymin=lwr, ymax=upr, col=type), width=0, position=position_dodge(width=4)) +
  geom_point(aes(tmeasure, cmean, col=type, shape=type), position=position_dodge(width=4)) +
  scale_color_manual(values=c("#0072B2", "#009E73")) +
  scale_x_continuous("Time of measurement (days)", expand=c(0, 0), limits=c(36, 124)) +
  scale_y_continuous("Mean serial interval (days)", expand=c(0, 0), limits=c(0, NA)) +
  ggtitle("D") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.position = "none"
  )

gtot <- arrangeGrob(g1, g2, g3, g4, nrow=2)

ggsave("figure_seir2.pdf", gtot, width=8, height=6)
