library(tidyverse)
library(MASS)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
load("baker_cohort.rda")
load("baker_bias.rda")

data <- read_tsv(file = "BAKER_Supplementary material S1_data.tsv")

data <- data %>% 
  mutate(tReport = as.integer((`reporting date` %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tSymptomOnset = as.integer((symptom_onset %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tStartExposure = as.integer((exposure_start %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31")),
         tEndExposure = as.integer((exposure_end %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))) %>%
  # if no start of exposure (i.e. Wuhan residents) use arbitrarily chosen start exposure in far away in past (here half December 2019)
  mutate(tStartExposure = ifelse(is.na(tStartExposure), min(tSymptomOnset)-21, tStartExposure)) %>%
  # if symptom onset in Wuhan, exposure ends at symptom onset
  mutate(tEndExposure = ifelse(tEndExposure > tSymptomOnset, tSymptomOnset, tEndExposure))

tt <- table(data$tSymptomOnset)

day <- as.numeric(names(tt))
allday <- seq(min(day), max(day), by=1)

alldate <- paste0("Jan ", allday)

dd <- tibble(
  day=allday,
  cases=0
)

dd$cases[match(day, dd$day)] <- as.numeric(tt)

gfit <- glm.nb(cases~day, data=dd)

pp <- predict(gfit, se.fit=TRUE)

preddata <- data.frame(
  time=dd$day,
  est=exp(pp$fit),
  lwr=exp(pp$fit-2*pp$se.fit),
  upr=exp(pp$fit+2*pp$se.fit)
)

g1 <- ggplot(dd) +
  geom_ribbon(data=preddata, aes(time, ymin=lwr, ymax=upr), alpha=0.3) +
  geom_line(data=preddata, aes(time, est)) +
  geom_point(aes(day, cases), col="red") +
  geom_line(aes(day, cases), col="red") +
  ggtitle("A") +
  scale_x_continuous("Symptom onset dates",
                     breaks=1:6*4,
                     label=alldate[match(1:6*4, allday)]) +
  scale_y_continuous("Number of cases", expand=c(0, 0)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

baker_cohort2 <- baker_cohort %>%
  ungroup %>%
  mutate(
    type=factor(type, levels=c("cohort-based adjustment", "naive"),
                labels=c("cohort-based mean", "observed mean"))
  )

g2 <- ggplot(baker_cohort2) + 
  geom_hline(yintercept=median(expvec), col='gray', lwd=2, lty=1) +
  geom_hline(yintercept=6.4, col='gray', lwd=2, lty=2) +
  geom_point(aes(tmeasure, est, col=type, shape=type), position=position_dodge(width=0.4)) +
  geom_errorbar(aes(tmeasure, ymin=lwr, ymax=upr, col=type), width=0, position=position_dodge(width=0.4)) +
  # geom_line(aes(tmeasure, est, col=type, lty=type)) +
  # geom_ribbon(aes(tmeasure, ymin=lwr, ymax=upr, fill=type, col=type, lty=type), alpha=0) +
  ggtitle("A") +
  scale_x_continuous("Date of measurement",
                     breaks=1:6*4,
                     label=alldate[match(1:6*4, allday)]) +
  scale_y_continuous("Mean incubation period (days)", expand=c(0, 0)) +
  scale_color_manual(values=c("red", "black")) +
  scale_fill_manual(values=c("red", "black")) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    legend.title = element_blank(),
    legend.position = c(0.7, 0.13)
  )

g3 <- ggplot(baker_prop) + 
  geom_hline(yintercept=0, col='gray', lwd=2, lty=2) +
  geom_hline(yintercept=median(biasvec), col='gray', lwd=2, lty=1) +
  geom_point(aes(tmeasure, est-1)) +
  geom_errorbar(aes(tmeasure, ymin=lwr-1, ymax=upr-1), width=0) +
  ggtitle("B") +
  # geom_line(aes(tmeasure, est)) +
  # geom_ribbon(aes(tmeasure, ymin=lwr, ymax=upr), alpha=0, col="black") +
  scale_x_continuous("Date of measurement",
                     breaks=1:6*4,
                     label=alldate[match(1:6*4, allday)]) +
  scale_y_continuous("Relative bias", expand=c(0, 0),
                     limits=c(-0.5, 0.02)) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line()
  )

gtot <- arrangeGrob(g2, g3, nrow=1)

ggsave("figure_baker.pdf", gtot, width=8, height=4)
