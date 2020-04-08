library(dplyr)
library(readxl)
library(ggplot2); theme_set(theme_bw())
library(mvtnorm)
library(gridExtra)
library(sn)
load("serial_estimate.rda")

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

rr <- read_xlsx("Table S5.xlsx", skip=1) %>%
  mutate(
    serial=`Seconday - symptom onset date`-`Index - symptom onset date`
  )

dd <- data.frame(
  id=c(rr$`Index ID`, rr$`Secondary ID`),
  onset=c(rr$`Index - symptom onset date`, rr$`Seconday - symptom onset date`),
  province=rr$Province
) %>%
  filter(!duplicated(id))

mm <- max(rr$`Seconday - symptom onset date`)

cdata <- data_frame(
  x=-5:mm,
  y=mm-x
)

g1 <- ggplot(dd) +
  geom_bar(aes(onset), alpha=0.1, col=1) +
  scale_x_continuous("Symptom onset date", expand=c(0, 0), limits=c(-4, 30)) +
  scale_y_continuous("Daily number of symptomatic cases", expand=c(0, 0), limits=c(0, 67)) +
  ggtitle("A") +
  theme(
    panel.grid = element_blank()
  )

rr2 <- rr %>%
  group_by(`Index - symptom onset date`) %>%
  summarize(
    mean=mean(serial),
    lwr=min(serial),
    upr=max(serial)
  )

rr_summ <- rr %>%
  arrange(`Seconday - symptom onset date`) %>%
  mutate(
    csum=cumsum(serial),
    csum2=cumsum(serial^2),
    cmean=csum/1:n(),
    csd=sqrt((csum2-csum^2/1:n())/(1:n()-1)),
    cse=csd/sqrt(1:n()),
    lwr=cmean-2*cse,
    upr=cmean+2*cse
  ) %>%
  group_by(`Seconday - symptom onset date`) %>%
  summarize(
    csum=tail(csum, 1),
    cmean=tail(cmean, 1),
    lwr=tail(lwr, 1),
    upr=tail(upr, 1)
  ) %>%
  filter(`Seconday - symptom onset date` > 0) %>%
  mutate(
    upr=ifelse(upr > 6.5, 6.5, upr),
    lwr=ifelse(lwr < -4.5, -4.5, lwr)
  )
 
g2 <- ggplot(rr_summ) +
  geom_ribbon(aes(`Seconday - symptom onset date`, ymin=lwr, ymax=upr), alpha=0.2) +
  geom_line(aes(`Seconday - symptom onset date`, cmean)) +
  scale_x_continuous("Date of measurement", limits=c(-4, 30), expand=c(0, 0)) +
  scale_y_continuous("Observed mean serial interval (days)", expand=c(0, 0), limits=c(-4.5, 6.5)) +
  ggtitle("B") +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.85)
  )

fitdata <- data_frame(
  time=seq(-12, 26, by=0.01),
  density=dsn(time, xi=coef(serial_estimate_fixed)[1], omega=exp(coef(serial_estimate_fixed)[2]), alpha=coef(serial_estimate_fixed)[3])
)

g3 <- ggplot(fitdata) +
  geom_histogram(data=rr, aes(serial, y=..density..), fill=1, alpha=0.1, col=1, breaks=-30:30) +
  stat_function(fun=dnorm, args=list(mean=3.96, sd=4.75), aes(col="Du et al., 2020")) +
  geom_line(aes(time, density, col="Right-censor adjusted")) +
  scale_x_continuous("Serial interval (days)", limits=c(-12, 26), expand=c(0, 0)) +
  scale_y_continuous("Probability density", limits=c(0, 0.118), expand=c(0, 0)) +
  scale_color_manual(values=c(1, 2)) + 
  ggtitle("C") +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.75, 0.83),
    legend.title = element_blank()
  )

set.seed(101)
tdsamp <- rmvnorm(1000, mean=coef(serial_estimate_td), sigma=vcov(serial_estimate_td))

tdest <- apply(tdsamp, 1, function(x) {
  x[1] + x[2] * -4:30 + exp(x[3]) * x[4]/sqrt(1+x[4]^2)*sqrt(2/pi)
})

fitdata2 <- data_frame(
  time=-4:30,
  est=apply(tdest, 1, median),
  lwr=apply(tdest, 1, quantile, 0.025),
  upr=apply(tdest, 1, quantile, 0.975),
  type="time-dependent"
)

g4 <- ggplot(rr2) +
  geom_errorbar(aes(`Index - symptom onset date`, ymin=lwr, ymax=upr), width=0) +
  geom_point(aes(`Index - symptom onset date`, mean)) +
  geom_line(data=cdata, aes(x, y), lty=2) +
  # geom_ribbon(data=fitdata, aes(time, ymin=lwr, ymax=upr), alpha=0.3) +
  # geom_line(data=fitdata, aes(time, est)) +
  geom_ribbon(data=fitdata2, aes(time, ymin=lwr, ymax=upr), alpha=0.3) +
  geom_line(data=fitdata2, aes(time, est)) +
  scale_x_continuous("Symptom onset date of index case", expand=c(0, 0), limits=c(-4, 30)) +
  scale_y_continuous("Mean forward serial interval (days)", limits=c(-11, 21), expand=c(0, 0)) +
  scale_size_area(max_size = 4, guide=FALSE) +
  ggtitle("D") +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.8)
  )

gtot <- arrangeGrob(g1, g2, g3, g4, nrow=2)

ggsave("serial.pdf", gtot, width=8, height=6)
