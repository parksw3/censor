library(rstan)
library(tidyr)
library(dplyr)
load("stanfit.rda")
# options(mc.cores=1)

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

cc <- confint(gfit)

# pgamma(cc[2,2], 55, 55/coef(gfit)[2]) - pgamma(cc[2,1], 55, 55/coef(gfit)[2])

ee <- extract(stanfit)

nsample <- length(ee$alphaInc)

expvec <- meanvec <- biasvec <- rep(NA, nsample)

set.seed(101)
for (i in 1:nsample) {
  if (i %% 1000 == 0) print(i)
  alpha <- ee$alphaInc[i]
  sigma <- ee$sigmaInc[i]
  
  r <- rgamma(1, 55, 55/coef(gfit)[2])
  
  ww <- rweibull(10000, alpha, sigma)
  
  meanvec[i] <- mean(ww)
  expvec[i] <- weighted.mean(ww, exp(r*ww))
  
  biasvec[i] <- meanvec[i]/expvec[i]-1
}

save("biasvec", "meanvec", "expvec", file="baker_bias.rda")
