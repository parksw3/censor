library(rstan)
library(tidyr)
library(dplyr)
load("stanfit.rda")

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

ee <- extract(stanfit)

nsamp <- length(ee$alphaInc)

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

reslist <- vector('list', length(nsamp))

for (i in 1:nsamp) {
  if (i %% 1000 == 0) print(i)
  
  infected <- ee$tE[i,]
  symptomatic <- data$tSymptomOnset
  
  d1 <- data.frame(
    tmeasure=symptomatic[order(symptomatic)],
    tdiff=(symptomatic-infected)[order(symptomatic)]
  ) %>%
    mutate(
      csum=cumsum(tdiff),
      cmean=csum/1:n()
    ) %>%
    group_by(tmeasure) %>%
    summarize(cmean=tail(cmean,1)) %>%
    mutate(
      type="naive"
    )
  
  d2 <- data.frame(
    tmeasure=symptomatic[order(infected)],
    tdiff=(symptomatic-infected)[order(infected)]
  ) %>%
    mutate(
      tmeasure=stepfun(tmeasure),
      csum=cumsum(tdiff),
      cmean=csum/1:n()
    ) %>%
    group_by(tmeasure) %>%
    summarize(cmean=tail(cmean, 1))  %>%
    mutate(
      type="cohort-based adjustment"
    )
  
  reslist[[i]] <- bind_rows(d1, d2)
}

baker_cohort <- reslist %>%
  bind_rows(.id="sim") %>%
  group_by(type, tmeasure) %>%
  summarize(
    est=median(cmean),
    lwr=quantile(cmean, 0.025),
    upr=quantile(cmean, 0.975)
  )

baker_prop <- reslist %>%
  bind_rows(.id="sim") %>%
  group_by(sim, tmeasure) %>%
  summarize(
    prop=ifelse(length(cmean)==1, NA, cmean[type=="naive"]/cmean[type=="cohort-based adjustment"])
  ) %>%
  filter(!is.na(prop)) %>%
  group_by(tmeasure) %>%
  summarize(
    est=median(prop),
    lwr=quantile(prop, 0.025),
    upr=quantile(prop, 0.975)
  )

save("baker_cohort", "baker_prop", file="baker_cohort.rda")
