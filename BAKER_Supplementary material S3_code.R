## slightly edited because I'm getting erros...

# This supplementary material is hosted by Eurosurveillance as supporting information
# alongside the article “The incubation period of 2019-nCoV infections among travellers from
# Wuhan, China”, on behalf of the authors, who remain responsible for the accuracy and
# appropriateness of the content. The same standards for ethics, copyright, attributions and
# permissions as for the article apply. Supplements are not edited by Eurosurveillance and the
# journal is not responsible for the maintenance of any links or email addresses provided
# therein.

#########################################################################################
# Incubation period estimation for 2019-nCoV
# based on 88 confirmed cases outside of Wuhan in the period 20-28 January 2020
# with known travel history from Wuhan,China and symptom onset date
# 29 January 2020
# Jantien Backer (jantien.backer@rivm.nl)
#########################################################################################

library(tidyverse)
library(rstan)
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

input.data <- list(
  N = nrow(data),
  tStartExposure = data$tStartExposure,
  tEndExposure = data$tEndExposure,
  tSymptomOnset = data$tSymptomOnset)


# compile model

model <- stan_model(model_code = "
data{
  int <lower = 1> N;
  vector[N] tStartExposure;
  vector[N] tEndExposure;
  vector[N] tSymptomOnset;
}

parameters{
  real<lower = 0> alphaInc; 	// Shape parameter of weibull distributed incubation period
  real<lower = 0> sigmaInc; 	// Scale parameter of weibull distributed incubation period
  vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
}

transformed parameters{
  vector[N] tE; 	// infection moment
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
}

model{
  // Contribution to likelihood of incubation period
  target += weibull_lpdf(tSymptomOnset -  tE  | alphaInc, sigmaInc);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = weibull_lpdf(tSymptomOnset[i] -  tE[i]  | alphaInc, sigmaInc);
  }
}
")

stanfit <- sampling(model, data = input.data, 
                init = "random",
                warmup = 4000,
                iter = 14000, 
                chains = 8)

save("stanfit", file="stanfit.rda")

# check results and convergence
print(stanfit)

# modelfit
LL = extract_log_lik(stanfit, parameter_name = 'log_lik')
loo(LL)

# results
alpha <- rstan::extract(stanfit)$alphaInc
sigma <- rstan::extract(stanfit)$sigmaInc

# posterior median and 95%CI of mean
quantile(sigma*gamma(1+1/alpha), probs = c(0.025,0.5,0.975))

# posterior median and 95%CI of sd
quantile(sqrt(sigma^2*(gamma(1+2/alpha)-(gamma(1+1/alpha))^2)), probs = c(0.025,0.5,0.975))

# posterior median and 95%CI of percentiles
percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qweibull(p = p, shape = alpha, scale = sigma), probs = c(0.025, 0.5, 0.975)))
colnames(percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
