require(tidyverse)
require(lubridate)
require(cmdstanr)
require(rstan)

#--------------------------------------------------
# Data
#--------------------------------------------------

cases <- read_csv("cases.csv") %>%
         rename(Day = ymd, Province = province) %>%
         filter(!(Province %in% c("aomen", "taiwan", "xianggang"))) %>%
         filter(Day >= ymd("2022-11-01"), Day <= ymd("2022-11-11")) %>%
         select(c(Day, Province, locAsymNum, locIncrNum)) %>%
         mutate(Cases = locAsymNum + locIncrNum) %>%
         group_by(Day) %>%
         summarize(Cases = sum(Cases)) %>%
         arrange(Day)

survey <- read_csv("survey.csv")
dec26 <- c()
dec26[1] <- survey %>%
            mutate(S = (Uninfected + Asymptomatic) * PopulationSize2021) %>%
            pull(S) %>% sum()
dec26[2] <- survey %>%
            mutate(I = Symptomatic * PopulationSize2021) %>%
            pull(I) %>% sum()
dec26[3] <- survey %>%
            mutate(R = Recovered * PopulationSize2021) %>%
            pull(R) %>% sum()

N <- sum(survey$PopulationSize2021)

data_list <- list(N = N,
                  days1 = 1:21,  # Oct 22 - Nov 11
                  days2 = 21:47, # Nov 11 - Dec 07
                  days3 = 47:66, # Dec 07 - Dec 26
                  daysC = 11:21, # Nov 01 - Nov 11
                  days4 = 66:93, # Dec 26 - Jan 22
                  cases = as.integer(cases$Cases),
                  dec26 = as.integer(dec26))
data_list[["n_days"]] <- c(length(data_list$days1), length(data_list$days2),
                           length(data_list$days3), length(data_list$daysC),
                           length(data_list$days4))

#--------------------------------------------------
# Run Stan
#--------------------------------------------------

model <- "seir"
label <- "base"

sm <- cmdstan_model(paste0(model, ".stan"))

sam <- sm$sample(data = data_list,
                 chains = 4,
                 parallel_chains = 4,
                 iter_warmup = 1000,
                 iter_sampling = 1000,
                 output_dir = "stan-cache",
                 refresh = 200
                )

sam_rstan <- read_stan_csv(sam$output_files())
chain <- rstan::extract(sam_rstan)

stan_dens(sam_rstan, pars=c("beta1", "beta2", "beta3", "f", "E0", "phi"), separate_chains=T)

save(data_list, chain, file=paste0(model, "-", label, ".rda"))
