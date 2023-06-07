library(PDSIR)
library(data.table) # rbindlist
library(tidyverse)

a <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(is.na(a)) a <- 0

set.seed(a)

## Set up ####
S0 <- 1000  # initial number of susceptible individuals
I0 <- 10    # initial number of infectious individuals

t_end <- 6 # end of observation period

iota_dist <- "weibull" # distribution of the infection periods
theta   <- list(R0 = 2  , lambda = 1  , shape = 2)
theta_0 <- list(R0 = 0.2, lambda = 0.1, shape = 2)

N      <- 1e6
thin   <- 1e1
burnin <- min(1e5, N / 2)

## Simulate artificial data ####
SIR <- simulate_SEM(S0, I0, t_end, theta, iota_dist)

## Observed data ####
Y <- observed_data(SIR, K = 10)

## Run MCMC ####
out <- run_DAMCMC(Y, N, iota_dist = iota_dist, theta_0 = theta_0, rho = 0.2, thin = thin)

## Coverage ####
theta_true <- complete_theta(theta, iota_dist = iota_dist, S0=S0)

THETA_summary <- out$theta %>%
      rbindlist() %>%
      tibble() %>%
      select(-shape) %>%
      mutate(Iteration = (1:nrow(.))*thin) %>%
      filter(Iteration > burnin) %>% # remove burnin
      pivot_longer(R0:beta, names_to = "var", values_to = "draws") %>%
      group_by(var) %>%
      summarize(
            q_low  = quantile(draws, probs = 0.05),
            q_high = quantile(draws, probs = 0.95),
            posterior_mean = mean(draws)
            ) %>%
      mutate(
            a = a, N = N,
            rate_accept = out$rate_accept,
            theta_true = c(
                  theta_true$R0, theta_true$beta, theta_true$gamma, theta_true$lambda
            ),
            coverage = between(theta_true, q_low, q_high)
      )

# output
save(THETA_summary, file = paste0("figures/coverage/a=", a, ".RDATA"))
