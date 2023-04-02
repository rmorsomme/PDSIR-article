#devtools::install_github("rmorsomme/PDSIR", force = TRUE)
library(PDSIR)
library(tidyverse)
library(data.table) # rbindlist
library(coda) # effectiveSize

theme_set(theme_bw())
theme_update(
      #axis.text.x = ggplot2::element_text(size = 25),
      text = element_text(size = 25)
)


#
# Proof of concept ####

## Set up ####
set.seed(1)
S0 <- 2000  # initial number of susceptible individuals
I0 <- 10    # initial number of infectious individuals

t_end <- 6 # end of observation period

iota_dist <- "weibull" # distribution of the infection periods
theta <- list(
      R0 = 2,               # basic reproduction number
      lambda = 1, shape = 2 # parameters of the Weibull distribution for the infection periods
)


## Simulate artificial data ####
SIR <- simulate_SEM(S0, I0, t_end, theta, iota_dist)

## Trajectories of compartments ####
draw_trajectories(SIR, t_end, include_legend = FALSE, text_size = 25)

ggsave("trajectories.jpeg", path = "figures", width = 1.61803, height = 1, scale = 5)

## Observed data ####
K <- 10 # number of observation intervals
Y <- observed_data(SIR, K)
sum(Y$I_k) # total number of infections


## Run MCMC ####
theta_0 <- list(R0 = 0.2, lambda = 0.1, shape = 2) # initial values for theta
out <- run_DAMCMC(Y, N = 1e6, iota_dist = iota_dist, theta_0 = theta_0, rho = 0.1, print_i = T)

print(out$run_time)    # run time in seconds
print(out$rate_accept) # acceptance rate in the Metropolis-Hastings step for the latent data

## figures ####
THETA <- out$theta %>%
      rbindlist() %>%
      tibble() %>%
      mutate(Iteration = 1:nrow(.))

THETA %>%
      slice(1:2e4) %>%
      ggplot(aes(x = Iteration, y = lambda)) +
      geom_line() +
      geom_hline(yintercept = theta$lambda, col = "red", linewidth = 1.5, linetype = 2) +
      labs(y = expression(lambda))

ggsave("traceplot-transient.jpeg", path = "figures", width = 1.61803, height = 1, scale = 5)

THETA %>%
      slice(2e4:nrow(.)) %>%
      ggplot(aes(x = Iteration, y = lambda)) +
      geom_line() +
      geom_hline(yintercept = theta$lambda, col = "red", linewidth = 1.5, linetype = 2) +
      labs(y = expression(lambda))

ggsave("traceplot-recurrent.jpeg", path = "figures", width = 1.61803, height = 1, scale = 5)



#
# Impact of rho ####

## Set up ####
rm(list=ls())

S0s  <- c(1e2, 5e2, 1e3)
R0s  <- c(2  , 2.5, 3  )
rhos <- c(0.02, 0.05, 0.1, 0.25, 0.5, 1)

results <- tibble(
      rho = numeric(), S0 = numeric(),
      accept_rate = numeric(), run_time = numeric(),
      ESS_R0 = numeric(),ESSsec_R0 = numeric()
      )

## Run simulations ####
for(k in 1 : length(S0s)) {

      # Parameters
      S0    <- S0s[k]
      theta <- list(R0 = R0s[k], lambda = 1, shape = 2)

      # SEM
      SEM    <- simulate_SEM(S0, I0 = 10, t_end = 6, theta, iota_dist = "weibull")
      Y      <- observed_data(SEM, K = 10)

      for(rho in rhos) {  print(paste0("S0=",S0, " - rho=", rho))

            # MCMC
            out <- run_DAMCMC(Y, N = 1e5, rho, theta_0 = theta, iota_dist = "weibull")

            ESS_R0 <- out$theta %>%
                  map_dbl("R0") %>%
                  effectiveSize()

            results  <- tibble::add_row(
                  results,
                  S0            = S0,
                  rho           = rho,
                  accept_rate   = out$rate_accept,
                  run_time      = out$run_time,
                  ESS_R0        = ESS_R0,
                  ESSsec_R0     = ESS_R0 / out$run_time
                  )

      } # end-for rho

} # end-for S0

results %>%
      ggplot(aes(rho, run_time , linetype = factor(S0))) +
      geom_line() + geom_point() +
      labs(x = expression(rho), y = "Run time") +
      expand_limits(y=0) +
      theme(legend.position = "none")
ggsave("rho-time.jpeg", path = "figures", width = 1.61803, height = 1, scale = 5)

results %>%
      ggplot(aes(rho, accept_rate , linetype = factor(S0))) +
      geom_line() + geom_point() +
      labs(x = expression(rho), y = "Acceptance rate") +
      expand_limits(y=0) +
      theme(legend.position = "none")
ggsave("rho-accept.jpeg", path = "figures", width = 1.61803, height = 1, scale = 5)

results %>%
      ggplot(aes(rho, ESSsec_R0)) +
      geom_line() + geom_point() +
      labs(x = expression(rho), y = "Effective sample size per second") +
      expand_limits(y=0) +
      facet_grid(S0~., labeller = label_both, scales = "free")
ggsave("rho-ess.jpeg", path = "figures", width = 1.61803, height = 1, scale = 5)
