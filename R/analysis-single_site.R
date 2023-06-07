
# Set up ####
library(PDSIR)
set.seed(1)

S0 <- 1000  # initial number of susceptible individuals
I0 <- 10    # initial number of infectious individuals

t_end <- 6 # end of observation period

iota_dist <- "weibull" # distribution of the infection periods
theta <- list(
      R0 = 2,               # basic reproduction number
      lambda = 1, shape = 2 # parameters of the Weibull distribution for the infection periods
      )


# Simulate artificial data ####
SIR <- simulate_SEM(S0, I0, t_end, theta, iota_dist)


# Observed data ####
K <- 10 # number of intervals
Y <- observed_data(SIR, K)
n_I <- sum(Y$I_k)


# DA-MCMC ####
N       <- 1e6
thin    <- 100

out_block       <- run_DAMCMC(Y, N, iota_dist = iota_dist, theta_0 = theta, rho = 0.1    , thin = thin)
out_single_site <- run_DAMCMC(Y, N, iota_dist = iota_dist, theta_0 = theta, rho = 1/(I0+n_I), thin = thin)


# Output ####
save(out_block, out_single_site, thin, file = "mcmc-draws/single_site.RDATA")
