
## Set up ####
library(PDSIR)
set.seed(1)

S0 <- 2500  # initial number of susceptible individuals
I0 <- 10    # initial number of infectious individuals

t_end <- 6 # end of observation period

iota_dist <- "weibull" # distribution of the infection periods
theta <- list(
      R0 = 2,               # basic reproduction number
      lambda = 1, shape = 2 # parameters of the Weibull distribution for the infection periods
      )
print(complete_theta(theta, iota_dist, S0)) # value of beta


## Simulate artificial data ####
SIR <- simulate_SEM(S0, I0, t_end, theta, iota_dist)
SIR$I[SIR$t < t_end][sum(SIR$t < t_end)] # number of infectious individuals at time t=t_end



## Trajectories of compartments ####
g_trajectories <- draw_trajectories(SIR, t_end)


## Observed data ####
K <- 10 # number of intervals
Y <- observed_data(SIR, K)
print(Y$I_k) # number of infections per interval
sum(Y$I_k) # total number of infections


## DA-MCMC ####
N       <- 1e6
thin    <- 10
theta_0 <- list(R0 = 0.2, lambda = 0.1, shape = 2) # initialize theta in a low density region

out <- run_DAMCMC(Y, N, iota_dist = iota_dist, theta_0 = theta_0, rho = 0.1, thin = thin)


## Output ####
save(out, SIR, g_trajectories, thin, theta, file = "mcmc-draws/PoC.RDATA")
