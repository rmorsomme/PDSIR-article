library(tidyverse)
library(PDSIR)
library(coda)

a <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## Set up ####
S0s  <- c(1e2, 5e2, 1e3)
R0s  <- c(1.25, 1.5, 2 )
rhos <- c(0.02, 0.05, 0.1, 0.25, 0.5, 1)

results <- tibble(
      rho = numeric(), S0 = numeric(),
      accept_rate = numeric(), run_time = numeric(),
      ESS_R0 = numeric(),ESSsec_R0 = numeric(),
      ESS_beta = numeric(),ESSsec_beta = numeric()
)

# Parameters
S0    <- S0s[a]
theta <- list(R0 = R0s[a], lambda = 1, shape = 2)

# SEM
set.seed(1)
SEM    <- simulate_SEM(S0, I0 = 10, t_end = 6, theta, iota_dist = "weibull")
Y      <- observed_data(SEM, K = 10)

for(rho in rhos) {  print(paste0("S0=",S0, " - rho=", rho))

      # MCMC
      out <- run_DAMCMC(Y, N = 1e6, rho, theta_0 = theta, iota_dist = "weibull", thin = 1e1)

      ESS_R0 <- out$theta %>%
            map_dbl("R0") %>%
            effectiveSize()
      ESS_beta <- out$theta %>%
            map_dbl("beta") %>%
            effectiveSize()

      results  <- tibble::add_row(
            results,
            S0            = S0,
            rho           = rho,
            accept_rate   = out$rate_accept,
            run_time      = out$run_time,
            ESS_R0        = ESS_R0,
            ESSsec_R0     = ESS_R0 / out$run_time,
            ESS_beta      = ESS_beta,
            ESSsec_beta   = ESS_beta / out$run_time
      )

} # end-for rho

save(results, file = paste0("mcmc-draws/rho-S0=",S0, ".RDATA"))
