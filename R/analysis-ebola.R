library(PDSIR)
set.seed(1)


# set up ####
S0 <- 292e3 # https://en.wikipedia.org/wiki/Prefectures_of_Guinea
I0 <- 5

ebola_gueckedou$S0    <- S0
ebola_gueckedou$I0    <- I0
ebola_gueckedou$t_end <- max(ebola_gueckedou$ts)

iota_dist <- "weibull"
theta_0 <- list(R0 = 1, lambda = 0.01, shape = 2)


# run MCMC ####
N    <- 1e6
thin <- 1e1
out  <- run_DAMCMC(ebola_gueckedou, N = N, rho = 0.1, iota_dist = iota_dist, thin = thin, theta_0 = theta_0)

save(out, thin, file = "mcmc-draws/ebola.RDATA")
