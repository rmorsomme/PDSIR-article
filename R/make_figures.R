
#
# setup ####

library(default) # default()
library(data.table) # rbindlist()
library(coda) # effectiveSize()
library(tidyverse)
library(PDSIR)

rm(list=ls())

# default configuration for figures
theme_set(theme_bw())
theme_update(text = element_text(size = 25))


save_figures <- function(name, path = "figures", scale){

   # golden ratio
   phi = (sqrt(5)+1)/2

   # save .jpg (for presentation)
   ggsave(
      paste0(name, ".jpg"), path = path,
      scale = scale, width = phi, height = 1, units = "cm"
   )

   # save .eps (for publication)
   ggsave(
      paste0(name, ".eps"), path = path,
      device = "eps", dpi = "retina",
      scale = scale, width = phi, height = 1, units = "cm"
   )
}


#
# proof of concept ####

load("mcmc-draws/PoC.RDATA")

print(out$run_time)    # run time in seconds
print(out$rate_accept) # acceptance rate in the Metropolis-Hastings step for the latent data

g_trajectories + labs(x="Time", y="Size of\ncompartments") + theme(legend.position="none",text = element_text(size = 25))
save_figures("trajectories", scale = 7.5)



THETA <- out$theta %>%
   rbindlist() %>%
   tibble() %>%
   mutate(iteration = (1:nrow(.))*thin)

burnin <- 25e3
THETA %>%
   filter(iteration < burnin) %>%
   mutate(iteration_1000 = iteration / 1e3) %>%
   ggplot(aes(x = iteration_1000, y = lambda)) +
   geom_line() +
   geom_hline(yintercept = theta$lambda, col = "red", linewidth = 1.5, linetype = 2) +
   labs(x = "Iteration (in 1,000)", y = expression(lambda)) +
   scale_x_continuous(breaks=c(0,10,20))
save_figures("traceplot_transient", scale = 7.5)

THETA %>%
   filter(iteration > burnin) %>%
   mutate(iteration_1000 = iteration / 1e3) %>%
   ggplot(aes(x = iteration_1000, y = lambda)) +
   geom_line() +
   geom_hline(yintercept = theta$lambda, col = "red", linewidth = 1.5, linetype = 2) +
   labs(x = "Iteration (in 1,000)", y = expression(lambda)) +
   scale_x_continuous(breaks=c(0,500,1e3))
save_figures("traceplot_recurrent", scale = 7.5)


#
# coverage ####

load("figures/coverage/a=1.RDATA")
THETA_summary_all <- THETA_summary

for(a in 2:1e3){
   load(paste0("figures/coverage/a=",a,".RDATA"))
   THETA_summary_all <- rbind(THETA_summary_all, THETA_summary)
}

THETA_summary_all %>%
   group_by(var) %>%
   summarize(
      coverage   = mean(coverage),
      average    = mean(posterior_mean),
      sd         = sd(posterior_mean),
      theta_true = unique(theta_true)
   )

summary(THETA_summary_all$rate_accept)


#
# rho ####

load("mcmc-draws/rho-S0=100.RDATA")
results_all <- results

for(S0 in c(500, 1000)){
   load(paste0("mcmc-draws/rho-S0=",S0,".RDATA"))
   results_all <- rbind(results_all, results)
}

results_all %>%
   ggplot(aes(rho, run_time , linetype = factor(S0))) +
   geom_line() + geom_point() +
   labs(x = expression(rho), y = "Run time (seconds)") +
   expand_limits(y=0) +
   theme(legend.position = "none") +
   theme(axis.title.y = ggplot2::element_text(size = 23))
save_figures("rho_time", scale = 10)

results_all %>%
   ggplot(aes(rho, accept_rate , linetype = factor(S0))) +
   geom_line() + geom_point() +
   labs(x = expression(rho), y = "Acceptance rate") +
   expand_limits(y=0) +
   theme(legend.position = "none")
save_figures("rho_accept", scale = 10)

results_all %>%
   ggplot(aes(rho, log(ESSsec_beta), linetype = factor(S0))) +
   geom_line() + geom_point() +
   labs(x = expression(rho), y = "Effective sample size\n per second (log)") +
   expand_limits(y=0) + #+
#   facet_grid(S0~., labeller = label_both, scales = "free")
   theme(legend.position = "none") +
   theme(axis.title.y = ggplot2::element_text(size = 21))
save_figures("rho_ess", scale = 10)

#
# Ebola ####
load("mcmc-draws/ebola.RDATA")

sum(ebola_gueckedou$I_k) # total number of observed infections

tibble(I_k = ebola_gueckedou$I_k, time = ebola_gueckedou$dates) %>%
   ggplot(aes(time, I_k)) +
   geom_col() +
   labs(x = "Time (weeks)", y = "Count of\ninfections")
save_figures("ebola_Ik", scale = 10)


print(out$run_time / 3600) # hours
print(out$rate_accept)

THETA <- out$theta %>%
   data.table::rbindlist() %>%
   tibble() %>%
   mutate(
      iteration = (1:nrow(.) * thin),
      iteration_1000 = iteration / 1000,
      EIP = lambda^{-1/shape}*gamma(1+1/shape) # expectation of weibull distribution
   ) %>%
   filter(iteration > 100e3) # discard burnin

print(effectiveSize(THETA))

THETA %>%
   ggplot(aes(x = EIP)) +
   geom_histogram(aes(y = after_stat(density))) +
   labs(x = "Expected duration\nof infection (days)", y = "Density")
save_figures("ebola_histogram", scale = 10)

THETA %>%
   ggplot(aes(x = iteration_1000, y = EIP)) +
   geom_line() +
   labs(x = "Iteration (in 1,000)", y = "Expected duration\nof infection (days)")
save_figures("ebola_traceplot", scale = 10)


#
# single-site versus block sampler ####

load("mcmc-draws/single_site.RDATA")

print(out_single_site$rate_accept)
print(out_block$rate_accept      )

print(out_single_site$run_time / 60) # minutes
print(out_block$run_time       / 60) # minutes

THETA_single_site <- out_single_site$theta %>%
   rbindlist() %>%
   tibble() %>%
   mutate(iteration = (1:nrow(.))*thin)
THETA_block <- out_block$theta %>%
   rbindlist() %>%
   tibble() %>%
   mutate(iteration = (1:nrow(.))*thin)

## traceplots
THETA_single_site %>%
   mutate(iteration_1000 = iteration / 1e3) %>%
   ggplot(aes(x = iteration_1000, y = beta)) +
   geom_line() +
   labs(x = "Iteration (in 1,000)", y = expression(beta)) +
   scale_x_continuous(breaks=c(0,500,1000))
save_figures("single_site_traceplot", scale = 8)

THETA_block %>%
   mutate(iteration_1000 = iteration / 1e3) %>%
   ggplot(aes(x = iteration_1000, y = beta)) +
   geom_line() +
   labs(x = "Iteration (in 1,000)", y = expression(beta)) +
   scale_x_continuous(breaks=c(0,500,1000))
save_figures("block_traceplot", scale =  8)

## acf
n_lag <- 50
tibble(
   acf = acf(THETA_block$beta, lag.max = n_lag)$acf %>% as.vector(),
   lag = 0 : n_lag
   ) %>%
   ggplot(aes(lag, acf)) +
   geom_col(width = 0.25) +
   labs(x = "Lags", y = "Auto-correlation\nfunction")
save_figures("block_acf", scale = 8)

tibble(
   acf = acf(THETA_single_site$beta, lag.max = n_lag)$acf %>% as.vector(),
   lag = 0:n_lag
) %>%
   ggplot(aes(lag, acf)) +
   geom_col(width = 0.25) +
   labs(x = "Lags", y = "Auto-correlation\nfunction")
save_figures("single_site_acf", scale = 8)

## ESS/sec
THETA_single_site %>%
   pivot_longer(cols=R0:beta, names_to = "var", values_to = "draws") %>%
   group_by(var) %>%
   summarize(ESS_sec = effectiveSize(draws) / out_single_site$run_time)

THETA_block %>%
   pivot_longer(cols=R0:beta, names_to = "var", values_to = "draws") %>%
   group_by(var) %>%
   summarize(ESS_sec = effectiveSize(draws) / out_block$run_time)


#
# misc. ####


## mu ####

set.seed(1)

S0    <- 10
I0    <- 2
t_end <- 0.6
iota_dist <- "weibull"
theta <- list(R0 = 2, lambda = 5, shape = 1) %>% complete_theta(iota_dist, S0)
K     <- 3

SIR   <- simulate_SEM(S0, I0, t_end, theta, iota_dist)

tmp <- SIR$x$tau_I
df_infections <- tibble(Time = tmp[ 0 < tmp & tmp < t_end])
tmp <- SIR$x$tau_R
df_removals   <- tibble(Time = tmp[ 0 < tmp & tmp < t_end])

df_SIR <- tibble(
   Time = SIR$t,
   mu   = theta$beta * SIR$I
   ) %>%
   filter(Time < t_end) %>%
   add_row(Time = t_end, mu = 4) # add last row manually so that the last step do not stop with vertical line.

df_PD <- tibble(
   Time  = seq(0, t_end, length.out = K + 1),
   Time_SIR = Time %>% map_dbl(~ max(df_SIR$Time[df_SIR$Time <= .]))
   ) %>%
   left_join(df_SIR, by = c("Time_SIR" = "Time"))

df_PD[K + 1, 3] <- df_PD[K, 3] # change the last value of mu_PD so that the last step do not stop with vertical line.

df_SIR %>%
   ggplot(aes(x = Time, y = mu), alpha = 0.5, size = 1) +
   geom_step(linetype = 2) +
   geom_step(data = df_PD) +
   geom_point(data = df_infections, y = -0.075, size = 3, shape = 17) +
   geom_point(data = df_removals  , y =  0.075, size = 3, shape = 15) +
   expand_limits(x = 0, y = 0) +
   labs(x = "Time", y = expression("Infection rate "*mu))
save_figures("mu", scale = 15)


#
## SIR versus PDSIR ####

S0    <- 1e4
I0    <- 10
theta <- list(R0 = 2.5, lambda = 1, shape = 2)
t_end <- 6
iota_dist <- "weibull"
SIR   <- simulate_SEM(S0, I0, t_end, theta, iota_dist)

df_SIR <- tibble(
   Time = SIR$t,
   S    = SIR$S,
   I    = SIR$I,
   R    = (S0+I0) - SIR$S - SIR$I
) %>%
   pivot_longer(cols=S:R, names_to="compartment", values_to = "count") %>%
   filter(Time < t_end) %>%
   mutate(compartment = factor(compartment, levels = c("S", "I", "R")))

for(k in c(5, 10, 50, 1e3)){

   Y <- observed_data(SIR, k)

   PDSIR <- rprop_x(theta, Y, iota_dist = "weibull", gener = F, approx = "ldp")
   SS <- suff_stat(PDSIR, Y, return_SI = TRUE, gener=F, b=1)

   df_PDSIR <- tibble(
      Time = SS$t,
      S    = SS$S,
      I    = SS$I,
      R    = (S0+I0) - SS$S - SS$I
   ) %>%
      pivot_longer(cols=S:R, names_to="compartment", values_to = "count") %>%
      filter(Time < t_end) %>%
      mutate(compartment = factor(compartment, levels = c("S", "I", "R")))

   df_SIR %>%
      ggplot(aes(Time, count, group = compartment, color = compartment), size = 1.25) +
      geom_line(alpha = 0.25, color = "black") +
      geom_line(data = df_PDSIR) +
      theme(legend.position = "none") +
      labs(y = "Size of\ncompartment")

   save_figures(paste0("comparison_k=", k), scale = 10)

}


#
## gamma minorization ####

# Setup
x <- seq(0, 8, by = 0.01)

a <- 2
b <- 0.5
A <- 1
B <- 1

xa <- a/B*log(1+B/b)
xb <- 1/b*(gamma(a+A)/gamma(a))^(1/A)
xA <- (a+A)/B*log(1+B/b)
xB <- 1/(b+B)*(gamma(a+A)/gamma(a))^(1/A)


# Gamma 1 (alpha in (0,A), beta=0)
df1 <- expand.grid(alpha = seq(0, A, by = 0.25), beta = 0, x = x) %>%
   mutate(
      Density = list(x, alpha, beta) %>% pmap(~ dgamma(..1, a + ..2, b + ..3))
   ) %>%
   unnest(Density)

ggplot(df1, aes(x = x, y = Density)) +
   geom_line(aes(col = as.factor(alpha)), linewidth = 1) +
   scale_color_discrete(name = expression(Parameter~alpha)) +
   geom_vline(xintercept = xb) +
   annotate(geom = "text", x = xb + 0.5, y = 0.18, label = expression(x[b]), size = 6)
save_figures("gamma1a", scale = 10)


# Gamma 2 (alpha=0, beta in (0,B))
df2 <- expand.grid(alpha = 0, beta = seq(0, B, by = 0.25), x = x) %>%
   mutate(
      Density = list(x, alpha, beta) %>% purrr::pmap(~ dgamma(..1, a + ..2, b + ..3))
      ) %>%
   unnest(Density)

ggplot(df2, aes(x = x, y = Density)) +
   geom_line(aes(col = as.factor(beta)), linewidth = 1) +
   scale_color_discrete(name = expression(Parameter~beta)) +
   geom_vline(xintercept = xa) +
   annotate(geom = "text", x = xa + 0.5, y = 0.525, label = expression(x[a]), size = 6)
save_figures("gamma1b", scale = 10)


# Gamma joint (alpha in {0,1}, beta in {0,1})
df3 <- expand.grid(alpha = c(0,A), beta = c(0,B), x = x) %>%
   mutate(
      Density = list(x, alpha, beta) %>% pmap(~ dgamma(..1, a + ..2, b + ..3))) %>%
   unnest(Density) %>%
   mutate(code = paste0(alpha, beta))

ggplot(df3, aes(x = x, y = Density)) +
   geom_line(aes(col = code), linewidth = 1) +
   scale_color_discrete(
      name = substitute(paste("(", alpha, ", ", beta, ")")),
      breaks = c("00", "01", "10", "11"),
      labels = c("(0, 0)", "(0, 1)", "(1, 0)", "(1, 1)")
   ) +
   geom_vline(xintercept = c(xa, xA, xb, xB)) +
   annotate(geom = "text", x = xa + 0.35, y = 0.525, label = expression(x[a]), size = 6) +
   annotate(geom = "text", x = xb + 0.35, y = 0.525, label = expression(x[b]), size = 6) +
   annotate(geom = "text", x = xA + 0.35, y = 0.525, label = expression(x[A]), size = 6) +
   annotate(geom = "text", x = xB + 0.35, y = 0.525, label = expression(x[B]), size = 6)
save_figures("gamma2", scale = 10)
