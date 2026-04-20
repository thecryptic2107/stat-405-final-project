library(dplyr)
library(ggplot2)
library(tidyr)
suppressPackageStartupMessages(library(cmdstanr))

set.seed(405)

# Pre-processing the data

raw_data <- read.csv("diabetic_data.csv", header=FALSE)

colnames(raw_data) <- raw_data[1, ]

hospice_death_ids <- c(11, 3, 14, 19, 20, 21)

diabetic_data <- raw_data %>%
  filter(!discharge_disposition_id %in% hospice_death_ids)

diabetic_data <- diabetic_data %>%
  arrange(patient_nbr) %>%
  group_by(patient_nbr) %>%
  slice(1) %>%
  ungroup

diabetic_data <- diabetic_data %>%
  mutate(readmitted_30 = ifelse(readmitted == "<30", 1, 0))

diabetic_data <- diabetic_data %>%
  mutate(hba1c_tested = ifelse(is.na(A1Cresult) | A1Cresult == "None", 0, 1))

diabetic_data <- diabetic_data %>%
  mutate(
    age_group  = case_when(
      age %in% c("[0-10)", "[10-20)", "[20-30)") ~ "young",
      age %in% c("[30-40)", "[40-50)", "[50-60)") ~ "middle",
      TRUE                                        ~ "old"
    ),
    age_middle = ifelse(age_group == "middle", 1, 0),
    age_old    = ifelse(age_group == "old",    1, 0)
  )

diabetic_data$time_in_hospital <- as.numeric(diabetic_data$time_in_hospital)

diabetic_data <- diabetic_data %>%
  mutate(time_in_hospital = scale(time_in_hospital)[, 1])

model1_vars <- c("readmitted_30", "hba1c_tested", "age_middle", "age_old", "time_in_hospital")

diabetic_model <- diabetic_data %>%
  select(all_of(model1_vars)) %>%
  na.omit()

y <- diabetic_model$readmitted_30

# Defining all shared functions

log_posterior <- function(beta) {
  eta <- as.numeric(X %*% beta)
  ll  <- sum(y * eta - log1p(exp(eta)))
  lp  <- sum(dnorm(beta, 0, prior_sd, log = TRUE))
  ll + lp
}

component_mh <- function(n_iter,
                          prop_sd     = c(0.30, 0.15, 0.15, 0.15, 0.10),
                          burn_in     = 3000,
                          param_names = c("beta0", "beta1", "beta2", "beta3", "beta4")) {
  chain    <- matrix(NA_real_, nrow = n_iter, ncol = p)
  colnames(chain) <- param_names[seq_len(p)]
  beta_cur <- rep(0, p)
  lp_cur   <- log_posterior(beta_cur)
  accepted <- integer(p)

  for (t in seq_len(n_iter)) {
    for (j in seq_len(p)) {
      beta_prop    <- beta_cur
      beta_prop[j] <- beta_cur[j] + rnorm(1, 0, prop_sd[j])
      lp_prop <- log_posterior(beta_prop)
      if (log(runif(1)) < lp_prop - lp_cur) {
        beta_cur    <- beta_prop
        lp_cur      <- lp_prop
        accepted[j] <- accepted[j] + 1L
      }
    }
    chain[t, ] <- beta_cur
  }

  cat("Per-component acceptance rates:\n")
  rates <- accepted / n_iter
  names(rates) <- param_names[seq_len(p)]
  print(round(rates, 3))

  posterior <- chain[(burn_in + 1):n_iter, , drop = FALSE]
  list(chain = chain, posterior = posterior, accept_rates = rates)
}

ess_geyer <- function(x) {
  ac    <- acf(x, plot = FALSE, lag.max = 500)$acf[-1]
  pairs <- ac[seq(1, length(ac)-1, 2)] + ac[seq(2, length(ac), 2)]
  k     <- which(pairs <= 0)[1]
  if (is.na(k)) k <- length(pairs)
  length(x) / (1 + 2 * sum(ac[seq_len(2*k - 1)]))
}

# 
# MODEL 1
# 

stan_data <- list(
  N             = nrow(diabetic_model),
  hba1c_tested  = diabetic_model$hba1c_tested,
  readmitted_30 = y
)

# Method 1a: Component-wise MH

X        <- cbind(intercept = 1, hba1c_tested = diabetic_model$hba1c_tested)
p        <- ncol(X)
prior_sd <- c(2, 1)

set.seed(405)
result_m1 <- component_mh(
  n_iter      = 20000,
  prop_sd     = c(0.04, 0.09),
  burn_in     = 3000,
  param_names = c("beta0", "beta1")
)
post_m1       <- result_m1$posterior
post_means_m1 <- colMeans(post_m1)
post_ci_m1    <- apply(post_m1, 2, quantile, probs = c(0.025, 0.975))
ess_m1        <- apply(post_m1, 2, ess_geyer)

print(round(post_means_m1, 4))
print(round(post_ci_m1, 4))
print(round(ess_m1, 1))
cat(sprintf("\nP(β₁ < 0 | y) = %.3f\n", mean(post_m1[, "beta1"] < 0)))
cat(sprintf("Posterior mean OR = %.3f\n", exp(post_means_m1["beta1"])))

# Trace plot
chain_df_m1 <- as.data.frame(result_m1$chain)
chain_df_m1$iteration <- seq_len(20000)
chain_df_m1$phase     <- ifelse(chain_df_m1$iteration <= 3000, "Burn-in", "Sampling")

param_labels_m1 <- c(beta0 = "β₀ (intercept)", beta1 = "β₁ (hba1c_tested)")

chain_long_m1 <- chain_df_m1 %>%
  pivot_longer(cols = c(beta0, beta1), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = names(param_labels_m1)))

p_trace_m1 <- ggplot(chain_long_m1, aes(x = iteration, y = value, colour = phase)) +
  geom_line(linewidth = 0.25, alpha = 0.8) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 1,
             labeller = as_labeller(param_labels_m1)) +
  scale_colour_manual(values = c("Burn-in" = "grey60", "Sampling" = "steelblue")) +
  labs(title    = "Model 1 – Component-wise MH: Trace Plots",
       subtitle = "Grey = discarded burn-in  |  Blue = post-burn-in samples used for inference",
       x = "Iteration", y = "Value", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave("model1_mh_trace.png", p_trace_m1, width = 9, height = 6, dpi = 150)

# Method 1b: ADVI (using logistic_m1.stan)

mod_m1_vi <- cmdstan_model("logistic_m1.stan")

fit_vi_m1 <- mod_m1_vi$variational(
  data           = stan_data,
  seed           = 405,
  algorithm      = "meanfield",
  iter           = 10000,
  output_samples = 2000,
  adapt_engaged  = FALSE,
  eta            = 1
)
print(fit_vi_m1$summary(c("beta0", "beta1")))

vi_m1_draws <- fit_vi_m1$draws(variables = c("beta0", "beta1"), format = "df")
cat(sprintf("\nVI P(β₁ < 0 | y) = %.3f\n",    mean(vi_m1_draws$beta1 < 0)))
cat(sprintf("VI Posterior mean OR = %.3f\n", exp(mean(vi_m1_draws$beta1))))

# Model 1: MH vs ADVI comparison

mh_b1_q <- quantile(post_m1[, "beta1"], probs = c(0.025, 0.975))
vi_b1_q <- quantile(vi_m1_draws$beta1,  probs = c(0.025, 0.975))

m1_comparison <- data.frame(
  Method     = c("Component-wise MH", "ADVI (mean-field)"),
  beta1_mean = round(c(post_means_m1["beta1"],   mean(vi_m1_draws$beta1)), 4),
  beta1_sd   = round(c(sd(post_m1[, "beta1"]),   sd(vi_m1_draws$beta1)),   4),
  CI_lower   = round(c(mh_b1_q[1], vi_b1_q[1]), 4),
  CI_upper   = round(c(mh_b1_q[2], vi_b1_q[2]), 4),
  OR         = round(exp(c(post_means_m1["beta1"], mean(vi_m1_draws$beta1))), 4),
  stringsAsFactors = FALSE
)

print(m1_comparison)

# Density comparison — MH vs ADVI

mh_df_m1 <- data.frame(beta0 = post_m1[, "beta0"],  beta1 = post_m1[, "beta1"], method = "Component-wise MH")
vi_df_m1 <- data.frame(beta0 = vi_m1_draws$beta0,   beta1 = vi_m1_draws$beta1,  method = "ADVI (mean-field)")

m1_all_long <- bind_rows(mh_df_m1, vi_df_m1) %>%
  pivot_longer(cols = c(beta0, beta1), names_to = "parameter", values_to = "value") %>%
  mutate(
    method    = factor(method, levels = c("Component-wise MH", "ADVI (mean-field)")),
    param_lbl = param_labels_m1[parameter]
  )

p_m1_compare <- ggplot(m1_all_long, aes(x = value, fill = method, colour = method)) +
  geom_density(alpha = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  facet_wrap(~ param_lbl, scales = "free", ncol = 2) +
  scale_fill_manual(values   = c("Component-wise MH" = "forestgreen", "ADVI (mean-field)" = "darkorange")) +
  scale_colour_manual(values = c("Component-wise MH" = "darkgreen",   "ADVI (mean-field)" = "darkorange3")) +
  labs(title    = "Model 1 — Posterior Comparison: Component-wise MH vs ADVI",
       subtitle = "Green = Component-wise MH  |  Orange = ADVI (mean-field)",
       x = "Parameter value", y = "Density", fill = "Method", colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("model1_method_comparison.png", p_m1_compare, width = 10, height = 5, dpi = 150)

# 
# MODEL 2
# 

X <- cbind(
  intercept        = 1,
  hba1c_tested     = diabetic_model$hba1c_tested,
  age_middle       = diabetic_model$age_middle,
  age_old          = diabetic_model$age_old,
  time_in_hospital = diabetic_model$time_in_hospital
)
N        <- nrow(X)
p        <- ncol(X)
prior_sd <- c(2, 1, 1, 1, 1)

# Method 2a: Component-wise MH 

n_iter  <- 20000
burn_in <- 3000

set.seed(405)
result <- component_mh(n_iter = n_iter, burn_in = burn_in)
post   <- result$posterior

post_means <- colMeans(post)
print(round(post_means, 4))

post_ci <- apply(post, 2, quantile, probs = c(0.025, 0.975))
print(round(post_ci, 4))

cat(sprintf("P(β₁ < 0 | y) = %.3f  [HbA1c testing is protective]\n",    mean(post[, "beta1"] < 0)))
cat(sprintf("P(β₄ > 0 | y) = %.3f  [longer stay → higher readmission]\n", mean(post[, "beta4"] > 0)))
cat(sprintf("Posterior mean OR for HbA1c = %.3f\n", exp(post_means["beta1"])))

ess_vals <- apply(post, 2, ess_geyer)
print(round(ess_vals, 1))

param_labels <- c(
  beta0 = "β₀ (intercept)",
  beta1 = "β₁ (hba1c_tested)",
  beta2 = "β₂ (age_middle)",
  beta3 = "β₃ (age_old)",
  beta4 = "β₄ (time_in_hospital)"
)

chain_df <- as.data.frame(result$chain)
colnames(chain_df) <- names(param_labels)
chain_df$iteration <- seq_len(n_iter)
chain_df$phase     <- ifelse(chain_df$iteration <= burn_in, "Burn-in", "Sampling")

chain_long <- pivot_longer(chain_df,
                           cols = starts_with("beta"), names_to = "parameter", values_to = "value")
chain_long$parameter <- factor(chain_long$parameter, levels = names(param_labels))

p_trace <- ggplot(chain_long, aes(x = iteration, y = value, colour = phase)) +
  geom_line(linewidth = 0.25, alpha = 0.8) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 1,
             labeller = as_labeller(param_labels)) +
  scale_colour_manual(values = c("Burn-in" = "grey60", "Sampling" = "steelblue")) +
  labs(title    = "Model 2 – Component-wise MH: Trace Plots",
       subtitle = "Grey = discarded burn-in  |  Blue = post-burn-in samples used for inference",
       x = "Iteration", y = "Value", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave("model2_trace_plots.png", p_trace, width = 9, height = 11, dpi = 150)

post_df   <- as.data.frame(post)
post_long <- pivot_longer(post_df, cols = everything(),
                          names_to = "parameter", values_to = "value")
post_long$parameter <- factor(post_long$parameter, levels = names(param_labels))

p_dens <- ggplot(post_long, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.45, colour = "steelblue4") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.7) +
  facet_wrap(~ parameter, scales = "free", ncol = 3,
             labeller = as_labeller(param_labels)) +
  labs(title    = "Model 2 – Posterior Densities (Component-wise MH)",
       subtitle = "Dashed red = 0",
       x = "Parameter value", y = "Density") +
  theme_bw(base_size = 11)

ggsave("model2_posterior_densities.png", p_dens, width = 10, height = 6, dpi = 150)

running_mean_b1 <- cumsum(post[, "beta1"]) / seq_along(post[, "beta1"])

p_lln <- ggplot(data.frame(t = seq_along(running_mean_b1), rm = running_mean_b1),
                aes(x = t, y = rm)) +
  geom_line(colour = "steelblue", linewidth = 0.6) +
  geom_hline(yintercept = post_means["beta1"],
             linetype = "dashed", colour = "darkred", linewidth = 0.8) +
  annotate("text",
           x     = length(running_mean_b1) * 0.55,
           y     = post_means["beta1"] + 0.003,
           label = sprintf("Posterior mean = %.4f", post_means["beta1"]),
           colour = "darkred", size = 3.5) +
  labs(title    = "Markov LLN: Running Mean of β₁ (hba1c_tested)",
       subtitle = "Time-average converges to posterior expectation — topic03",
       x = "Post-burn-in iteration t",
       y = "(1/t) Σ β₁") +
  theme_bw(base_size = 12)

ggsave("model2_markov_lln.png", p_lln, width = 8, height = 4, dpi = 150)

summary_df <- data.frame(
  Parameter = unname(param_labels),
  Mean      = round(post_means, 4),
  `2.5%`    = round(post_ci[1, ], 4),
  `97.5%`   = round(post_ci[2, ], 4),
  ESS       = round(ess_vals, 0),
  check.names = FALSE
)
print(summary_df, row.names = FALSE)

cat(sprintf("  β₁ (HbA1c):        mean = %.3f, OR = %.3f, P(protective) = %.1f%%\n",
            post_means["beta1"], exp(post_means["beta1"]), 100 * mean(post[, "beta1"] < 0)))
cat(sprintf("  β₃ (age_old):       mean = %.3f\n", post_means["beta3"]))
cat(sprintf("  β₄ (time in hosp):  mean = %.3f\n", post_means["beta4"]))

# Method 2b: ADVI (using logistic_full.stan)

stan_data_full <- list(
  N                = nrow(diabetic_model),
  hba1c_tested     = diabetic_model$hba1c_tested,
  age_middle       = diabetic_model$age_middle,
  age_old          = diabetic_model$age_old,
  time_in_hospital = diabetic_model$time_in_hospital,
  readmitted_30    = diabetic_model$readmitted_30
)

mod_full <- cmdstan_model("logistic_full.stan")

fit_vi <- mod_full$variational(
  data           = stan_data_full,
  seed           = 405,
  algorithm      = "meanfield",
  iter           = 10000,
  output_samples = 2000
)

vi_params <- c("beta0", "beta1", "beta2", "beta3", "beta4")

print(fit_vi$summary(vi_params))

vi_draws_raw <- fit_vi$draws(variables = vi_params, format = "df")

cat(sprintf("\nVI  P(β₁ < 0 | y) = %.3f\n",       mean(vi_draws_raw$beta1 < 0)))
cat(sprintf("VI  P(β₄ > 0 | y) = %.3f  [longer stay increases readmission risk]\n", mean(vi_draws_raw$beta4 > 0)))
cat(sprintf("VI  Posterior mean OR (HbA1c) = %.3f\n", exp(mean(vi_draws_raw$beta1))))

# Model 2: MH vs VI comparison

vi_mat <- as.matrix(vi_draws_raw[, vi_params])
mh_mat <- post[, vi_params]

vi_q <- apply(vi_mat, 2, quantile, probs = c(0.025, 0.975))
mh_q <- apply(mh_mat, 2, quantile, probs = c(0.025, 0.975))

comparison_df <- data.frame(
  Parameter = vi_params,
  MH_mean   = round(colMeans(mh_mat),  4),
  VI_mean   = round(colMeans(vi_mat),  4),
  MH_sd     = round(apply(mh_mat, 2, sd), 4),
  VI_sd     = round(apply(vi_mat, 2, sd), 4),
  MH_95lo   = round(mh_q[1, ], 4),
  MH_95hi   = round(mh_q[2, ], 4),
  VI_95lo   = round(vi_q[1, ], 4),
  VI_95hi   = round(vi_q[2, ], 4),
  stringsAsFactors = FALSE
)

print(comparison_df)

param_labels_full <- c(
  beta0 = "β₀ (Intercept)",
  beta1 = "β₁ (HbA1c tested)",
  beta2 = "β₂ (Age: middle)",
  beta3 = "β₃ (Age: old)",
  beta4 = "β₄ (Time in hospital)"
)

mh_long <- as.data.frame(mh_mat) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(method = "Component-wise MH")

vi_long <- as.data.frame(vi_mat) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(method = "VI (mean-field)")

all_draws_long <- bind_rows(mh_long, vi_long) %>%
  mutate(
    method    = factor(method, levels = c("Component-wise MH", "VI (mean-field)")),
    param_lbl = param_labels_full[parameter]
  )

density_plot <- ggplot(all_draws_long, aes(x = value, fill = method, colour = method)) +
  geom_density(alpha = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  facet_wrap(~ param_lbl, scales = "free", ncol = 3) +
  scale_fill_manual(values   = c("Component-wise MH" = "steelblue",  "VI (mean-field)" = "darkorange")) +
  scale_colour_manual(values = c("Component-wise MH" = "steelblue4", "VI (mean-field)" = "darkorange3")) +
  labs(
    title    = "Model 2 — Posterior Comparison: Component-wise MH vs ADVI",
    subtitle = "Blue = Component-wise MH  |  Orange = ADVI (mean-field)",
    x = "Parameter value", y = "Density", fill = "Method", colour = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "top")

ggsave("model2_mh_vs_vi_densities.png", density_plot, width = 10, height = 7, dpi = 150)

forest_df <- bind_rows(
  data.frame(parameter = vi_params, method = "Component-wise MH",
             mean = comparison_df$MH_mean, lo = comparison_df$MH_95lo, hi = comparison_df$MH_95hi,
             stringsAsFactors = FALSE),
  data.frame(parameter = vi_params, method = "VI (mean-field)",
             mean = comparison_df$VI_mean, lo = comparison_df$VI_95lo, hi = comparison_df$VI_95hi,
             stringsAsFactors = FALSE)
) %>%
  mutate(
    method    = factor(method, levels = c("Component-wise MH", "VI (mean-field)")),
    param_lbl = factor(param_labels_full[parameter], levels = rev(unname(param_labels_full)))
  )

forest_plot <- ggplot(forest_df, aes(x = mean, y = param_lbl, colour = method, shape = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.5), size = 0.6) +
  scale_colour_manual(values = c("Component-wise MH" = "steelblue", "VI (mean-field)" = "darkorange")) +
  labs(title  = "Forest Plot: Posterior Means + 95% CIs — MH vs VI",
       x = "Parameter value", y = NULL, colour = "Method", shape = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("model2_mh_vs_vi_forest.png", forest_plot, width = 8, height = 5, dpi = 150)

beta1_mh <- mh_mat[, "beta1"]
beta1_vi <- vi_draws_raw$beta1

beta1_summary <- data.frame(
  Method   = c("Component-wise MH", "VI (mean-field)"),
  Mean     = round(c(mean(beta1_mh),  mean(beta1_vi)),  4),
  SD       = round(c(sd(beta1_mh),    sd(beta1_vi)),    4),
  CI_lower = round(c(quantile(beta1_mh, 0.025, names = FALSE),
                     quantile(beta1_vi, 0.025, names = FALSE)), 4),
  CI_upper = round(c(quantile(beta1_mh, 0.975, names = FALSE),
                     quantile(beta1_vi, 0.975, names = FALSE)), 4)
)

print(beta1_summary)
