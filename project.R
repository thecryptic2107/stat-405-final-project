library(dplyr)
library(ggplot2)
suppressPackageStartupMessages(library(cmdstanr))

set.seed(405)

raw_data <- read.csv("~/Downloads/diabetes+130-us+hospitals+for+years+1999-2008/diabetic_data.csv", header=FALSE)

colnames(raw_data) <- raw_data[1, ]
# remove the patients that could not be included because of death as shown in th epaper

hospice_death_ids <- c(11, 3, 14, 19, 20, 21)

diabetic_data <- raw_data %>%
                  filter(!discharge_disposition_id %in% hospice_death_ids)

# to ensure that observation are independent, we only use ONe encounter per patient

diabetic_data <- diabetic_data %>%
  arrange(patient_nbr) %>%
  group_by(patient_nbr) %>%
  slice(1) %>%
  ungroup


# for our simple model we create a binary repsonse variable
# 1: if the readmission withing 30 days
# 0: if there was no readmission or >30 days

diabetic_data <-  diabetic_data %>%
  mutate(readmitted_30 = ifelse(readmitted == "<30", 1, 0))


# for Model 1 we are creating a binary encoding for HbA1c testing for simplicity
# None = not tested at all (0)
# ">7" ">8" = tested (1)

diabetic_data <- diabetic_data %>%
  mutate(hba1c_tested = ifelse(is.na(A1Cresult) | A1Cresult == "None", 0 , 1))


# for the age group we make three distinct age intervals
# "young" = < 30
# "middle" = 30-60
# old" = 60+

diabetic_data <- diabetic_data %>%
  mutate(age_group = case_when(
    age %in% c("[0-10)" , "[10-20)", "[20-30)") ~ "young",     # < 30
      age %in% c("[30-40)", "[40-50)", "[50-60)") ~ "middle",    # 30-60
      TRUE                                        ~ "old"      # >60
  ),
  # Reference level: "young" - encoded as two dummies
  age_middle = ifelse(age_group == "middle", 1, 0),
  age_old    = ifelse(age_group == "old",    1, 0)
  )

# standardize time in hospital
diabetic_data$time_in_hospital <- as.numeric(diabetic_data$time_in_hospital)

diabetic_data <- diabetic_data %>%
  mutate(
    time_in_hospital = scale(time_in_hospital)[, 1]
  )


# finally we only select the variables we plan to use for the model

model1_vars <- c("readmitted_30", "hba1c_tested", "age_middle",
                 "age_old", "time_in_hospital")

diabetic_model <- diabetic_data %>%
  select(all_of(model1_vars)) %>%
  na.omit()

# --------DATA CLEANING DONE FOR MODEL 1 -------

# -------MODEL 1: MCMC ------------

y <- diabetic_model$readmitted_30
x <- model.matrix(~ hba1c_tested, data = diabetic_model)
n <- nrow(x)
p <- ncol(x)


mod <- cmdstan_model("project.stan")

stan_data <- list(
  N             = n,
  hba1c_tested  = diabetic_model$hba1c_tested,
  readmitted_30 = y
)

fit <- mod$sample(
  data         = stan_data,
  seed         = 405,
  chains       = 1,
  iter_sampling= 2500,
  refresh      = 500
)


cat("=== Posterior Summary ===\n")
print(fit$summary(c("beta0", "beta1")))

draws_df   <- fit$draws(format = "df")
beta1_draws <- draws_df$beta1

cat(sprintf("\nP(β₁ < 0 | y) = %.3f  [probability testing is protective]\n",
            mean(beta1_draws < 0)))
cat(sprintf("P(β₁ > 0 | y) = %.3f  [probability testing increases risk]\n",
            mean(beta1_draws > 0)))
cat(sprintf("Posterior mean OR = %.3f  (< 1 = protective)\n",
            exp(mean(beta1_draws))))

# CONCLUSION: "Bayesian logistic regression via MCMC showed strong evidence 
# that HbA1c testing during hospitalization reduces 30-day readmission risk 
# (β₁ = -0.103, OR = 0.903, 95% CI: [-0.161, -0.041], P(protective) = 99.7%)."


