# STAT 405 — Final Project
## Bayesian Inference for Diabetes Risk Prediction

---
Project repository:
https://github.com/thecryptic2107/stat-405-final-project.git
## 👥 Team Members

| Name | Student ID | Contributions |
|---|---|---|
| Maahi Gumber | 76901438 | TO BE DONE!!!!!! |
| Rishabh Mathur | 71497275 | TO BE DONE!!!!!! |
| Ayushman | — |  |


## 🗄️ Datasets

### Primary Dataset — Diabetes 130-US Hospitals (1999–2008)

> **101,766 patient encounters × 47 features**

**Source:** UCI Machine Learning Repository — CC BY 4.0 License  
**URL:** https://archive.ics.uci.edu/dataset/296/diabetes+130-us+hospitals+for+years+1999-2008  
**Access:** Direct CSV download, no registration required.

### Backup Dataset — CDC Diabetes Health Indicators (BRFSS 2015)

> **253,680 survey responses × 22 features**

**Source:** CDC Behavioral Risk Factor Surveillance System (BRFSS) 2015, cleaned by Alex Teboul  
**URL:** https://www.kaggle.com/datasets/alexteboul/diabetes-health-indicators-dataset 

**Description:**  
A cleaned subset of the CDC's 2015 BRFSS — the largest continuously conducted health survey
in the world. Contains 253,680 US adult respondents with 21 lifestyle, demographic, and
clinical indicator features. Outcome: **binary diabetes diagnosis**. Class distribution:
86.1% non-diabetic, 13.9% diabetic/prediabetic.

## 🔬 Scientific Questions & Methodology

---

### Scientific Question 1 — Primary Dataset (130-US Hospitals)

**Question:** What is the posterior probability of 30-day hospital readmission for a
diabetic patient given their clinical profile — and does this probability differ credibly
across racial and age subgroups?

This is a strictly Bayesian question. We want the full posterior distribution over
readmission risk, not a point prediction. Specifically: does the effect of HbA1c
measurement on readmission risk vary across subgroups in a way that is credibly
non-zero — i.e., is P(β_HbA1c < 0 | data) large enough to support clinical action?

---
## Methodology – Primary Dataset (130 US Hospitals)

We model the probability that a diabetic patient is readmitted within 30 days using Bayesian logistic regression.

Model 1 is a baseline pooled logistic regression where all patients are treated as exchangeable. The probability of readmission for patient i is

p_i = sigmoid(β0 + Σ βj x_ij)

and the observed outcome follows

readmitted_i ~ Bernoulli(p_i)

We place weakly informative priors on the coefficients

β0 ~ Normal(0, 2)  
βj ~ Normal(0, 1)

Posterior inference for this model is performed using Metropolis Hastings MCMC. We run multiple chains and evaluate convergence using trace plots, R-hat, and effective sample size.

Model 2 extends this to a hierarchical logistic regression to capture subgroup variation across race and age groups. Each subgroup has its own intercept

β_g ~ Normal(μβ, τ²)

with hyperpriors

μβ ~ Normal(0,1)  
τ ~ HalfNormal(0.5)

This allows partial pooling so that small subgroups borrow strength from the overall population.

Posterior inference for the hierarchical model is performed using variational inference with stochastic gradient optimization. The variational distribution approximates the true posterior by minimizing the KL divergence. Gradients are estimated using the reparameterization trick and mini batch sampling.

We compare the two methods by examining posterior distributions, predictive performance, and calibration of credible intervals.

### Scientific Question 2 — Backup Dataset (CDC BRFSS 2015)

**Question:** What is the posterior distribution of diabetes risk as a function of BMI
and age — and does the strength of the BMI effect differ credibly across income and
education subgroups?

## Methodology – Backup Dataset (CDC BRFSS 2015)

For the BRFSS dataset we model the probability of diabetes diagnosis as a function of BMI, age, and other health indicators using Bayesian logistic regression.

The baseline model is a pooled logistic regression

p_i = sigmoid(β0 + Σ βj x_ij)

Diabetes_i ~ Bernoulli(p_i)

where x_ij includes predictors such as BMI, age, income, education, and lifestyle indicators.

Priors are chosen to be weakly informative

β0 ~ Normal(0,2)  
βj ~ Normal(0,1)

Posterior sampling for the baseline model is performed using Metropolis Hastings MCMC to obtain samples from the posterior distribution of the regression coefficients.

To study subgroup heterogeneity we build a hierarchical model where the BMI effect can vary across demographic groups defined by age and income.

β_BMI_g ~ Normal(μβ, τ²)

with hyperpriors

μβ ~ Normal(0,1)  
τ ~ HalfNormal(0.5)

This hierarchical structure allows different demographic groups to have different BMI effects while still sharing information across groups.

Posterior inference for this model is performed using variational inference with stochastic gradient descent. The variational approximation provides a scalable alternative to MCMC for large datasets.

Results from both approaches are compared using posterior summaries, predictive accuracy, and posterior predictive checks.
