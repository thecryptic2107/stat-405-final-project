# STAT 405 — Final Project
## Bayesian Inference for Diabetes Risk Prediction
### Hierarchical Modeling | SGD-VI vs. MCMC | Hospital & Population-Level Datasets

---

## 👥 Team Members

| Name | Student ID | Contributions |
|---|---|---|
| Maahi Gumber | 76901438 | TO BE DONE!!!!!! |
| Rishabh | — |  |
| Ayushman | — |  |


## 🗄️ Datasets

### Primary Dataset — Diabetes 130-US Hospitals (1999–2008)

> **101,766 patient encounters × 47 features**

**Source:** UCI Machine Learning Repository — CC BY 4.0 License  
**URL:** https://archive.ics.uci.edu/dataset/296/diabetes+130-us+hospitals+for+years+1999-2008  
**Access:** Direct CSV download, no registration required.

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

### Scientific Question 1 — Backup Dataset (CDC BRFSS 2015)

**Question:** What is the posterior distribution of diabetes risk as a function of BMI
and age — and does the strength of the BMI effect differ credibly across income and
education subgroups?
