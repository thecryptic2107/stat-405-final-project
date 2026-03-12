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

### Scientific Question 2 — Backup Dataset (CDC BRFSS 2015)

**Question:** What is the posterior distribution of diabetes risk as a function of BMI
and age — and does the strength of the BMI effect differ credibly across income and
education subgroups?
