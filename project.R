library(dplyr)
library(ggplot2)

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
# ">7" ">8" = teste (1)

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
  # Reference level: "young" — encoded as two dummies
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