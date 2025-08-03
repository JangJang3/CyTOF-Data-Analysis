################################################################
# Set working environment 
################################################################

current_dir <- "/mnt/volumeb/jin/survival_analysis"
getwd()
setwd(dir= current_dir)


################################################################
# Loading Packages 
################################################################

library(DT)
library(dplyr)
library(tidyr)
library(openxlsx) 
library(ggplot2)


################################################################
# Loading Clinical File
################################################################

# Cohort with all HCC patients collected from Tübingen 

clinical_data <- openxlsx::read.xlsx("/mnt/volumeb/jin/HCC patient clinical data July 2025 updated.xlsx", sheet = 2, rows = 1:41, colNames = TRUE, skipEmptyRows = TRUE, detectDates = TRUE)

head(clinical_data)


################################################################
# Check File for type date values 
################################################################

# Because some dates are converted into an integer, we need to check 

View(clinical_data)

# Date_first_application, death_date and progression_date include date values 

is_date_format <- function(col){
  all(grepl("^\\d{4}-\\d{2}-\\d{2}$", col)) # check if the values match YYYY-MM-DD
}

# Check within the columns 
is_date_format(clinical_data$Date_first_application) # TRUE
is_date_format(clinical_data$death_date)  # FALSE
is_date_format(clinical_data$progression_date)  # FALSE
is_date_format(clinical_data$last_contact) # TRUE

clinical_data$death_date

################################################################
# Convert the Values to Date values
################################################################

# Change ? to NA, because there are no further information about the patients death and progression 

clinical_data <- clinical_data %>% 
  mutate(death_date = na_if(death_date, "?")) # Replace ? with NA


clinical_data$death_date # check results 


#' Convert the date values e.g. "44228" to actual dates
#' For that we have to return every value of class "Date"
#' Calculates the square of the input
#'
#' @param x a value of type integer or Date
#' @return a Date object 

  
convert_values_date <- function(x) {
  if (is.na(x)) {
    return(as.Date(NA))
  } else if (grepl("^\\d+$", x)) {
    print(x)
    print(convertToDate(x))
    date_value <- convertToDate(as.numeric(x))
    return (date_value)
  } else {
    return (as.Date(x))
  }
}

clinical_data$death_date <- lapply(clinical_data$death_date, convert_values_date) 

# To ensure that the "date" types are kept, I will use do.call to 
# Because with unlist(...) we get an vector
# And a vector only differentiates between character, integers 

clinical_data$death_date <- do.call(c, clinical_data$death_date)


clinical_data$death_date # check results 


# now for progression_date 
clinical_data$progression_date

# clinical_data <- clinical_data %>% 
#   mutate(progression_date = na_if(progression_date, "?")) # Replace ? with NA

clinical_data$progression_date <- lapply(clinical_data$progression_date, convert_values_date)

clinical_data$progression_date <- do.call(c, clinical_data$progression_date)

clinical_data$progression_date


################################################################
# Extract important information 
################################################################

# Columns to keep:
# PatID, Age, Sex, BMI, viral, alcohol, diabetes, Date_first_application
# Staging_3months, alt_staging_3months, Staging_6months
# death_date, progression_date, OS_11_24, PFS_11_24

keep_features <- c("PatID", "Age", "Sex", "BMI_kg_m2", "viral", "Alcohol", "Diabetes_mellitus_I", "Diabetes_mellitus_II", "Date_first_application",
                   "Staging_3months", "altStaging_3months", "Staging_6months", "OS_11_24", "PFS_11_24", "death_date", "progression_date", "last_contact")

clinical_data <- clinical_data[keep_features]


################################################################
# Change Patient ID
################################################################

# We need to change it to be compatible to the metadata of the samples

clinical_data$PatID <- paste("Pat_ID_", clinical_data$PatID, sep = "")

clinical_data$PatID


################################################################
# Change column PatID to patient_id
################################################################

# To keep the same key for meta data, change column PatID 
clinical_data <- clinical_data %>% 
  rename(patient_id = PatID)


################################################################
# Extra specific PatIDs 
################################################################

# used_IDs <- c(1, 5, 7, 8, 11, 12, 13, 9, 14, 15, 18, 19, 21, 22, 23, 24, 25, 26, 29, 30, 33, 36, 40, 43, 16, 42, 44)

# test <- c()
# for(i in used_IDs){
#  test <- c(test, (paste("PatID_", i, sep = "")))
# }

used_IDs <- c("Pat_ID_1", "Pat_ID_5", "Pat_ID_7", "Pat_ID_8", "Pat_ID_11", "Pat_ID_12", "Pat_ID_13", "Pat_ID_9", "Pat_ID_14", "Pat_ID_15", "Pat_ID_18", "Pat_ID_19", "Pat_ID_21", "Pat_ID_22", 
              "Pat_ID_23", "Pat_ID_24", "Pat_ID_25", "Pat_ID_26", "Pat_ID_29", "Pat_ID_30", "Pat_ID_33", "Pat_ID_36", "Pat_ID_40", "Pat_ID_43", "Pat_ID_16", "Pat_ID_42", "Pat_ID_44")

clinical_data <- clinical_data[clinical_data$patient_id %in% used_IDs, ]


################################################################
# Add survival information 
################################################################

# Each HCC patient is grouped into survival short term, mid term or long term. 

patient_survival <- c(
  "Pat_ID_1"  = "Long term",
  "Pat_ID_5"  = "Long term",
  "Pat_ID_12" = "Long term",
  "Pat_ID_13" = "Long term",
  "Pat_ID_14" = "Long term",
  "Pat_ID_29" = "Long term",
  "Pat_ID_33" = "Long term",
  "Pat_ID_36" = "Long term",
  "Pat_ID_43" = "Long term",
  
  "Pat_ID_8"  = "Mid term",
  "Pat_ID_15" = "Mid term",
  "Pat_ID_16" = "Mid term",
  "Pat_ID_23" = "Mid term",
  "Pat_ID_24" = "Mid term",
  "Pat_ID_25" = "Mid term",
  "Pat_ID_30" = "Mid term",
  "Pat_ID_42" = "Mid term",
  "Pat_ID_44" = "Mid term",
  
  "Pat_ID_7"  = "Short term",
  "Pat_ID_9"  = "Short term",
  "Pat_ID_11" = "Short term",
  "Pat_ID_18" = "Short term",
  "Pat_ID_19" = "Short term",
  "Pat_ID_21" = "Short term",
  "Pat_ID_22" = "Short term",
  "Pat_ID_26" = "Short term",
  "Pat_ID_40" = "Short term"
)

#setdiff(names(patient_survival), clinical_data$patient_id)
clinical_data$survival_status <- patient_survival[clinical_data$patient_id]


################################################################
# Save Clinical data 
################################################################
saveRDS(clinical_data, "new_clinical_data_filtered_July.rds")
#clinical_data <- readRDS("new_clinical_data_filtered_July.rds")


################################################################
# Create Plots
################################################################

# Histogram of Age Distribution
hist(clinical_data$Age, main ="Age Distribution of Patients", xlab = "Age", col = "skyblue", border = "black")

# There is a "m " with space in the data 
clinical_data$Sex <- trimws(clinical_data$Sex)

# Boxplot of Age by Gender 
boxplot(Age ~ Sex, data=clinical_data, main = "Age by Sex", xlab = "Sex", ylab = "Age", col=c("pink", "skyblue"))


# Boxplot of BMI by Gender
boxplot(BMI_kg_m2 ~ Sex, data=clinical_data, main = "BMI by Sex", xlab = "Sex", ylab = "BMI", col=c("pink", "skyblue"))



# Grouped barchart of gender, conditions at Staging 6 months 
conditions <- clinical_data %>% 
  pivot_longer(cols = c(viral, Alcohol, Diabetes_mellitus_I, Diabetes_mellitus_II), names_to = "Preconditions", values_to = "Status") 

ggplot(conditions, aes(fill = Preconditions, x = Staging_6months)) +
  geom_bar(position = "dodge", stat = "count")


# Responder vs Non-Responders by Preconditions (only HCC patients with preconditions)
conditions_yes <- clinical_data %>% 
  pivot_longer(cols = c(viral, Alcohol, Diabetes_mellitus_I, Diabetes_mellitus_II), names_to = "Preconditions", values_to = "Status") %>%
  filter(Status == "y")

ggplot(conditions_yes, aes(fill = Preconditions, x = Staging_6months)) +
  geom_bar(position = "dodge", stat = "count") +
  geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 1), vjust = -0.5, size = 2) + 
  scale_fill_discrete(name = "Preconditions") +
  labs(title = "Responder vs Non-Responder by preconditions", x = "Status", y = "Count")


# Responder vs Non-Responders by Preconditions (only HCC patients with preconditions) by Gender
ggplot(conditions_yes, aes(fill = Preconditions, x = Sex)) +
  geom_bar(position = "dodge", stat = "count") +
  geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 1), vjust = -0.5, size = 2) + 
  scale_fill_discrete(name = "Preconditions") +
  labs(title = "Responder vs Non-Responder by preconditions", x = "Gender", y = "Count")



# Histogram of Responder vs Non-Responder 
ggplot(clinical_data, aes(x = Staging_6months, fill = Staging_6months)) + 
  geom_bar() + 
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 2) +
  labs(title = "Responder Vs Non-Responder Distribution of Patients", x = "Status", y = "Count")
  


################################################################
# Survival Analysis
################################################################
#install.packages(c("survival", "survminer"))

library(survival)
library(survminer)
library(lubridate) 

# Create Survival Status at 6 months 
ggplot(clinical_data, aes(x = Staging_6months, fill = survival_status)) +
  geom_bar(stat = "count", position = position_dodge(width = 0.9))  +
  geom_text(stat = "count",
            aes(label = ..count..),
            position = position_dodge(width = 0.9),
            vjust = -0.5,  
            size = 3.5    
  ) +
  labs(
    x = "Staging at 6 Months",
    y = "Number of Patients",
    fill = "Surival Status",
    title = "Survival Status: Staging at 6 Months"
  ) +
  theme_minimal()
ggsave("survival_status.pdf", width = 6, height = 4)


# # leave patient 22 out, because missing 
# clinical_data_filtered <- clinical_data[clinical_data$PatID != "Pat_ID_22", ]
# 
# 
# # Calculate the survival days
# # Censor with last contact if there is no death date 
# clinical_data_filtered <- clinical_data_filtered %>% 
#   mutate(death_date_new = ifelse(is.na(death_date), last_contact, death_date))
# 
# clinical_data_filtered$death_date_new <- as_date(clinical_data_filtered$death_date_new)
# 
# # calculate time (as in days)
# clinical_data_filtered$time <- difftime(clinical_data_filtered$death_date_new, clinical_data_filtered$Date_first_application)



# Add event status ( 0 = censored (alive), and 1 = died) 
# censored <- c(1, 5, 13, 25, 29, 33, 43, 44, 45, 47)
# 
# test <- c()
# for(i in censored){
#  test <- c(test, (paste("PatID_", i, sep = "")))
# }

censored_ids <- c("PatID_1", "PatID_5", "PatID_13", "PatID_25", "PatID_29", "PatID_33", "PatID_43", "PatID_44", "PatID_45", "PatID_47")
clinical_data <- clinical_data %>%
  mutate(event_status = if_else(patient_id %in% censored_ids, 0, 1))


# Compare overall survival vs gender
fit <- survfit(Surv(OS_11_24, event_status) ~ Sex, data = clinical_data)
print(fit)

# Summary of survival curves
summary(fit)

# Access to the sort summary table
summary(fit)$table

# Survival information in dataframe 
# d <- data.frame(time = fit$time,
#                 n.risk = fit$n.risk,
#                 n.event = fit$n.event,
#                 n.censor = fit$n.censor,
#                 surv = fit$surv,
#                 upper = fit$upper,
#                 lower = fit$lower
# )
# head(d)


# Change color, linetype by strata, risk.table color by strata
surv_plot_os <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata",
           linetype = "strata", 
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"),
           title = "Overall survival for HCC Patients [Gender] ")
combined_plot <- arrange_ggsurvplots(list(surv_plot), print = FALSE)

ggsave("OS_KP_curve_gender.pdf", width = 20 , height = 8, plot = combined_plot)



# Compare overall survival vs staging status 
fit <- survfit(Surv(OS_11_24, event_status) ~ Staging_6months, data = clinical_data)
print(fit)


# Change color, linetype by strata, risk.table color by strata
surv_plot_os <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata",
           linetype = "strata", 
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"),
           title = "Overall survival for HCC Patients [Responder vs Non-Responder] ")


# Compare progression-free survival vs staging status 
fit <- survfit(Surv(PFS_11_24, event_status) ~ Staging_6months, data = clinical_data)
print(fit)

# Change color, linetype by strata, risk.table color by strata
surv_plot_pfs <- ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata",
           linetype = "strata", 
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           palette = c("#E7B800", "#2E9FDF"),
           title = "Progression-free survival for HCC Patients [Responder vs Non-Responder] ")

combined_plot <- arrange_ggsurvplots(list(surv_plot_os, surv_plot_pfs), print = FALSE)
ggsave("OS_PFS.pdf", width = 15 , height = 8, plot = combined_plot)


# Cox regression overall survival 
res.cox <- coxph(Surv(OS_11_24, event_status) ~ Age + Sex + viral + BMI_kg_m2 + Staging_6months, data =  clinical_data)
summary(res.cox)

forest_os <- ggforest(
  res.cox,
  main = "Hazard ratio with OS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.7,
  refLabel = "reference",
  noDigits = 2
)



# Cox regression progression-free survival
res.cox <- coxph(Surv(PFS_11_24, event_status) ~ Age + Sex + viral + BMI_kg_m2 + Staging_6months, data =  clinical_data)
summary(res.cox)

forest_pfs <- ggforest(
  res.cox,
  main = "Hazard ratio with PFS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.7,
  refLabel = "reference",
  noDigits = 2
)


ggsave("OS_PFS_forest.pdf", width = 15 , height = 8, plot = grid.arrange(forest_os, forest_pfs, ncol = 2))