# ===============================
# NHANES Prediabetes Subtype Analysis (Fixed)
# ===============================

# ---- Load packages ----
library(dplyr)
library(survival)
library(survminer)
library(broom)
library(ggplot2)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(xtable)
library(tableone)

# ---- Set working directory and read dataset ----
setwd("/Users/avitaneja/Documents/emory/extracirriculars/research")
mydata <- readRDS("rth01_selected df.RDS")

# ---- Identify prediabetes ----
no_diabetes <- (mydata$dm_doc_told == 2 | is.na(mydata$dm_age))
prediabetes <- no_diabetes & (
  (mydata$fasting_glucose >= 100 & mydata$fasting_glucose <= 125) |
  (mydata$glycohemoglobin >= 5.6 & mydata$glycohemoglobin <= 6.4)
)
mydata$prediabetes <- prediabetes

cat("Prediabetes cases:", sum(prediabetes, na.rm = TRUE), "\n")

# ---- Define cluster centroids ----
cluster_means <- data.frame(
  row.names = c("1_Overweight_normotensive", "2_Dysglycemia", "3_Overweight_hypertensive",
                "4_Obese_early_onset", "5_Obese_moderate_IR", "6_Older_severe_IR", "7_Obese_severe_IR"),
  age = c(59.64709118, 59.78944708, 54.89094371, 44.08253801, 59.49204648, 65.64438561, 55.36587772),
  glycohemoglobin = c(5.503029412, 6.204511533, 5.63236755, 5.76603968, 5.753035033, 5.766283716, 5.821352478),
  bmi = c(24.98430676, 29.06433963, 26.27414524, 33.06441474, 33.17694144, 29.98786514, 32.18790268),
  homa2b = c(74.63511765, 92.78880597, 79.0343819, 120.8730482, 113.3552203, 252.0642857, 159.3466883),
  homa2ir = c(1.062294118, 1.627645862, 1.262444812, 2.032620459, 1.940062435, 4.853146853, 3.481333951),
  egfr = c(74.51230988, 72.67762125, 72.81202589, 107.5359492, 80.17968411, 71.51735816, 68.67514723),
  sbp = c(134.7091678, 126.4644328, 145.4261477, 115.1221316, 136.6411015, 123.0394516, 145.8778439),
  dbp = c(85.86461382, 78.06768554, 98.86771362, 76.87933071, 88.07331892, 70.14724069, 96.93581813),
  ldl = c(113.5063959, 126.7215984, 147.7487464, 121.9864682, 122.3806313, 106.6826828, 146.5486526),
  hdl = c(65.66297647, 38.74463704, 32.6306239, 47.26272481, 51.53480923, 51.8545989, 24.81419592)
)

# ---- Variables to scale ----
vars <- c("age", "glycohemoglobin", "bmi", "homa2b", "homa2ir", "egfr", "sbp", "dbp", "ldl", "hdl")

# ---- Scale data ----
data_scaled <- scale(mydata[, vars], center = TRUE, scale = TRUE)
cluster_scaled <- scale(cluster_means[, vars],
                        center = attr(data_scaled, "scaled:center"),
                        scale = attr(data_scaled, "scaled:scale"))

# ---- Assign closest cluster safely ----
assign_cluster <- function(row) {
  dists <- apply(cluster_scaled, 1, function(c) sqrt(sum((row - c)^2)))
  names(which.min(dists))
}

mydata$closest_cluster <- NA  # initialize
non_missing <- complete.cases(mydata[, vars])
mydata$closest_cluster[non_missing] <- apply(data_scaled[non_missing, ], 1, assign_cluster)

# ---- Map clusters to subphenotypes ----
cluster_to_subtype <- c(
  "1_Overweight_normotensive" = "MARD",
  "2_Dysglycemia" = "SIDD",
  "3_Overweight_hypertensive" = "MARD",
  "4_Obese_early_onset" = "SIDD",
  "5_Obese_moderate_IR" = "MOD",
  "6_Older_severe_IR" = "SIRD",
  "7_Obese_severe_IR" = "SIRD"
)

mydata$subphenotype <- factor(
  cluster_to_subtype[mydata$closest_cluster],
  levels = c("MARD", "SIDD", "MOD", "SIRD")
)

# ---- Clean data for survival ----
mydata_clean <- mydata[!is.na(mydata$mortstat) &
                         !is.na(mydata$permth_int) &
                         !is.na(mydata$subphenotype), ]

# ---- Cox proportional hazards model ----
cox_adj <- coxph(Surv(permth_int, mortstat) ~ subphenotype + age + gender,
                 data = mydata_clean)
hr_df <- tidy(cox_adj, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(term = gsub("subphenotype", "", term)) %>%
  rename(Subtype = term, HR = estimate,
         Lower_95_CI = conf.low, Upper_95_CI = conf.high, p_value = p.value)

write.csv(hr_df, "Table3_HazardRatios.csv", row.names = FALSE)
hr_df
# ---- Figure 1: Flowchart ----
n_total <- nrow(mydata)
n_pred <- sum(mydata$prediabetes, na.rm = TRUE)
n_clean <- nrow(mydata_clean)

flow <- grViz(paste0("
digraph flowchart {
  node [shape=box, style=filled, fillcolor=lightgrey];
  A [label='Total NHANES participants (N=", n_total, ")'];
  B [label='Prediabetes identified (N=", n_pred, ")'];
  C [label='Complete survival data (N=", n_clean, ")'];
  A -> B -> C;
}
"))
rsvg_png(charToRaw(export_svg(flow)), "Figure1_flowchart.png")

# ---- Figure 2: Cumulative incidence ----
fit <- survfit(Surv(permth_int / 12, mortstat) ~ subphenotype, data = mydata_clean)
ggsurvplot(fit, data = mydata_clean, fun = "event",
           risk.table = TRUE, xlab = "Years", ylab = "Cumulative Incidence",
           legend.title = "Subtype",
           palette = c("red", "green3", "deepskyblue", "purple"))
ggsave("Figure2_CumulativeIncidence.png", width = 7, height = 5)

# ---- Table 1: Descriptive statistics ----
vars_numeric <- vars
vars_categorical <- c("gender", "mortstat")

NHANEStable <- tableone::CreateTableOne(vars = vars_numeric,
                                        strata = "subphenotype",
                                        data = mydata_clean,
                                        factorVars = vars_categorical)
print(NHANEStable, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)

# ---- Table 2: Death summary by year ----
mydata_clean$time_years <- mydata_clean$permth_exm / 12
mydata_clean$year_interval <- floor(mydata_clean$time_years)
table2 <- aggregate(mortstat ~ year_interval + subphenotype, mydata_clean, sum)
total_counts <- aggregate(mortstat ~ year_interval + subphenotype, mydata_clean, length)
table2$N_total <- total_counts$mortstat
table2$percent_dead <- round((table2$mortstat / table2$N_total) * 100, 2)
write.csv(table2, "Table2_DeathsByYear.csv", row.names = FALSE)
table2

# ---- Save workspace ----
saveRDS(list(mydata_clean = mydata_clean, hr_df = hr_df, table2 = table2),
        "NHANES_prediabetes_results.RDS")

cat("✅ Analysis complete — all tables and figures saved.\n")
