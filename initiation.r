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
library(missRanger)
library(tableone)
library(xtable)

# ---- Set working directory and load data ----
setwd("/Users/avitaneja/Documents/emory/extracirriculars/research")
mydata <- readRDS("rth01_selected df.RDS")


# ---- Filter to prediabetes cases ----
eligible_people <- mydata %>%
  filter(
    dm_doc_told != 1 &    # exclude those told they have diabetes
    (fasting_glucose < 126 | is.na(fasting_glucose)) &
    (glycohemoglobin < 6.5 | is.na(glycohemoglobin)) &
    (
      (fasting_glucose >= 100 & fasting_glucose <= 125) |
      (glycohemoglobin >= 5.6 & glycohemoglobin <= 6.4)
    )
  )

  View(eligible_people)
eligible_people %>% summarise(count = n())
write.csv(eligible_people, "eligible_people.csv", row.names = FALSE)
# ---- Impute missing data ----
missing_summary <- sapply(eligible_people, function(x) sum(is.na(x)))
missing_percent <- (missing_summary / nrow(eligible_people)) * 100
missing_percent
eligible_people_clean <- eligible_people[, missing_percent < 70]
set.seed(42)
imputed_eligible_people <- missRanger(
  data = eligible_people_clean,
  pmm.k = 3,             # predictive mean matching (more realistic imputations)
  num.trees = 100,       # number of trees in each random forest
  maxiter = 10,          # iterations until convergence
  verbose = TRUE
)
write.csv(imputed_eligible_people, "eligible_people_imputed.csv", row.names = FALSE)
# ---- Define pooled mean & SD ----
pooled_mean <- c(
  age = 56.06254821,
  glycohemoglobin = 5.765865914,
  bmi = 29.63809744,
  homa2b = 117.6386563,
  homa2ir = 2.091170981,
  egfr = 79.98366775,
  sbp = 132.0693503,
  dbp = 85.2937305,
  ldl = 126.8683825,
  hdl = 45.22730234
)

pooled_sd <- c(
  age = 10.49421806,
  glycohemoglobin = 0.444204112,
  bmi = 6.187436173,
  homa2b = 63.87755264,
  homa2ir = 1.473068288,
  egfr = 21.03660163,
  sbp = 22.1010336,
  dbp = 15.94766959,
  ldl = 36.89972999,
  hdl = 19.10282818
)

# ---- Define cluster means ----
cluster_means <- data.frame(
  row.names = c(
    "1_Overweight_normotensive", "2_Dysglycemia",
    "3_Overweight_hypertensive", "4_Obese_early_onset",
    "5_Obese_moderate_IR", "6_Older_severe_IR",
    "7_Obese_severe_IR"
  ),
  age = c(59.64709118, 59.78944708, 54.89094371, 44.08253801,
          59.49204648, 65.64438561, 55.36587772),
  glycohemoglobin = c(5.503029412, 6.204511533, 5.63236755,
                      5.76603968, 5.753035033, 5.766283716, 5.821352478),
  bmi = c(24.98430676, 29.06433963, 26.27414524,
          33.06441474, 33.17694144, 29.98786514, 32.18790268),
  homa2b = c(74.63511765, 92.78880597, 79.0343819,
             120.8730482, 113.3552203, 252.0642857, 159.3466883),
  homa2ir = c(1.062294118, 1.627645862, 1.262444812,
              2.032620459, 1.940062435, 4.853146853, 3.481333951),
  egfr = c(74.51230988, 72.67762125, 72.81202589,
           107.5359492, 80.17968411, 71.51735816, 68.67514723),
  sbp = c(134.7091678, 126.4644328, 145.4261477,
          115.1221316, 136.6411015, 123.0394516, 145.8778439),
  dbp = c(85.86461382, 78.06768554, 98.86771362,
          76.87933071, 88.07331892, 70.14724069, 96.93581813),
  ldl = c(113.5063959, 126.7215984, 147.7487464,
          121.9864682, 122.3806313, 106.6826828, 146.5486526),
  hdl = c(65.66297647, 38.74463704, 32.6306239,
          47.26272481, 51.53480923, 51.8545989, 24.81419592)
)

# ---- Variables to scale ----
vars <- c("age", "glycohemoglobin", "bmi", "homa2b", "homa2ir",
          "egfr", "sbp", "dbp", "ldl", "hdl")

# ---- Scale NHANES and cluster centroids using pooled mean & SD ----
data_scaled <- scale(imputed_eligible_people[, vars], center = pooled_mean, scale = pooled_sd)
cluster_scaled <- scale(cluster_means[, vars], center = pooled_mean, scale = pooled_sd)

# ---- Assign closest cluster ----
assign_cluster <- function(row) {
  dists <- apply(cluster_scaled, 1, function(c) sqrt(sum((row - c)^2)))
  names(which.min(dists))
}

imputed_eligible_people[rowSums(is.na(eligible_people)) > 0, ]
imputed_eligible_people$closest_cluster <- NA
imputed_eligible_people$closest_cluster <- apply(data_scaled, 1, assign_cluster)

View(imputed_eligible_people %>% select(closest_cluster, all_of(vars)))

# ---- Cox model comparing clusters to cluster 1 ----
imputed_eligible_people$closest_cluster <- factor(imputed_eligible_people$closest_cluster,
                                       levels = rownames(cluster_means))
cox_adj <- coxph(Surv(permth_int, mortstat) ~ closest_cluster + age + gender,
                 data = imputed_eligible_people)
hr_df <- tidy(cox_adj, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(term = gsub("closest_cluster", "", term)) %>%
  rename(Cluster = term, HR = estimate,
         Lower_95_CI = conf.low, Upper_95_CI = conf.high, p_value = p.value)
write.csv(hr_df, "Table3_HazardRatios.csv", row.names = FALSE)
hr_df

# ---- Figure 1: Flowchart ----
n_total <- nrow(mydata)
n_pred <- nrow(imputed_eligible_people)
flow <- grViz(paste0("
digraph flowchart {
  node [shape=box, style=filled, fillcolor=lightgrey];
  A [label='Total NHANES participants (N=", n_total, ")'];
  B [label='Prediabetes identified (N=", n_pred, ")'];
  A -> B;
}"))
rsvg_png(charToRaw(export_svg(flow)), "Figure1_flowchart.png")

# ---- Figure 2: Cumulative incidence ----
fit <- survfit(Surv(permth_int / 12, mortstat) ~ closest_cluster, data = imputed_eligible_people)
ggsurvplot(fit, data = imputed_eligible_people, fun = "event",
           risk.table = TRUE, xlab = "Years", ylab = "Cumulative Incidence",
           legend.title = "Cluster",
           palette = c("red", "green3", "blue", "purple", "orange", "brown", "cyan"))
ggsave("Figure2_CumulativeIncidence.png", width = 7, height = 5)

# ---- Table 1: Descriptive statistics ----
vars_categorical <- c("gender", "mortstat")
NHANEStable <- tableone::CreateTableOne(
  vars = vars,
  strata = "closest_cluster",
  data = imputed_eligible_people,
  factorVars = vars_categorical,
)
print(NHANEStable, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE, 
nonnormal = c("glycohemoglobin", "bmi", "homa2b", "homa2ir", "sbp", "ldl", "hdl")
  ) 

# ---- Table 2: Death summary year ----
View(imputed_eligible_people)
imputed_eligible_people$time_years <- imputed_eligible_people$permth_exm / 12
imputed_eligible_people$year_interval <- floor(imputed_eligible_people$time_years)
table2 <- aggregate(mortstat ~ year_interval + closest_cluster, data = imputed_eligible_people, sum)
total_counts <- aggregate(mortstat ~ year_interval + closest_cluster, data = imputed_eligible_people, length)
table2$N_total <- total_counts$mortstat
table2$percent_dead <- round((table2$mortstat / table2$N_total) * 100, 2)
table2
latex_table <- xtable(table2[, c("year_interval", "closest_cluster", "mortstat", "N_total", "percent_dead")])
print(latex_table, include.rownames = FALSE)
# ---- Save results ----
saveRDS(list(imputed_eligible_people = imputed_eligible_people, hr_df = hr_df, table2 = table2),
        "NHANES_prediabetes_results.RDS")

cat("✅ Analysis complete — all tables and figures saved.\n")



