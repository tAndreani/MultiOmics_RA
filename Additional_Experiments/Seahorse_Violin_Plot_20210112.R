# Introduction ------------------------------------------------------------

# Project: CIA (09/2020 - 12/2020)
# Experiment: Seahorse assay
# This R script re-makes some plots done on 20210110

# Load packages -----------------------------------------------------------

library(tidyverse)
library(data.table)
cat("\f") # Clear console

# Add OCR:ECAR ratio ------------------------------------------------------

# Import calculated data
df_Cal <- fread("./Seahorse/data/Seahorse_Combined_All_20210110.csv", quote = "")
View(df_Cal)

# Calculate OCR:ECAR ratio (After glucose, before Oligomycin)
df_Cal2 <- df_Cal %>%
  mutate(OCR_ECAR_Ratio = Avg_p456/AvgE_p456)
View(df_Cal2)

# Save data
fwrite(df_Cal2, file = "./Seahorse/data/Seahorse_Calculation_All_20210112.csv")

# Import saved data -------------------------------------------------------

df_Cal2 <- fread("./Seahorse/data/Seahorse_Calculation_All_20210112.csv", quote = "")
View(df_Cal2)
glimpse(df_Cal2)

# Add score data to Seahorse calculation ----------------------------------

# Import score data
scoreData <- fread("./Phenotype/data/CIA_Ctrl_Score_Sum_Sacrifice_20210118.csv", quote = "")
glimpse(scoreData)

# Merge score data to Seahorse calculation
df_Cal2 <- df_Cal2 %>% left_join(scoreData, by = c("mouseID", "Treatment", "absolutDay"))
glimpse(df_Cal2)

sum(is.na(df_Cal2$Score)) # No NA values in scores

# Save Seahorse data with scores
write.csv(df_Cal2, "./Phenotype/data/Seahorse_Calculation_with_Score_20210210.csv")

# Rank levels
df_Cal2$Treatment <- factor(df_Cal2$Treatment, levels = c("Ctrl", "CIA"))

df_Cal2$absolutDay <- as.factor(df_Cal2$absolutDay)

# Violin plot -------------------------------------------------------------

# Load saved Seahorse data with scores
df_Cal2 <- fread("./Phenotype/data/Seahorse_Calculation_with_Score_20210210.csv")
View(df_Cal2)
glimpse(df_Cal2)

# Change to long format
df_long <- df_Cal2 %>%
  select(absolutDay, Treatment, mouseID, Score, BR, ATP, MR, Glc) %>%
  pivot_longer(
    cols = 5:8,
    names_to = "Parameter",
    values_to = "Value",
    values_drop_na = TRUE
  )
View(df_long)

# Violin plot
df_long$Parameter <- factor(df_long$Parameter, levels = c("Glc", "BR", "ATP", "MR"))

labelNames <- c(
  "Glc" = "Glycolysis (mpH/min)",
  "BR" = "Basal respiration (pmol/min)",
  "ATP" = "ATP-coupled respiration (pmol/min)",
  "MR" = "Maximal respiration (pmol/min)"
)

# Metabolic calculation results
p_base <- df_long %>%
  filter(!absolutDay %in% c(-7, 23, 42)) %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  ggplot(aes(x = absolutDay, y = Value, fill = Treatment)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = .4, scale = "width", draw_quantiles = 0.5, size = 1, color = "#999999") +
  geom_point(aes(color = Score), position = position_jitterdodge(dodge.width = 0.9, seed = 2019), size = 1.5)

p_cal <- p_base +
  labs(x = "Day after immunization", y = NULL) +
  scale_fill_brewer(palette = "Dark2") +
  scale_colour_viridis_c(limits = c(0,12), breaks = c(0, 3, 6, 9, 12)) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  facet_wrap(vars(Parameter), labeller = as_labeller(labelNames), scales = "free")

ggsave("Calculation_Label_With_Score_20210118.png", plot = p_cal, width = 18, height = 16, units = "cm")

# OCR:ECAR ratio
p_ratio <- df_Cal2 %>%
  filter(!absolutDay %in% c(-7, 23, 42)) %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  ggplot(aes(x = absolutDay, y = OCR_ECAR_Ratio, fill = Treatment)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = .5, scale = "width", draw_quantiles = 0.5, size = 1, color = "#999999") +
  geom_point(aes(color = Score), position = position_jitterdodge(dodge.width = 0.9, seed = 2019), size = 1.5) +
  scale_colour_viridis_c(limits = c(0,12), breaks = c(0, 3, 6, 9, 12)) + 
  labs(x = "Day after immunization", y = "OCR:ECAR ratio") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw()

ggsave("OCR_ECAR_Ratio_Label_With_Score_20210118.png", plot = p_ratio, width = 12, height = 8, units = "cm")

# Plot for department meeting ---------------------------------------------

matchTable <- fread("./Phenotype/data/matchTable_absolutDay_Timepoint.csv")
matchTable

df_long$Treatment <- factor(df_long$Treatment, levels = c("Ctrl", "CIA"))

p_base <- df_long %>%
  left_join(matchTable[, 1:2], by = "absolutDay") %>%
  filter(Parameter %in% c("Glc", "BR"), absolutDay != -7) %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  mutate(Timepoint = factor(Timepoint, levels = c(
    "2 days", "2 weeks", "3 weeks",
    "4 weeks", "6 weeks", "7 weeks", "8 weeks", "10 weeks"
  ))) %>%
  ggplot(aes(x = Timepoint, y = Value, fill = Treatment)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = .4, scale = "width", draw_quantiles = 0.5, size = 1, color = "#999999") +
  geom_point(aes(color = Score), position = position_jitterdodge(dodge.width = 0.9, seed = 2019), size = 2.5, alpha = .8)

labelNames <- c("Glc" = "Glycolysis (mpH/min)", "BR" = "Basal respiration (pmol/min)")

p_cal <- p_base +
  labs(x = NULL, y = NULL) +
  scale_fill_brewer(palette = "Dark2") +
  scale_colour_viridis_c(limits = c(0,12), breaks = c(0, 3, 6, 9, 12)) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  facet_grid(Parameter ~ ., labeller = as_labeller(labelNames), scales = "free", switch = "y")

ggsave("./Seahorse/plot/Calculation_Label_With_Score_Dept_Meeting_20210903.png", plot = p_cal, width = 8, height = 6)

# Plot for TA review: Only early timepoints -------------------------------

# Metabolic calculation results
p_base <- df_long %>%
  filter(Parameter %in% c("Glc", "BR")) %>%
  filter(absolutDay %in% c(2, 14, 30)) %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  ggplot(aes(x = absolutDay, y = Value, fill = Treatment)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = .4, scale = "width", draw_quantiles = 0.5, size = 1, color = "#999999") +
  geom_point(aes(color = Score), position = position_jitterdodge(dodge.width = 0.9, seed = 2019), size = 1.5)

labelNames <- c("Glc" = "Glycolysis (mpH/min)", "BR" = "Basal respiration (pmol/min)")

p_cal <- p_base +
  labs(x = "Day after immunization", y = NULL) +
  scale_fill_brewer(palette = "Dark2") +
  scale_colour_viridis_c(limits = c(0,12), breaks = c(0, 3, 6, 9, 12)) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold")) +
  facet_grid(Parameter ~ ., labeller = as_labeller(labelNames), scales = "free", switch = "y")

ggsave("./Seahorse/plot/Calculation_Label_With_Score_TA_Review_2panels_20210511.png", plot = p_cal, width = 12, height = 16, units = "cm")

# Unpaired two-samples Mann-Whitney test ----------------------------------

# Load saved Seahorse data with scores
df_Cal2 <- fread("./Phenotype/data/Seahorse_Calculation_with_Score_20210210.csv")

# Change to long format
df_long <- df_Cal2 %>%
  filter(!mouseID %in% c(177, 115, 144, 200, 45, 79, 95, 117, 132, 155)) %>%
  select(absolutDay, Treatment, mouseID, BR, ATP, MR, Glc, OCR_ECAR_Ratio) %>%
  pivot_longer(
    cols = 4:8,
    names_to = "Parameter",
    values_to = "Value",
    values_drop_na = TRUE
  )
View(df_long)

# Define the timepoints to calculate p values
listTime <- c(2, 14, 30, 49, 56, 70)

# Create a variable for data frame
pVal <- NULL

for (comp in unique(df_long$Parameter)) {
  
  # Subset for parameter
  df_comp <- df_long %>%
    filter(Parameter == comp)
  
  for (i in listTime) {
    df_subset <- df_comp %>%
      filter(absolutDay == i)
    
    # U test
    res <- wilcox.test(Value ~ Treatment, data = df_subset, exact = FALSE)
    
    print(paste0(comp, " at Timepoint: ", i, " p-val: ", res$p.value))
    
    # Put into a data frame
    pVal <- rbind(pVal, data.frame(comp, i, res$p.value))
  }
}

View(pVal)

# Add stars to data frame pVal

pVal <- pVal %>%
  mutate(Star = case_when(
    res.p.value > 0.01 & res.p.value < 0.05 ~ "*",
    res.p.value > 0.001 & res.p.value < 0.01 ~ "**",
    res.p.value < 0.001 ~ "***",
    TRUE ~ "Not Sig."
  ))
  
pVal

# Save p values
fwrite(pVal, "./Phenotype/data/Seahorse_pValue_20210210.csv")
